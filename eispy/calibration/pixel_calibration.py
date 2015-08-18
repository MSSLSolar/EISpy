# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
This module deals with pixel-wise calibration for eis_prep - zero values, dark
current, hot, warm and dusty pixels.
"""
import datetime as dt
import numpy as np
from scipy.io import readsav
import locale
import heapq
from bs4 import BeautifulSoup
import urllib
from eispy.calibration.constants import missing, darts

__all__ = ['remove_zeros_saturated', 'remove_dark_current', 'calibrate_pixels']

__pix_memo__ = {}


def remove_zeros_saturated(*data_and_errors):
    """
    Finds pixels in the data where the data numbers are zero or saturated,
    sets them to zero and marks them as bad in the error array. Note that this
    method modifies arrays in-place, and does not create or return new arrays.

    Parameters
    ----------
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be stripped of invalid data
    """
    for data, err, _ in data_and_errors:
        zeros = data <= 0
        saturated = data >= 16383  # saturated pixels have a value of 16383.
        zeros[saturated] = True  # equivalent to |=
        data[zeros] = 0
        err[zeros] = missing


def remove_dark_current(meta, *data_and_errors):
    # TODO: support 40" slot
    # TODO: support dcfiles
    """
    Calculates and subtracts dark current and CCD pedestal from the data. If
    the retain keyword is set to True then values less than zero are kept,
    otherwise they are floored at zero and marked as missing.

    Parameters
    ----------
    meta: dict
        observation metadata
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be corrected
    """
    for data, err, idx in data_and_errors:
        ccd_xwidth = meta['TDETXW' + str(idx)]
        if ccd_xwidth == 1024:
            _remove_dark_current_full_ccd(data, meta, idx)
        else:
            _remove_dark_current_part_ccd(data)
        negatives = data <= 0
        data[negatives] = 0
        err[negatives] = missing


def calibrate_pixels(meta, *data_and_errors, **kwargs):
    """
    Fetches and marks as missing hot, warm and dusty pixels present in the
    observation. If there is no data available for the exact date of a
    particular observation, the closest available one is used.

    Parameters
    ----------
    meta: dict
        observation metadata
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be corrected
    """
    date = dt.datetime.strptime(meta['DATE_OBS'][:10], "%Y-%m-%d")
    y_window = (meta['YWS'], meta['YWS'] + meta['YW'])
    for data, err, index in data_and_errors:
        x_window = (meta['TDETX' + str(index)], (meta['TDETX' + str(index)] +
                                                 meta['TDETXW' + str(index)]))
        detector = meta['TWBND' + str(index)].lower()
        hots, warms, dusties = _cal_arrays(date, detector, y_window, x_window,
                                           data.shape[1:], **kwargs)
        locations = hots == 1
        locations[warms == 1] = True
        locations[dusties == 1] = True
        for x_slice in range(err.shape[0]):
            err[x_slice][locations] = missing
            data[x_slice][locations] = 0


# ==========================    Dark current utils    =========================
def _remove_dark_current_full_ccd(data, meta, window):
    """
    Remove the dark current and CCD pedestal from a data array that encompasses
    the entire CCD.

    Parameters
    ----------
    data: 3D numpy array
        The CCD data
    meta: dict
        observation metadata
    window: int
        data window from FITS file
    """
    ccd_x_start = meta['TDETX' + str(window)]
    if ccd_x_start == 1024:
        pixels = (944, 989)
    elif ccd_x_start == 3072:
        pixels = (926, 971)
    else:
        pixels = (39, 84)
    quiet_vals = data[:, :, pixels[0]:pixels[1]]
    data -= np.median(quiet_vals)


def _remove_dark_current_part_ccd(data):
    """
    Remove the dark current and CCD pedestal from a data array that takes up
    only part of the CCD

    Parameters
    ----------
    data: 3D numpy array
        The CCD data
    """
    flatarr = data.flatten()
    lowest = 0.01 * flatarr.shape[0]
    low_value = heapq.nsmallest(int(lowest), flatarr)[-1]
    data -= low_value


# =======================    Pixel calibration utils    =======================
def _download_calibration_data(date, pix_type, detector, top_bot, left_right):
    """
    Downloads the requested calibration data from the DARTS repository. If the
    required data is not present, it looks for the nearest one that fits the
    requirements.

    Parameters
    ----------
    date: datetime object
        Date of the observation.
    pix_type: 'dp', 'hp', 'wp'
        Whether to download data for dusty, hot or warm pixels.
    detector: 'a', 'b'
        Long- or Short-wave detector
    top_bot: 'top', 'bottom', 'middle', 'both'
        Y location of the observation on the CCD
    left_right: 'left', 'right', 'both'
        X location of the observation on the CCD (left, right, or both)
    """
    tb_tuple = ('top', 'bottom') if top_bot == 'both' else (top_bot,)
    lr_tuple = ('left', 'right') if left_right == 'both' else (left_right,)
    retfiles = {}
    for vert in tb_tuple:
        for horiz in lr_tuple:
            key = vert + horiz
            arr = _try_download_nearest_cal(date, pix_type, detector, vert,
                                            horiz)
            retfiles.update({key: arr})
    return retfiles


def _get_dusty_array(y_window, x_window):
    """
    Returns the sliced array of dusty pixels
    """
    url = darts + 'data/cal/dp/dusty_pixels.sav'
    http_down = urllib.urlretrieve(url)
    dusties = readsav(http_down[0]).dp_data
    return dusties[y_window[0]:y_window[1], x_window[0]: x_window[1]]


def _try_download_nearest_cal(date, pix_type, detector, top_bot, left_right):
    """
    Tries to download the requested calibration data, looking for up to one
    month before and after to do so.
    """
    key = _construct_hot_warm_pix_url(date, pix_type, detector, top_bot,
                                      left_right)
    if key in __pix_memo__:
        return __pix_memo__[key]
    dates = _get_cal_dates(pix_type)
    dates.sort(key=lambda d: d - date if d > date else date - d)
    # dates is now a sorted list of the dates closest to the specified date
    for cal_date in dates:
        url = _construct_hot_warm_pix_url(cal_date, pix_type, detector,
                                          top_bot, left_right)
        http_response = urllib.urlopen(url)
        http_response.close()
        if http_response.code == 200:  # HTTP OK
            http_down = urllib.urlretrieve(url)
            arr = readsav(http_down[0]).ccd_data
            __pix_memo__[key] = arr
            return arr


def _get_cal_dates(pix_type):
    """
    Retrieves the list of available dates for a given pixel type.
    """
    url = darts + 'data/cal/' + pix_type + '/'
    http_request = urllib.urlopen(url)
    soup = BeautifulSoup(http_request.read())
    http_request.close()
    links = soup.find_all('a')
    date_str = [link.get('href') for link in links]
    dates = [dt.datetime.strptime(date, '%Y-%m-%d/') for date in date_str[5:]]
    return dates  # This isn't a numpy array so we can sort by keys.


def _construct_hot_warm_pix_url(date, pix_type, detector, top_bot, left_right):
    """
    Constructs a DARTS URL to download hot or warm pixels given the relevant
    parameters.
    """
    url = darts + 'data/cal/'
    url += pix_type + '/'
    url += date.strftime("%Y-%m-%d") + '/'
    url += 'coords_' + detector + '_' + left_right + '_'
    loc = locale.getlocale()
    locale.setlocale(locale.LC_ALL, 'en_GB')
    datestr = date.strftime("%d%b%y").lower()
    locale.setlocale(locale.LC_ALL, loc[0])
    if pix_type == 'wp':
        url += top_bot + '_'
        url += datestr + '_100s.sav'
    else:
        url += datestr
        if top_bot != 'middle':
            url += '_' + top_bot
        url += '.sav'
    return url


def _calculate_detectors(date, y_window_start, n_y_pixels, x_start, x_width):
    """
    Calculates what area of the detector the observation is in.
    """
    if date <= dt.datetime(2008, 1, 18):
        top_bot = 'middle'
    else:
        top_bot = 'top' if y_window_start >= 512 else \
                  'bottom' if y_window_start + n_y_pixels <= 512 else \
                  'both'
    # XXX: Warning: there may be an off-by-one error here!
    left_right = 'right' if x_start - 50 >= 1024 else \
                 'left' if (x_start - 50 + x_width) < 1024 else \
                 'both'  # CCD starts at pixel 50
    return top_bot, left_right


def _get_pixel_map(date, pix_type, detector, y_window, x_window):
    """
    Returns the pixel calibration map for the specified pixel type and detector
    at the given date.
    """
    y_window_start = y_window[0]
    x_start = x_window[0] % 2048
    n_y_pixels = y_window[1] - y_window[0]
    x_width = x_window[1] - x_window[0]
    detector_areas = _calculate_detectors(date, y_window_start, n_y_pixels,
                                          x_start, x_width)
    arrays = _download_calibration_data(date, pix_type, detector,
                                        detector_areas[0], detector_areas[1])
    glued_array = np.zeros((1024, 2048))
    zero_arr = np.zeros((512, 1024))
    glued_array[0:512, 0:1024] = arrays.get('topleft', zero_arr)
    glued_array[512:1024, 1024:2048] = arrays.get('bottomright', zero_arr)
    glued_array[512:1024, 0:1024] = arrays.get('bottomleft', zero_arr)
    glued_array[0:512, 1024:2048] = arrays.get('topright', zero_arr)
    return glued_array[y_window_start:y_window[1],
                       x_start:(x_window[1] % 2048)]


def _cal_arrays(date, detector, y_window, x_window, shape, **kwargs):
    """
    Checks for keyword arguments and if some pixels are disabled, simply
    returns a zero array for them.
    """
    verbose = kwargs.get('verbose', False)
    if kwargs.get('calhp', True):
        if verbose:
            print "Fetching hot pixel map..."
        hots = _get_pixel_map(date, 'hp', detector, y_window, x_window)
    else:
        hots = np.zeros(shape)
    if kwargs.get('calwp', True):
        if verbose:
            print "Fetching warm pixel map..."
        warms = _get_pixel_map(date, 'wp', detector, y_window, x_window)
    else:
        warms = np.zeros(shape)
    if kwargs.get('caldp', True):
        if verbose:
            print "Fetching dusty pixel map..."
        dusties = _get_dusty_array(y_window, x_window)
    else:
        dusties = np.zeros(shape)
    return hots, warms, dusties
