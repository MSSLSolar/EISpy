# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
Calibration and error calculation for EIS level 0 files.
This module calls several corrections, then attempts to interpolate missing or
damaged data, and calculates the 1-sigma errors of the good data.
"""

from __future__ import absolute_import
from astropy.io import fits
from scipy.io import readsav
import numpy as np
import datetime as dt
import locale
import urllib
__missing__ = -100
__darts__ = "http://darts.jaxa.jp/pub/ssw/solarb/eis/data/cal/"


def _read_fits(filename, **kwargs):
    """
    Reads a FITS file and returns two things: a dictionary of wavelengths to a
    2-tuple of 3D ndarrays (data, error) and a dictionary containing all the l0
    metadata found in the file. The error array is initailly set to 0, and has
    the same shape as the data array. Extra keyword arguments are passed on to
    the astropy fits reader.

    Parameters
    ----------
    filename: str
        Location of the file to be opened.
    """
    hdulist = fits.open(filename, **kwargs)
    header = dict(hdulist[1].header)
    wavelengths = [c.name for c in hdulist[1].columns if c.dim is not None]
    data = {wav: hdulist[1].data[wav] for wav in wavelengths}
    data_with_errors = {k: (data[k], np.zeros(data[k].shape)) for k in data}
    return data_with_errors, header


def _remove_zeros_saturated(*data_and_errors):
    """
    Finds pixels in the data where the data numbers are zero or saturated,
    sets them to zero and marks them as bad in the error array. Note that this
    method modifies arrays in-place, and does not create or return new arrays.

    Parameters
    ----------
    data_and_errors: one or more 2-tuples of ndarrays
        tuples of the form (data, error) to be stripped of invalid data
    """
    for data, err in data_and_errors:
        zeros = data <= 0
        saturated = data >= 16383  # saturated pixels have a value of 16383.
        zeros[saturated] = True  # equivalent to |=
        data[zeros] = 0
        err[zeros] = __missing__


def _remove_dark_current(meta, *data_and_errors, **kwargs):
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
    data_and_errors: one or more 2-tuples of ndarrays
        tuples of the form (data, error) to be corrected
    retain: bool
        If True, data less than zero will be retained.
    """
    window = 1
    retain = kwargs.pop('retain', False)
    for data, err in data_and_errors:
        ccd_xwidth = meta['TDETXW' + str(window)]
        if ccd_xwidth == 1024:
            _remove_dark_current_full_ccd(data, meta, window)
        else:
            _remove_dark_current_part_ccd(data)
        if not retain:
            negatives = data <= 0
            data[negatives] = 0
            err[negatives] = __missing__
        window += 1


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
    flatarr.sort()
    low_value = np.median(flatarr[:0.02 * flatarr.shape[0]])
    data -= low_value


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
    for tb in tb_tuple:
        for lr in lr_tuple:
            key = tb + lr
            arr = _try_download_nearest_cal(date, pix_type, detector, tb, lr)
            retfiles.update({key: arr})
    return retfiles


def _try_download_nearest_cal(date, pix_type, detector, top_bot, left_right):
    """
    Tries to download the requested calibration data, looking for up to one
    month before and after to do so.
    """
    day = dt.timedelta(1)
    for i in range(30):  # Search one month before and after
        datebefore = date - (day * i)
        url = _construct_hot_warm_pix_url(datebefore, pix_type, detector,
                                          top_bot, left_right)
        http_response = urllib.urlopen(url)
        http_response.close()
        if http_response.code != 200:
            dateafter = date + (day * i)
            url = _construct_hot_warm_pix_url(dateafter, pix_type, detector,
                                              top_bot, left_right)
            http_response = urllib.urlopen(url)
            http_response.close()
        if http_response.code == 200:
            http_down = urllib.urlretrieve(url)
            return readsav(http_down[0]).ccd_data
    raise UserWarning("No data found within one month of the specified date")
    return np.zeros((512, 1024))


def _construct_hot_warm_pix_url(date, pix_type, detector, top_bot, left_right):
    """
    Constructs a DARTS URL to download hot or warm pixels given the relevant
    parameters.
    """
    url = __darts__
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
        url += datestr + '_'
        url += top_bot + '.sav'
    return url
