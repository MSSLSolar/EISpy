# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101, E0611
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
from bs4 import BeautifulSoup
from astroscrappy import detect_cosmics
from astropy import units as u
from astropy import constants
import re
from scipy.interpolate import UnivariateSpline
from math import exp

__missing__ = -100
__darts__ = "http://darts.jaxa.jp/pub/ssw/solarb/eis/"
__eff_area_ver_a__ = '004'
__eff_area_ver_b__ = '005'
__eff_areas_a__ = {}
__eff_areas_b__ = {}
__CCD_gain__ = 6.93
__pix_memo__ = {}
__sensitivity_tau__ = 1894.0
# TODO: general cleanup, maybe split into different files?


def eis_prep(filename, **kwargs):
    # TODO: Write docstring
    # TODO: FITS output
    # TODO: verbose?
    data_and_errors, meta = _read_fits(filename, **kwargs)
    if kwargs.get('zeros', True):
        _remove_zeros_saturated(*data_and_errors)
    if kwargs.get('darkcur', True):
        _remove_dark_current(meta, *data_and_errors)
    _calibrate_pixels(meta, *data_and_errors, **kwargs)
    if kwargs.get('interp', True):
        _interpolate_missing_pixels(*data_and_errors)
    if kwargs.get('cosmics', True):
        _remove_cosmic_rays(*data_and_errors, **kwargs)
    if kwargs.get('sens', True):
        _correct_sensitivity(meta, *data_and_errors)
    _radiometric_calibration(meta, *data_and_errors, **kwargs)
    return data_and_errors, meta


# /===========================================================================\
# |                          Methods for main steps                           |
# \===========================================================================/
def _read_fits(filename, **kwargs):
    """
    Reads a FITS file and returns two things: a dictionary of wavelengths to a
    3-tuple of 3D ndarrays (data, error, index) and a dictionary containing all
    the l0 metadata found in the file. The error array is initailly set to 0,
    and has the same shape as the data array. Extra keyword arguments are
    passed on to the astropy fits reader.

    Parameters
    ----------
    filename: str
        Location of the file to be opened.
    """
    hdulist = fits.open(filename, **kwargs)
    header = dict(hdulist[1].header)
    waves = [c.name for c in hdulist[1].columns if c.dim is not None]
    data = {wav: np.array(hdulist[1].data[wav], dtype=u.Quantity)
            for wav in waves}
    data_with_errors = [(data[k], np.zeros(data[k].shape),
                         waves.index(k) + 1) for k in data]
    return data_with_errors, header


def _remove_zeros_saturated(*data_and_errors):
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
        err[zeros] = __missing__


def _remove_dark_current(meta, *data_and_errors):
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
        err[negatives] = __missing__


def _calibrate_pixels(meta, *data_and_errors, **kwargs):
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
        for x_slice in range(err.size[0]):
            err[x_slice][locations] = __missing__
            data[x_slice][locations] = 0


def _interpolate_missing_pixels(*data_and_errors):
    """
    Interpolates missing pixels, marking them as corrected if this is the case.
    Error calculation is different when a pixel has been corrected.
    """
    for data, err, _ in data_and_errors:
        missing = np.array(np.where(err == __missing__)).T
        ymax = err.shape[1]
        for x, y, z in missing:
            y_p, y_n, n_weight = _get_neighbors(y, ymax, err, x, z)
            if n_weight is not None:
                data[x, y, z] = y_n * n_weight + y_p * (1 - n_weight)
                err[x, y, z] = 0  # mark as corrected


def _remove_cosmic_rays(*data_and_errors, **kwargs):
    """
    Removes and corrects for cosmic ray damage in the measurements. This method
    uses astroscrappy, so refer to that documentation for fine-tuning the
    keyword arguments.

    Parameters
    ----------
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be corrected
    kwargs: dict-like, optional
        Extra arguments to be passed on to astroscrappy.
    """
    kwargs = _clean_kwargs(**kwargs)
    for data, err, _ in data_and_errors:
        slices = [detect_cosmics(data[i], inmask=(err[i] == __missing__),
                                 **kwargs) for i in range(data.shape[0])]
        data = np.array([ccd_slice[1] for ccd_slice in slices])


def _correct_sensitivity(meta, *data_and_errors):
    """
    Corrects for the loss of CCD sensitivity throughout time.

    Parameters
    ----------
    meta: dict
        observation metadata
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be corrected
    """
    launch = dt.datetime(2006, 9, 22, 21, 36, 0)  # Hinode launch date
    obs_start = dt.datetime.strptime(meta['DATE_OBS'],
                                     "%Y-%m-%dT%H:%M:%S.000")
    delta = launch - obs_start
    days = delta.days + delta.seconds / 86400
    factor = exp(days / __sensitivity_tau__)
    for data, _, _ in data_and_errors:
        data /= factor


def _radiometric_calibration(meta, *data_and_errors, **kwargs):
    """
    Performs the radiometric calculations and conversions from data numbers to
    units of spectral radiance.

    Parameters
    ----------
    meta: dict
        observation metadata
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be corrected
    """
    obs_start = dt.datetime.strptime(meta['DATE_OBS'],
                                     "%Y-%m-%dT%H:%M:%S.000")
    obs_end = dt.datetime.strptime(meta['DATE_END'],
                                   "%Y-%m-%dT%H:%M:%S.000")
    total_time = (obs_end - obs_start).total_seconds()
    for data, err, index in data_and_errors:
        wl_start = meta['TWMIN' + str(index)]
        wl_end = meta['TWMAX' + str(index)]
        wl_num = data.shape[2]
        wavelengths = np.linspace(wl_start, wl_end, wl_num)
        detector = meta['TWBND' + str(index)]
        slit = 2 if meta['SLIT_IND'] == 2 else 0  # 1" slit has index 0
        _conv_dn_to_number_of_photons(data, wavelengths)
        seconds_per_exposure = total_time / data.shape[0]
        data /= seconds_per_exposure
        data /= u.s
        _calculate_errors(data, err, meta)
        _conv_photon_rate_to_intensity(data, wavelengths, detector, slit)
        if kwargs.get('phot2int', True):
            _conv_phot_int_to_radiance(data, wavelengths)
        data = data.to(u.erg / ((u.cm**2) * u.Angstrom * u.s * u.sr))


# /===========================================================================\
# |                              Utility methods                              |
# \===========================================================================/

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
    flatarr.sort()
    low_value = np.median(flatarr[:0.02 * flatarr.shape[0]])
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
    url = __darts__ + 'data/cal/dp/dusty_pixels.sav'
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
    url = __darts__ + 'data/cal/' + pix_type + '/'
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
    url = __darts__ + 'data/cal/'
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
    if kwargs.get('calhp', True):
        hots = _get_pixel_map(date, 'hp', detector, y_window, x_window)
    else:
        hots = np.zeros(shape)
    if kwargs.get('calwp', True):
        warms = _get_pixel_map(date, 'wp', detector, y_window, x_window)
    else:
        warms = np.zeros(shape)
    if kwargs.get('caldp', True):
        dusties = _get_dusty_array(y_window, x_window)
    else:
        dusties = np.zeros(shape)
    return hots, warms, dusties


# =======================    Radio calibration utils    =======================
def _get_pixel_solid_angle(detector, slit):
    """
    Calculates and returns the solid angle for a pixel used in the radiometric
    calculations. Returns a dimensionless Astropy Quantity. This corresponds to
    omega-d in the EIS software note 02.
    """
    mirror_focus_length = 1938.68 * u.mm
    if slit == 1:
        width = 9.3 * u.um
        pix_factor = 1.067 if detector == 'A' else 1.087
    elif slit == 2:
        width = 19.2 * u.um
        pix_factor = 2.08 if detector == 'A' else 2.119
    return (width**2) / (mirror_focus_length**2 * pix_factor)


def _get_effective_areas(detector):
    """
    Returns the effective detector areas for EIS, in a dictionary from floats
    (in Angstroms) to Quantities.
    """
    areas_dic = __eff_areas_a__ if detector == 'A' else __eff_areas_b__
    if len(areas_dic) > 0:
        return areas_dic
    url = __darts__ + 'response/EIS_EffArea_' + detector + '.'
    url += __eff_area_ver_a__
    urlfile = urllib.urlopen(url)
    lines = [l.strip() for l in urlfile.readlines() if l[0] != '#']
    urlfile.close()
    areas = [re.findall(r'\d+\.\d+', l) for l in lines]
    areas = {float(a[0]): float(a[1]) for a in areas}
    areas_dic .update(areas)
    return areas_dic


def _get_eff_area_at_wl(detector, wavelengths):
    """
    Returns a spline interpolation of the effective area at the given lambdas.
    Wavelengths should be floats, in Angstroms.
    """
    dic = _get_effective_areas(detector)
    spl = UnivariateSpline(dic.keys(), dic.values())
    return spl(wavelengths) * u.cm**2


def _conv_dn_to_number_of_photons(data_numbers, wavelengths):
    """
    Returns the number of photons that correspond to a given data number at a
    given wavelength. The wavelength must be a float in Angstroms because those
    are the units the conversion expects (Sorry, no Quantities yet)
    """
    # TODO: Make this use Quantities.
    data_numbers *= __CCD_gain__
    electrons_per_photon = 12398.5 / (wavelengths * 3.65)  # 3.65 is the energy
    # of an electron-hole in Si, and 12398.5 is the conversion from A to eV.
    data_numbers /= electrons_per_photon
    return data_numbers


def _conv_photon_rate_to_intensity(photons_per_second, wavelengths, detector,
                                   slit):
    """
    Converts an array containing the rate of photon incidence into an array of
    photons.cm-2.s-1.sr-1. Wavelengths should be floats, in Angstroms.
    """
    pix_solid_angle = _get_pixel_solid_angle(detector, slit)
    areas = _get_eff_area_at_wl(detector, wavelengths)
    photons_per_second /= pix_solid_angle
    for wav in range(len(areas)):
        data_slice = photons_per_second[:, :, wav]
        data_slice /= areas[wav]
    photons_per_second /= areas[0].unit
    photons_per_second /= u.sr


def _get_radiance_factor(wavelengths):
    """
    Returns the radiance conversion factor to convert photon intensities to
    spectral radiance.
    """
    wavelengths *= u.Angstrom
    return constants.c * constants.h / (wavelengths**2)


def _conv_phot_int_to_radiance(photon_intensity, wavelengths):
    """
    Converts an array containing the photon intensity of a measurement into
    spectral radiance, in units of power.area-1.time-1.wavelength-1.sr-1
    """
    radiance_factors = _get_radiance_factor(wavelengths)
    for i in range(len(wavelengths)):
        data_slice = photon_intensity[:, :, i]
        data_slice *= radiance_factors[i]
    photon_intensity *= radiance_factors[0].unit


def _calculate_errors(data, err, index, meta):
    """
    Calculates the standard deviation of the data.
    """
    mask = err == 0
    err[mask] = data[mask]
    x_start = meta['TDETX' + str(index)]
    dark_current_err = (2.26 if x_start < 1074 else
                        2.29 if 1074 <= x_start <= 2098 else
                        2.37 if 2098 <= x_start <= 3222 else
                        2.24)
    err[mask] += dark_current_err**2
    err[mask] = np.sqrt(err[mask])


# ========================    Interpolation utils    =========================
def _get_neighbors(y, ymax, err, x, z):
    """
    Returns the locations of the pixels to use to interpolate a missing pixel
    and their weights.
    """
    y_p = y_n = 0
    for i in range(1, 4):
        if y_p is 0:
            y_p = i if y - i >= 0 and err[x, y - i, z] >= 0 else 0
        if y_n is 0:
            y_n = i if y + i < ymax and err[x, y + i, z] >= 0 else 0
    if y_p == 0:
        weight = 1.0 if y_n == 1 else None
    elif y_n == 0:
        weight = 0.0 if y_p == 1 else None
    elif y_p == 1:
        weight = 0.5 if y_n == 1 else 1/3 if y_n == 2 else 2/9
    elif y_n == 1:
        weight = 2/3 if y_p == 2 else 7/9
    elif y_p == 2:
        weight = 0.5 if y_p == 2 else None
    else:
        weight = None
    return y - y_p, y + y_n, weight


# =====================    Cosmic Ray Removal utils    =======================
def _clean_kwargs(**kwargs):
    """
    Returns a kwargs dictionary containing only the ones used by astroscrappy.
    We have to do this because astroscrappy doesn't properly sanitize its
    input
    """
    cosmic_args = ['sigclip', 'sigfrac', 'onjlim', 'pssl', 'gain', 'readnoise',
                   'satlevel', 'niter', 'sepmed', 'cleantype', 'fsmode',
                   'psfmodel', 'psfwhm', 'psfsize', 'psfk', 'psfbeta',
                   'verbose']
    cosmic_kwargs = {k: kwargs[k] for k in kwargs if k in cosmic_args}
    return cosmic_kwargs
