# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101, C0330
"""
Utilities used in EIS calculations, corrections and fits.
"""

from scipy.io import readsav
from scipy.interpolate import interp1d
import numpy as np
import datetime as dt
import warnings
import sunpy
import urllib
from math import sqrt
import astropy.constants as const
from astropy import units as u

__housekeeping_memo__ = {}
__dispersion_long__ = 0.022332  # Constant value for long wavelength
__dispersion_short__ = 0.022317  # Constant value for short wavelength
__disp_sq_long__ = -1.329e-8  # Quadratic term, long wavelength
__disp_sq_short__ = -1.268e-8  # Quadratic term, short wavelength


def get_dict_from_file(date, prefix="eis3"):
    """
    Reads an IDL .sav file containing EIS housekeeping data and returns its
    contents as a python dictionary. For speed, if the file has already been
    read, it may return the contents from a hidden memo. If the file is not
    found in the location specified it will attempt to download it once and
    save the file in the location originally specified.

    Parameters
    ----------
    date: date or datetime object
        Date of the observation required. If the file is present in the sunpy
        data directory, it will be read from there, or downloaded to that
        location if it isn't.
    prefix: str
        file prefix (eis3 for thermal correction, fpp1 for doppler shift)
    """
    key = date.year * 100 + date.month
    download_dir = sunpy.config.get("downloads", "download_dir")
    download_dir += "/EIS_corrections/" + prefix + "_" + str(key) + ".sav"
    if key in __housekeeping_memo__:
        file_dict = __housekeeping_memo__[key]
    else:
        try:
            file_dict = readsav(download_dir, python_dict=True)
        except IOError:
            url = "http://sdc.uio.no/eis_wave_corr_hk_data/eis3_" + \
                   str(key) + ".sav"
            urllib.urlretrieve(url, filename=download_dir)
            file_dict = readsav(download_dir, python_dict=True)
            warnings.warn("File was not found, so it was downloaded and " +
                          "placed at the given location", UserWarning)
        __housekeeping_memo__.update({key: file_dict})
    return file_dict


def get_hk_temperatures(time, _pos=None):
    """
    Given a housekeeping filename and a time, returns the array of temperature
    correction values for that time. If the time is out of range for that file
    the closest available time will be used (i.e. the very first or very last)

    Parameters
    ----------

    time: datetime object
        The date and time of the observation. Should be present even if it is
        overridden because it is used to read the appropriate file
    _pos: int
        The index of the desired time in the file's time field. This overrides
        the time argument, and should be used internally only.
    """
    file_dict = get_dict_from_file(time)
    if _pos is None:
        timestamp = datetime_to_ssw_time(time)
        position = np.argmin(file_dict['time'] - timestamp)
    else:
        position = _pos
    pos_before = position - 5 if position > 5 else 0
    times = file_dict['time'].shape[0]
    pos_after = position + 5 if position + 5 < times else times - 1
    temps = np.zeros(34)
    # Excuse the magic numbers, these are the temperatures we're interested in,
    # as defined in the IDL file eis_sts3_temp.pro. I have no idea why these
    # and not other values, but never mind...
    main_temps = [1, 5, 7, 10, 11, 13, 14, 15, 16, 17, 21, 22, 23, 24, 26, 28]
    aux_temps = [1, 10, 11, 13, 14, 22, 23, 24, 26]
    temps[:16] = [file_dict['data'][position].temp[i] for i in main_temps]
    temps[16:25] = [file_dict['data'][pos_before].temp[i] for i in aux_temps]
    temps[25:] = [file_dict['data'][pos_after].temp[i] for i in aux_temps]
    return temps


def correct_pixel(temps, time=None, slit2=False):
    """
    Calculates the orbital correction for a single pixel

    Parameters
    ----------
    temps: numpy ndarray
        The housekeeping temperature values at the given time
    time: datetime object
        The time at which the observation took place. If no time is given, then
        current system time is assumed (note this will most likely be wrong)
    slit2: boolean
        If set to true, the correction will assume the observation was done
        using the 2" slit
    """
    if time is None:
        time = dt.datetime.now()
        warnings.warn("Time not set, assuming current time.", UserWarning)
    correction_arr, pixel_ref = _get_corr_parameters(time)
    slit2_offset = -8.2
    pixel_ref += slit2_offset if slit2 else 0
    return np.sum(correction_arr * (temps - 15.0) / 10.0) + pixel_ref


def _get_corr_parameters(time):
    """
    Returns the correct correction parameters for the given time. They are
    different because of three adjustments that have been made to the device.
    The coefficients were calculated from S. Kamio's neural network approach.

    Parameters
    ----------
    time: float
        The time of the observation, in SSW format
    """
    # Heater adjustment time
    adj1 = datetime_to_ssw_time(dt.datetime(2007, 11, 29, 00, 00, 00))
    # slit focus adjustment
    adj2 = datetime_to_ssw_time(dt.datetime(2008, 8, 24, 00, 00, 00))
    # grating focus adjustment
    adj3 = datetime_to_ssw_time(dt.datetime(2008, 10, 21, 8, 00, 00))
    if time < adj1:
        correction_arr = np.array([4.10562e-01, 2.51204e+00, -7.03979e-01,
                                   1.21183e+00, -1.46165e+00, -2.03801e+00,
                                   -5.09189e+00, -3.31613e+00, 2.28654e-01,
                                   3.72455e+00, 8.19741e-01, 1.17212e+00,
                                   3.19226e+00, 2.21462e+00, -2.76307e+00,
                                   -7.75230e+00, 2.27707e+00, 8.62746e-02,
                                   -3.87772e+00, 8.50736e-01, 2.50457e-01,
                                   -4.62109e+00, -1.49986e+00, -9.98911e-01,
                                   -5.24012e+00, -4.88090e+00, 8.41629e-01,
                                   1.53231e+00, -5.56888e+00, 5.46359e+00,
                                   5.00476e+00, 6.83911e+00, 2.10491e+00,
                                   6.89056e+00])
        pixel_ref = 1.34524e+3
    elif adj1 < time < adj2:
        correction_arr = np.array([-7.60169e+00, -1.46383e+00, 3.64224e+00,
                                   6.22838e+00, 1.02071e+00, -5.87856e+00,
                                   -7.07813e+00, -3.29145e+00, -2.68002e+00,
                                   6.44214e+00, -5.64250e+00, 9.41400e+00,
                                   1.02490e+01, 1.00514e+00, 1.54987e+01,
                                   -2.43897e+01, 6.93774e+00, 7.99804e+00,
                                   -4.24839e+00, 1.94191e+00, -4.11472e+00,
                                   2.67682e+00, 2.63193e+00, -1.58034e+00,
                                   -1.36976e+01, -1.78314e+00, -3.97698e+00,
                                   -5.86437e+00, 2.30465e+00, 1.23473e+01,
                                   -1.35947e+00, 1.85987e+00, 4.27904e+00,
                                   -4.35809e+00])
        pixel_ref = 1.34915e+3
    else:
        correction_arr = np.array([-9.69118e-01, 2.12159e+00, -2.99428e+00,
                                   2.61100e+00, 1.41035e+00, -9.76397e-01,
                                   -1.61651e+01, -9.94312e-01, 1.04603e+00,
                                   8.57033e-01, 2.07951e+00, 4.80522e+00,
                                   8.65133e+00, -2.37848e-02, 1.09901e+00,
                                   -5.51204e+00, 1.58325e+00, 1.97708e+00,
                                   -3.42620e+00, 1.76606e+00, 6.50817e+00,
                                   -7.19983e+00, -3.21551e+00, -6.81840e-01,
                                   -5.75801e+00, -1.08458e-01, -3.76701e+00,
                                   -3.05294e+00, -4.01884e+00, 1.00570e+01,
                                   4.61089e-01, 6.69429e+00, -6.84122e-01,
                                   4.38880e+00])
        pixel_ref = 1.34281e+3
        pixel_ref += 4.88 if adj2 < time < adj3 else 0
    return correction_arr, pixel_ref


def calc_hk_thermal_corrections(times, slit2=False):
    """
    For a given filename (or month in the format 'yyyymm') and times of
    measurements, calculate the corrections needed on each of those times,
    interpolating if necessary when the file does not contain those exact times

    Parameters
    ----------
    times: numpy array of datetime objects
        Times the observations occurred
    slit2: boolean
        Whether the observation was made using the 2" slit
    """
    # TODO: include good and bad data samples
    measurement_times = get_dict_from_file(times[0])['time']
    times = [datetime_to_ssw_time(t) for t in times]
    min_wanted_index = np.argmin(measurement_times - np.min(times))
    max_wanted_index = np.argmax(measurement_times - np.max(times))
    pixels = np.zeros(max_wanted_index - min_wanted_index + 1)
    for i in range(min_wanted_index, max_wanted_index):
        temperatures = get_hk_temperatures(times[0], _pos=i)
        pixels[min_wanted_index - i] = correct_pixel(temperatures,
                                                     measurement_times[i],
                                                     slit2)
    shifted_corrections_fun = interp1d(measurement_times, pixels)
    return shifted_corrections_fun(times)


def datetime_to_ssw_time(time):
    """
    Converts a datetime oject into SSW-format timestamp, which is the number of
    seconds elapsed since 1979-01-01 00:00:00.

    Parameters
    ----------
    time: datetime object
        The datetime object to convert.
    """
    epoch = dt.datetime(1979, 1, 1, 0, 0, 0)   # Solarsoft epoch
    delta = time - epoch
    return delta.total_seconds()


def calc_slit_tilt(y_window_start, n_y_pixels, date, wavelength, slit):
    """
    Calculates the slit tilt correction, returning it as an array to be applied
    to each pixel in the observation

    Parameters
    ----------
    y_window_start: int
        The pixel where the observation starts
    n_y_pixels: int
        The number of y pixels from the observation
    date:
        The date of the observation. This is used because the slit focus was
        adjusted on Aug 24, 2008.
    wavelength: 'SHORT' or 'LONG'
        The corrections depend on whether the observation was done using the
        short or long wavelength modes of the instrument.
    slit: 1 or 2
        Which of the two slits was used in the observation.
    """
    slit_focus_adjustment = dt.date(2008, 8, 24)
    coefficients = np.array(  # Before adjustment, short wl, 2" slit
                           [[-2.7607165427059641e+00, 6.4390116579832180e-03,
                             -2.7122949886142483e-06, 1.3035120912928136e-09],
                              # Before adjustment, short wl, 1" slit
                            [-5.5413445912854886e-01, 1.6348272018403623e-03,
                             -1.1681674813813158e-06, 1.7382970471312863e-10],
                              # Before adjustment, long wl, 2" slit
                            [-2.5786595570389181e+00, 6.1481799132252490e-03,
                             -3.1526317889607469e-06, 1.9165497210094085e-09],
                              # Before adjustment, long wl, 1" slit
                            [-3.8458103837911040e-01, 1.3331898117030505e-03,
                             -1.5399093968859745e-06, 7.8727203402240153e-10],
                              # After adjustment, short wl, 2" slit
                            [-3.0444030416965404e+00, 6.4056986231720847e-03,
                             -8.0170597073123649e-07, -1.8739881646780448e-10],
                              # After adjustment, short wl, 1" slit
                            [-8.5068355624974856e-01, 2.2398405213417405e-03,
                             -1.0665536954454296e-06, -1.2311442746502073e-10],
                              # After adjustment, long wl, 2" slit
                            [-2.4247171434108168e+00, 5.0832726508360793e-03,
                             -7.9727705770693547e-07, 2.3158597348832410e-10],
                              # After adjustment, long wl, 1" slit
                            [-2.0225535365170799e-01, 7.4854735225926561e-04,
                             -7.6258247316829397e-07, 1.4085716859395248e-10]])
    coef_index = 0
    coef_index += 4 if date > slit_focus_adjustment else 0
    coef_index += 2 if wavelength is 'LONG' else 0
    coef_index += 1 if slit == 1 else 0
    poly_coefs = coefficients[coef_index]
    y_pixels = np.arange(y_window_start, y_window_start + n_y_pixels)
    y_polyval = np.polyval(poly_coefs[::-1], y_pixels)
    dispersion_factor = 0.0223
    return y_polyval * dispersion_factor


def calc_dispersion(wavelength):
    """
    Calculates dispersion at a given wavelength.

    Parameters
    ----------
    wavelength: Astropy Quantity
        The wavelength at which to calculate dispersion
    """
    ccd_pix = wavelength_to_ccd_pixel(wavelength)
    if wavelength > 230 * u.Angstrom:  # Long lambda detector
        return __dispersion_long__ + __disp_sq_long__ * ccd_pix
    else:
        return __dispersion_short__ + __disp_sq_short__ * ccd_pix


def wavelength_to_ccd_pixel(wavelength):
    """
    Converts a wavelength into a pixel position on the CCD

    Parameters
    ----------
    wavelength: Astropy Quantity
        The wavelength to convert
    """
    if wavelength > 230 * u.Angstrom:
        refwvl = 199.9389
        displin = __dispersion_long__
        dispsq = __disp_sq_long__
    else:
        refwvl = 166.131
        displin = __dispersion_short__
        dispsq = __disp_sq_short__
    pixel = refwvl - (wavelength.to(u.Angstrom)).value
    pixel *= 4 * dispsq
    pixel = sqrt(displin**2 - pixel)
    pixel -= displin
    pixel /= 2
    pixel /= dispsq
    return pixel


def calc_doppler_shift(times):
    """
    Calculates the offset to apply to the measurements due to line-of-sight
    velocity. Searches for the data file in the default sunpy data directory
    and downloads it to that location if it isn't found.

    Parameters
    ----------
    times: array of datetime objects
        Times at which to calculate Doppler shift
    """
    fpp_dict = get_dict_from_file(times[0], prefix="fpp")
    ssw_times = [datetime_to_ssw_time(t) for t in times]
    # File has data as an array of 1-tuples for some reason...
    data = np.array([v[0] for v in fpp_dict['data']])
    fpp_times = np.array(fpp_dict['time'])
    data_at_times_wanted = interp1d(fpp_times, data)(ssw_times)
    dispersion = calc_dispersion(195.12 * u.Angstrom)  # FeXII line
    doppler_shift = data_at_times_wanted / const.c.value * 195.12 / dispersion
    return doppler_shift
