# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
Utilities used in EIS calculations, corrections and fits.
"""

from scipy.io import readsav
from scipy.interpolate import interp1d
import numpy as np
import datetime as dt
import warnings
import os
import re
import urllib


__housekeeping_memo__ = {}


def get_dict_from_file(filename):
    """
    Reads an IDL .sav file containing EIS housekeeping data and returns its
    contents as a python dictionary. For speed, if the file has already been
    read, it may return the contents from a hidden memo. If the file is not
    found in the location specified it will attempt to download it once and
    save the file in the location originally specified.

    Parameters
    ----------
    filename: str
        Location of the file to be read. Note that as shorthand or simply to
        force a download (if the file hasn't been read already) one can use the
        shorthand "yyymm" (e.g. "200805") instead of the entire qualified
        location of the file.
    """
    key = int(re.findall(r"(\d{6})(\.sav)?$", filename)[0][0])
    if key in __housekeeping_memo__:
        file_dict = __housekeeping_memo__[key]
    else:
        filename = os.getcwd() + '/' + filename
        if not filename.endswith(".sav"):
                filename += ".sav"
        try:
            print filename
            file_dict = readsav(filename, python_dict=True)
        except IOError:
            url = "http://sdc.uio.no/eis_wave_corr_hk_data/eis3_" + \
                   str(key) + ".sav"
            urllib.urlretrieve(url, filename=filename)
            file_dict = readsav(filename, python_dict=True)
            warnings.warn("File was not found, so it was downloaded and " +
                          "placed at the given location", UserWarning)
        __housekeeping_memo__.update({key: file_dict})
    return file_dict


def get_hk_temperatures(filename, time, _pos=None):
    """
    Given a housekeeping filename and a time, returns the array of temperature
    correction values for that time. If the time is out of range for that file
    the closest available time will be used (i.e. the very first or very last)

    Parameters
    ----------
    filename: str
        The location of an IDL .sav file containing the housekeeping data for
        the particular month.
    time: datetime object
        The date and time of the observation
    _pos: int
        The index of the desired time in the file's time field. This overrides
        the time argument, and should be used internally only.
    """
    file_dict = get_dict_from_file(filename)
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


def calc_hk_orbital_corrections(filename, times, slit2=False):
    """
    For a given filename (or month in the format 'yyyymm') and times of
    measurements, calculate the corrections needed on each of those times,
    interpolating if necessary when the file does not contain those exact times

    Parameters
    ----------
    filename: str
        Location of the file to be used or date in the format 'yyyymm'
    times: numpy array of datetime objects
        Times the observations occurred
    slit2: boolean
        Whether the observation was made using the 2" slit
    """
    # TODO: include good and bad data samples
    measurement_times = get_dict_from_file(filename)['time']
    times = [datetime_to_ssw_time(t) for t in times]
    min_wanted_index = np.argmin(measurement_times - np.min(times))
    max_wanted_index = np.argmax(measurement_times - np.max(times))
    pixels = np.zeros(max_wanted_index - min_wanted_index + 1)
    for i in range(min_wanted_index, max_wanted_index):
        temperatures = get_hk_temperatures(filename, None, _pos=i)
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
