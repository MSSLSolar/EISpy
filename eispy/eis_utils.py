# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
Utilities used in EIS calculations, corrections and fits.
"""

from scipy.io import readsav
import numpy as np
import datetime as dt


def get_hk_temperatures_from_file(filename, time):
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
    """
    file_dict = readsav(filename, python_dict=True)
    position = time_to_index(time, file_dict['time'])
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


def time_to_index(time, times_arr):
    epoch = dt.datetime(1979, 01, 01, 00, 00, 00)  # Solarsoft epoch
    dif = time - epoch
    timestamp = dif.total_seconds()
    closest = np.argmin(times_arr - timestamp)
    return closest
