# -*- coding: utf-8 -*-
import pytest
import os
import datetime as dt
import time
import numpy as np
from eispy import eis_utils as eu


def test_get_dict_from_file():
    df = eu.get_dict_from_file(dt.date(2008, 5, 12))
    assert df['data'].shape == (42629,)
    assert df['time'].shape == (42629,)


@pytest.fixture(params=[1, 2])
def time_corr(request):
    slit = request.param
    location = os.getcwd() + "/eispy/tests/orbit_slit" + str(slit) + "_nointerp.txt"
    orbit_file = open(location)
    orbit_file.readline()
    dic = {dt_from_str(line): float(line[29:-1]) for line in orbit_file}
    orbit_file.close()
    return dic, slit == 2


def dt_from_str(line):
    timestamp = time.strptime(line[:19], "%Y-%m-%dT%H:%M:%S")
    micros = int(line[20:23]) * 1000
    delta = dt.timedelta(microseconds=micros)
    return dt.datetime.fromtimestamp(time.mktime(timestamp)) + delta


def test_get_hk_thermal_correction(time_corr):
    dic, slit2 = time_corr
    months = [(2007, 8), (2008, 5), (2008, 9), (2008, 10)]
    arrs = {m: {d: dic[d] for d in filter(lambda t: (t.year, t.month) == m, dic)} for m in months}
    for month in arrs:
        dates = np.array(arrs[month].keys())
        expecteds = np.array(arrs[month].values())
        actuals = eu.calc_hk_thermal_corrections(dates, slit2)
        assert np.allclose(expecteds, actuals)


def test_datetime_to_ssw_time():
    times = [dt.datetime(1979, 1, 1, 0, 0, 0),
             dt.datetime(1979, 1, 2, 0, 0, 0)]
    expecteds = [0, 86400]
    actuals = [eu.datetime_to_ssw_time(t) for t in times]
    assert np.allclose(expecteds, actuals)
