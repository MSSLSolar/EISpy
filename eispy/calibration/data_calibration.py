# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
This module deals with the data-wise interpolation, calibration and correction
for eis_prep
"""
import numpy as np
import datetime as dt
import astroscrappy as asc
import urllib
from astropy import units as u
from astropy import constants
import re
from scipy.interpolate import UnivariateSpline
from math import exp
from eispy.calibration.constants import missing, darts, eff_area_ver_a, \
                                        eff_area_ver_b, eff_areas_a, \
                                        eff_areas_b, CCD_gain, sensitivity_tau

__all__ = ['interpolate_missing_pixels', 'remove_cosmic_rays',
           'correct_sensitivity', 'radiometric_calibration']


def interpolate_missing_pixels(*data_and_errors):
    """
    Interpolates missing pixels, marking them as corrected if this is the case.
    Error calculation is different when a pixel has been corrected.

    Parameters
    ----------
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be corrected
    """
    for data, err, _ in data_and_errors:
        missing_pix = np.array(np.where(err == missing)).T
        ymax = err.shape[1]
        for x, y, z in missing_pix:
            y_p, y_n, n_weight = _get_neighbors(y, ymax, err, x, z)
            if n_weight is not None:
                data[x, y, z] = y_n * n_weight + y_p * (1 - n_weight)
                err[x, y, z] = 0  # mark as corrected


def remove_cosmic_rays(*data_and_errors, **kwargs):
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
    newargs = _clean_kwargs(**kwargs)
    for data, err, _ in data_and_errors:
        slices = [asc.detect_cosmics(data[i], inmask=(err[i] == missing),
                                     **newargs) for i in range(data.shape[0])]
        data = np.array([ccd_slice[1] for ccd_slice in slices])


def correct_sensitivity(meta, *data_and_errors):
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
    factor = exp(days / sensitivity_tau)
    for data, _, _ in data_and_errors:
        data /= factor


def radiometric_calibration(meta, *data_and_errors, **kwargs):
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
        wavelengths = _get_wavelengths(meta, index, data.shape[2])
        detector = meta['TWBND' + str(index)]
        slit = 2 if meta['SLIT_IND'] == 2 else 0  # 1" slit has index 0
        _conv_dn_to_number_of_photons(data, wavelengths)
        seconds_per_exposure = total_time / data.shape[0]
        data /= seconds_per_exposure
        data /= u.s
        _calculate_errors(data, err, index, meta)
        err /= u.s
        _conv_photon_rate_to_intensity(data, err, wavelengths, detector, slit)
        if kwargs.get('phot2int', True):
            _conv_phot_int_to_radiance(data, err, wavelengths)
        data = data.to(u.erg / ((u.cm**2) * u.Angstrom * u.s * u.sr))


# /===========================================================================\
# |                              Utility methods                              |
# \===========================================================================/
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
    areas_dic = eff_areas_a if detector == 'A' else eff_areas_b
    if len(areas_dic) > 0:
        return areas_dic
    url = darts + 'response/EIS_EffArea_' + detector + '.'
    url += eff_area_ver_a if detector == 'A' else eff_area_ver_b
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
    data_numbers *= CCD_gain
    electrons_per_photon = 12398.5 / (wavelengths * 3.65)  # 3.65 is the energy
    # of an electron-hole in Si, and 12398.5 is the conversion from A to eV.
    data_numbers /= electrons_per_photon
    return data_numbers


def _conv_photon_rate_to_intensity(photons_per_second, err, wavelengths,
                                   detector, slit):
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
    photons_per_second /= (areas[0].unit * u.sr)
    err /= (areas[0].unit * u.sr)


def _get_radiance_factor(wavelengths):
    """
    Returns the radiance conversion factor to convert photon intensities to
    spectral radiance.
    """
    wavelengths *= u.Angstrom
    return constants.c * constants.h / (wavelengths**2)


def _conv_phot_int_to_radiance(photon_intensity, err, wavelengths):
    """
    Converts an array containing the photon intensity of a measurement into
    spectral radiance, in units of power.area-1.time-1.wavelength-1.sr-1
    """
    radiance_factors = _get_radiance_factor(wavelengths)
    for i in range(len(wavelengths)):
        data_slice = photon_intensity[:, :, i]
        err_slice = err[:, :, i]
        goodmask = err_slice != missing
        data_slice *= radiance_factors[i]
        err_slice *= radiance_factors[i]
        err_slice[goodmask] = missing * data_slice[0, 0].unit
    err *= radiance_factors[0].unit
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


def _get_wavelengths(meta, index, size):
    """
    Returns an array of floats with the wavelength values in Angstroms.
    """
    wl_start = meta['TWMIN' + str(index)]
    wl_end = meta['TWMAX' + str(index)]
    return np.linspace(wl_start, wl_end, size)


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
