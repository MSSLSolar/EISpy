# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
Class for the EIS spectral cube.
"""

from sunpycube.spectra.spectral_cube import SpectralCube
from sunpycube.cube import cube_utils as cu
from sunpy.util.util import savitzky_golay
import numpy as np
from eispy import eis_utils as eu
import datetime as dt
from astropy import units as u
import time

__all__ = ['EISSpectralCube']


class EISSpectralCube(SpectralCube):
    """
    An EISSpectralCube is a subclass of SpectralCube with added methods for
    handling and correcting EIS's characteristics.
    """

    def __init__(self, spectra, wcs, meta, errors=None):
        SpectralCube.__init__(self, spectra, wcs, meta)
        self.errors = errors

    def orbital_correction_old(self, line_guess=None, yrange=None, **kwargs):
        """
        Calculates the corrections for orbital and slit tilt deviation.
        Calculates the average over a given y-range which ideally should be
        over quiet sun. If this is not given, then the whole cube will be used.
        A line guess may be given, the default is Fe XII. Extra arguments are
        given to the individual spectra's gaussian fit function.
        After this method is called, the corrections must be applied using the
        apply_corrections method. This is so that corrections are only
        calculated once and can be applied to cubes holding different peak
        wavelengths.

        Parameters
        ----------
        line_guess: tuple of three floats
            A guess for a strong line present in the data. The temperature
            drift is independent of wavelength, so this should simply be an
            easily fitted, strong, clear line present in the data.
        yrange: tuple of 2 Quantities
            The y-range to get the average from. This range should avoid active
            regions that may disrupt the fit. The default is to use the entire
            cube.
        """
        # TODO: handle multiple exposures in a single window
        line_guess = (1000, 195.12, 0.1) if not line_guess else line_guess
        if not yrange:
            yrange = (0, self.spectra.shape[1])
        else:
            ymin = cu.convert_point(yrange[0].value, yrange[0].unit,
                                    self.wcs, 1)
            ymax = cu.convert_point(yrange[1].value, yrange[1].unit,
                                    self.wcs, 1)
            yrange = (ymin, ymax)
        # First, get the centers of the line we're correcting for, throughout
        # the whole cube. This is a 2D array. Force a recalculation and clip
        # the fitting window as well.
        x_range = kwargs.get('x_range', (194.9, 195.3))
        kwargs.update({'recalc': True, 'x_range': x_range})
        centers = self._param_array(1, line_guess, **kwargs)
        # Now get a 1-D array of the average intensities along the y-axis.
        # Ideally the range should be chosen so that this is quiet sun.
        averages = [np.average(arr[yrange[0]:yrange[1]]) for arr in centers]
        # Calculate the offset from the centroid of the emission line over the
        # given data range.
        corrections = line_guess[1] - averages
        # Remove noise by appplying a smoothing filter and getting rid of
        # the less important frequencies.
        window_size = int(corrections.shape[0] / 3)
        window_size += 1 if window_size % 2 == 0 else 0
        corrections = savitzky_golay(corrections, window_size, 3)
        # Add the slit tilt component to all
        slit_tilt_corr = self._get_slit_tilt()
        corr = [slit_tilt_corr + a for a in corrections]
        return np.array(corr)

    def orbital_corrections_new(self):
        """
        Calculates the orbital corrections using the new housekeeping method,
        downloading the data if the fils are not found in the sunpy data folder
        """
        slit = 2 if self.meta['SLIT_IND'] == 2 else 0  # 1" slit has index 0
        # sit_and_stare = self.meta['FMIR_SS'] == 0
        slit_tilt_corr = self._get_slit_tilt()
        times = self.get_exposure_times()
        thermal_corr = eu.calc_hk_thermal_corrections(times, slit == 2)
        doppler_corr = eu.calc_doppler_shift(times)
        time_corr = thermal_corr + doppler_corr
        time_corr -= eu.wavelength_to_ccd_pixel(195.12 * u.Angstrom)
        corrections = [slit_tilt_corr + t for t in time_corr]
        return np.array(corrections)

    def _get_slit_tilt(self):
        """
        Calculates the corrections necessary due to slit tilt.
        """
        slit = 2 if self.meta['SLIT_IND'] == 2 else 0  # 1" slit has index 0
        obs_start_ts = time.strptime(self.meta['DATE_OBS'][:19],
                                     "%Y-%m-%dT%H:%M:%S")
        obs_start = dt.datetime.fromtimestamp(time.mktime(obs_start_ts))
        y_window_start = self.meta['TWS']
        n_y_pixels = self.meta['YW']
        wavelength = 'LONG' if self.meta['TWAVE'] > 230 else 'SHORT'
        slit_tilt_corr = eu.calc_slit_tilt(y_window_start, n_y_pixels,
                                           obs_start, wavelength, slit)
        return slit_tilt_corr

    def get_exposure_times(self):
        """
        Returns the exposure times for this SpectralCube
        """
        n_exposures = self.meta['NEXP']
        obs_start_ts = time.strptime(self.meta['DATE_OBS'][:19],
                                     "%Y-%m-%dT%H:%M:%S")
        obs_end_ts = time.strptime(self.meta['DATE_END'][:19],
                                   "%Y-%m-%dT%H:%M:%S")
        obs_start = dt.datetime.fromtimestamp(time.mktime(obs_start_ts))
        obs_end = dt.datetime.fromtimestamp(time.mktime(obs_end_ts))
        tdelta = (obs_end - obs_start) / n_exposures
        times = [obs_start + delt for delt in tdelta * np.arange(n_exposures)]
        return times

    def apply_corrections(self, corrections):
        """
        Applies the given corrections over the entire array of spectra.

        Parameters
        ----------
        corrections: 2D array
            The offsets to be applied.
        """
        for arr, corr in zip(self.spectra.T, corrections):
            for spec in arr:
                spec.shift_axis(corr)
