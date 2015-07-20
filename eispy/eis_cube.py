# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
'''EIS spectral cube definitions'''

from __future__ import absolute_import

from astropy.io import fits
from sunpy.wcs.wcs import WCS
from sunpycube.cube.datacube import Cube
from sunpycube.cube import cube_utils as cu
from sunpy.wcs import wcs_util as wu
import numpy as np
from sunpycube.spectra.spectrum import Spectrum
from eispy.eis_spectral_cube import EISSpectralCube
import re

__all__ = ['EISCube']


def _clean(header):
    # TODO: find a way to identify cubes containing time
    """ Fixes non-standard or deprecated CTYPEn FITS keywords.

    Parameters
    ----------
    header : astropy.io.fits.Header
        The header to be cleaned.
    """
    header['ctype1'] = 'HPLN-TAN'  # Helioprojective longitude, TAN projection
    header['ctype2'] = 'HPLT-TAN'  # Helioprojective latitude, TAN projection
    header['ctype3'] = 'WAVE   '  # Wavelength axis, default (TAB) projection
    header['naxis'] = 3
    return header


class EISCube(Cube):
    '''EIS Cube subclass.

    References
    ----------
    For an overview of the mission
    http://solarb.mssl.ucl.ac.uk/SolarB/
    '''
    def __init__(self,
                 data, wcs, window=1, dataHeader=None, primaryHeader=None):
        '''
        Constructor function.

        Parameters
        ----------
        data: numpy ndarray
            The cube containing the data
        wcs: sunpy.wcs.wcs.WCS object
            The world coordinate system for the array.
        window: int
            The window this cube belongs to in the file. Used to fetch the
            correct metadata from the header
        dataHeader: astropy.io.fits.Header object
            The header for the BINTableHDU section of the FITS file
        primaryHeader: astropy.io.fits.Header object
            The main header for the whole file.
        '''
        header = _dictionarize_header(dataHeader, primaryHeader, window)
        Cube.__init__(self, data.T, wcs, meta=header)
        # Data is transposed here because EIS orders (y, lambda) by x or time,
        # not (y, x) by lambda.

    @classmethod
    def read(cls, filename, **kwargs):
        """ Reads in a given FITS file and returns a dictionary of new
        EISSpectralCubes. Additional parameters are given to fits.open.

        Parameters
        ----------
        filename : string
            Complete location of the FITS file
        """
        hdulist = fits.open(name=filename, **kwargs)
        header = _clean(hdulist[0].header)
        # TODO: Make sure each cube has a correct wcs.
        w = WCS(header=header, naxis=3)
        wavelengths = [c.name for c in hdulist[1].columns if c.dim is not None]
        data = [hdulist[1].data[wav] for wav in wavelengths]
        cubes = [EISCube(data[i], w, i+1, dataHeader=hdulist[1].header,
                         primaryHeader=hdulist[0].header)
                 for i in range(len(data))]
        return dict(zip(wavelengths, cubes))

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        # TODO: This will have to be changed once other sources are added.
        return True

    def convert_to_spectral_cube(self):
        """
        Converts this cube into an EISSpectralCube. It will only work if the
        cube has exactly three dimensions and one of those is a spectral axis.
        """
        if self.data.ndim == 4:
            raise cu.CubeError(4, "Too many dimensions: Can only convert a " +
                               "3D cube. Slice the cube before converting")
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2, 'Spectral axis needed to create a spectrum')
        axis = 0 if self.axes_wcs.wcs.ctype[-1] == 'WAVE' else 1
        coordaxes = [1, 2] if axis == 0 else [0, 2]  # Non-spectral axes
        newwcs = wu.reindex_wcs(self.axes_wcs, np.array(coordaxes))
        time_or_x_size = self.data.shape[coordaxes[0]]
        y_size = self.data.shape[coordaxes[1]]
        spectra = np.empty((time_or_x_size, y_size), dtype=Spectrum)
        for i in range(time_or_x_size):
            for j in range(y_size):
                spectra[i][j] = self.slice_to_spectrum(i, j)
        return EISSpectralCube(spectra, newwcs, self.meta)


def _is_in_window(key, window):
    '''
    Checks if a given key forms part of the specified spectral window.

    Parameters
    ----------
    key: str
        The key to be validated
    window: int
        The desired window
    '''
    end = re.findall(r'\d+$', key)  # finds numbers at the end of the key
    if len(end) == 0:
        return True
    else:
        return window == int(end[0])


def _dictionarize_header(data_header, primary_header, window):
    '''
    Combines the given FITS primary header and the bintable header for a
    specified window into a dictionary.

    Parameters
    ----------
    data_header: astropy.io.fits.Header object, dict, or dict-like object.
        secondary header to be pruned for the specified window
    primary_header: astropy.io.fits.Header object, dict, or dict-like object.
        The main FITS file header
    window: int
        The window to be chosen out of the data header.
    '''
    ph = dict(primary_header)
    dh = {}
    for k in data_header:
        if _is_in_window(k, window):
            newkey = re.sub(r'\d+$', '', k)
            dh[newkey] = data_header[k]

    ph.update(dh)
    return ph
