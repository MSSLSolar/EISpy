'''EIS spectral cube definitions'''
import numpy as np

from astropy.io import fits
from astropy.nddata import StdDevUncertainty as sdu
from astropy.wcs import WCS

# from eispy.eis_spectral_cube import EISSpectralCube
# from eispy.calibration.constants import missing

from ndcube import NDCube

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
    header['CTYPE1'] = 'HPLT-TAN'  # Helioprojective longitude, TAN projection
    header['CTYPE2'] = 'HPLN-TAN'  # Helioprojective latitude, TAN projection
    header['CTYPE3'] = 'WAVE   '  # Wavelength axis, default (TAB) projection
    header['NAXIS'] = 3
    return header


class EISCube(NDCube):
    '''
    EIS Cube subclass.

    References
    ----------
    For an overview of the mission
    http://solarb.mssl.ucl.ac.uk/SolarB/
    '''
    @staticmethod
    def read(filename, er_filename=None):
        """
        Reads in a given FITS file and returns a dictionary of new
        EISSpectralCubes.

        Parameters
        ----------
        filename: string
            Location of the FITS file
        er_filename: string
            Location of the error FITS file
        """
        hdulist = fits.open(name=filename)
        errlist = fits.open(er_filename) if er_filename else None
        header = _clean(hdulist[0].header)
        # TODO: Make sure each cube has a correct wcs.
        wavelengths = [c.name for c in hdulist[1].columns if c.dim is not None]
        data = [hdulist[1].data[wav] for wav in wavelengths]
        errs = [errlist[1].data[wav] if errlist is not None else None
                for wav in wavelengths]
        cubes = []
        for i in range(len(data)):
            window = i + 1
            header = _dictionarize_header(hdulist[1].header,
                                          hdulist[0].header,
                                          window)
            uncertainty = sdu(errs[i]) if errlist else None
            header['NAXIS1'], header['NAXIS2'], header['NAXIS3'] = data[i].shape
            wcs = WCS(header=header, naxis=3)
            data[i] = data[i].T

            cubes += [EISCube(data[i], wcs, uncertainty=uncertainty)]
        return dict(zip(wavelengths, cubes))

    '''
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        # TODO: This will have to be changed once other sources are added.
        return True
    '''

    '''
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
        time_or_x_size = self.data.shape[coordaxes[1]]
        y_size = self.data.shape[coordaxes[0]]
        spectra = np.empty((time_or_x_size, y_size), dtype=Spectrum)
        for i in range(time_or_x_size):
            for j in range(y_size):
                spectra[i][j] = self.slice_to_spectrum(j, i)
        return EISSpectralCube(spectra, newwcs, self.meta)
    '''


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
    for key in data_header:
        if _is_in_window(key, window):
            newkey = re.sub(r'\d+$', '', key)
            dh[newkey] = data_header[key]

    dh.update(ph)
    dh['CRPIX3'] = 1
    dh['CRVAL3'] = dh['TWAVE']
    dh.pop('COMMENT', '')
    dh.pop('NAXIS1')
    return dh
