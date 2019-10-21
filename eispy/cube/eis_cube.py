'''EIS spectral cube definitions'''
import numpy as np

from astropy.io import fits
from astropy.nddata import StdDevUncertainty as sdu
from astropy.time import Time
from astropy.wcs import WCS

# from eispy.eis_spectral_cube import EISSpectralCube
# from eispy.calibration.constants import missing

from ndcube import NDCube

import re

__all__ = ['read', 'EISCube', 'EISObservation', 'EISObservationL2']


def _clean(header):
    # TODO: find a way to identify cubes containing time
    """
    Fixes non-standard or deprecated CTYPEn FITS keywords.

    Parameters
    ----------
    header : astropy.io.fits.Header
        The header to be cleaned.
    """
    header['CTYPE1'] = 'HPLN-TAN'  # Helioprojective longitude, TAN projection
    header['CTYPE2'] = 'HPLT-TAN'  # Helioprojective latitude, TAN projection
    header['CTYPE3'] = 'WAVE   '  # Wavelength axis, default (TAB) projection
    header['NAXIS'] = 3
    header['DATE-OBS'] = header.pop('DATE_OBS')
    # Drop the non-numbered keys that are already stored in the numbered keys
    for key in ['CRPIX', 'CRVAL', 'CDELT', 'CUNIT', 'CTYPE', 'CROTA']:
        if key in header:
            header.pop(key)
    return header


def read(filename, er_filename=None):
    """
    Reads in a given .fits file.

    Parameters
    ----------
    filename: string
        Location of the FITS file.
    er_filename: string
        Location of the error FITS file.

    Returns
    -------
    EISObservation
    """
    hdulist = fits.open(name=filename)
    errlist = fits.open(er_filename) if er_filename else None
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

        if data_level == 2:
            vmap = data[...]
            wmap = data[...]

        header['NAXIS1'], header['NAXIS2'], header['NAXIS3'] = data[i].shape
        wcs = WCS(header=header, naxis=3)
        data[i] = data[i].T

        cubes += [EISCube(data[i], wcs, uncertainty=uncertainty, meta=header)]

    primary_header = _clean(hdulist[0].header)
    return EISObservation(wavelengths, cubes, primary_header)


class EISObservation:
    """
    A single EIS observation (that comes from a single .fits file).

    This class stores a series of `EISCubes`, one for each wavelength that was
    observed.

    Parameters
    ----------
    wavelengths : list of str
        List of wavelengths observed.
    cubes : list of EISCube
        List of data cubes for each wavelength.
    primary_header : dict
        Primary data header for this observation.
    """
    def __init__(self, wavelengths, cubes, primary_header):
        self._cubes = dict(zip(wavelengths, cubes))
        self._header = primary_header

    def __getitem__(self, wavelength):
        """
        Return the cube corresponding to *wavelength*.
        """
        return self.cubes[wavelength]

    @property
    def wavelengths(self):
        """
        List of wavelengths.
        """
        return list(self._cubes.keys())

    @property
    def cubes(self):
        """
        List of data cubes.
        """
        return self._cubes

    @property
    def obs_starttime(self):
        return Time(self._header['DATE-OBS'])


class EISObservationL2(EISObservation):
    """
    An EIS observation derived from an L2 data product.

    In addition to spectral intensities, also contains rough estimates of the
    velocity and width maps.

    Parameters
    ----------
    wavelengths : list of str
        List of wavelengths observed.
    cubes : list of EISCube
        List of data cubes for each wavelength.
    """
    def __init__(self, wavelengths, cubes, primary_header, vmaps, wmaps):
        super.__init__(wavelengths, cubes, primary_header)
        self._vmaps = dict(zip(wavelengths, vmaps))
        self._wmaps = dict(zip(wavelengths, wmaps))

    def velocity_map(wavelength):
        return self._vmaps[wavelength]

    def width_map(wavelength):
        return self._wmaps[wavelength]


class EISCube(NDCube):
    '''
    EIS Cube subclass.

    References
    ----------
    For an overview of the mission
    http://solarb.mssl.ucl.ac.uk/SolarB/
    '''
    pass


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

    return _clean(dh)
