# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
This module handles file IO for eis_prep.
"""
import numpy as np
from astropy.io import fits
from astropy import units as u
import datetime as dt
from eispy.version import version as eispy_version
from eispy.calibration.constants import sensitivity_tau

__all__ = ['read_fits', 'write_to_fits']


def read_fits(filename, **kwargs):
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
    memmap = kwargs.pop('memmap', True)
    hdulist = fits.open(filename, memmap=memmap, **kwargs)
    header = dict(hdulist[1].header)
    waves = [c.name for c in hdulist[1].columns if c.dim is not None]
    data = {wav: np.array(hdulist[1].data[wav], dtype=u.Quantity)
            for wav in waves}
    data_with_errors = [(data[k], np.zeros(data[k].shape),
                         waves.index(k) + 1) for k in data]
    return data_with_errors, header


def write_to_fits(outdir, filename_in, *data_and_errors, **kwargs):
    """
    Writes the given data and errors to a new level 1 and error FITS files,
    with updated header values.

    Parameters
    ----------
    outdir: str
        Directory to output the files
    filename_in: str
        Location of the level 0 FITS file
    data_and_errors: one or more 3-tuples of ndarrays
        tuples of the form (data, error, index) to be corrected
    kwargs: dict
        Keyword arguments, used to update the header.
    """
    hdulist = fits.open(filename_in, memmap=True)
    main_header = hdulist[0].header
    data_header = hdulist[1].header
    _update_header(main_header, **kwargs)
    _update_header(data_header, **kwargs)
    win = len(data_and_errors)
    data_columns = [_get_column(hdulist[1], data, index)
                    for data, _, index in data_and_errors]
    err_columns = [_get_column(hdulist[1], err, index)
                   for _, err, index in data_and_errors]
    data_cdefs = fits.ColDefs(data_columns + list(hdulist[1].columns[win:]))
    err_cdefs = fits.ColDefs(err_columns)
    newdatahdu = fits.BinTableHDU.from_columns(data_cdefs, header=data_header)
    newerrhdu = fits.BinTableHDU.from_columns(err_cdefs, header=data_header)
    for data, error, index in data_and_errors:
        _clean_data_header(newdatahdu.header, data, index)
        _clean_data_header(newerrhdu.header, err, index)
    dataoutlist = fits.HDUList([hdulist[0], newdatahdu])
    erroutlist = fits.HDUList([hdulist[0], newerrhdu])
    clobber = kwargs.pop('clobber', True)
    dataoutlist.writeto(outdir + _filename(main_header, 'l1'), clobber=clobber)
    erroutlist.writeto(outdir + _filename(main_header, 'er'), clobber=clobber)


# =========================    FITS Output utils    ==========================
def _update_header(header, **kwargs):
    """
    Updates the HDU header to include relevant calibration data. This does not
    do data and data-header specific updates.
    """
    header['DATA_LEV'] = 1
    header['NAXIS'] = 3
    now = dt.datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
    _delete_cards(header)
    header.insert('TELESCOP', ('DATE_RF1', now,
                               'Date and time of Level 1 reformat'))
    header.insert('TELESCOP', ('ORIG_RF1',
                               kwargs.get('institute', "Unknown institute"),
                               'Institution where Level 1 reformat was done'))
    header.insert('TELESCOP', ('VER_RF1', "EISpy version " + eispy_version,
                               'EISpy version number'))
    header.insert('BITC_VER', ('CAL_DC', kwargs.get('darkcur', True),
                               "Dark current calibration"))
    header.insert('BITC_VER', ('CAL_HP', kwargs.get('calhp', True),
                               "Hot pixel calibration"))
    header.insert('BITC_VER', ('CAL_DP', kwargs.get('caldp', True),
                               "Dusty pixel calibration"))
    header.insert('BITC_VER', ('CAL_WP', kwargs.get('calwp', True),
                               "Warm pixel calibration"))
    header.insert('BITC_VER', ('CAL_CR', kwargs.get('cosmics', True),
                               "Cosmic ray calibration"))
    header.insert('BITC_VER', ('CAL_ABS', True, "Radiometric calibration"))
    header.insert('BITC_VER', ('CAL_PHOT', not kwargs.get('phot2int', True),
                               "T if units in photons"))
    header.insert('BITC_VER', ('CAL_SENS', kwargs.get('sens', True),
                               "Sensitivity correction"))
    header.insert('BITC_VER', ('CAL_RETA', False, "Negative values retained"))
    header.insert('BITC_VER', ('CAL_INT', kwargs.get('interp', True),
                               "Data interpolation done"))
    header.insert('BITC_VER', ('TAU_SENS', sensitivity_tau,
                               "Sensitivity value used"))


def _delete_cards(header):
    cards = ['DATE_RF1', 'ORIG_RF1', 'VER_RF1', 'CAL_DC', 'CAL_HP', 'CAL_WP',
             'CAL_DP', 'CAL_CR', 'CAL_ABS', 'CAL_PHOT', 'CAL_SENS', 'CAL_RETA',
             'CAL_INT', 'TAU_SENS', 'CAL_FF', 'CAL_WVL']
    for c in cards:
        header.pop(c, None)


def _get_column(hdu, data, index):
    """
    Updates the table and header values of the data HDU, and returns a new FITS
    Column object
    """
    name = hdu.columns[index - 1].name
    col = fits.Column(name=name, array=data, dim=str(data.shape[-1:0:-1]),
                      format=str(data.shape[1] * data.shape[2]) + 'E',
                      unit=data[0, 0, 0].unit.to_string('cds'))
    return col


def _clean_data_header(header, data, index):
    """
    Updates the min and max values in the new header.
    """
    header['TDMIN' + str(index)] = data.min().value
    header['TDMAX' + str(index)] = data.max().value


def _filename(header, filetype):
    """
    Returns the filename of the new file
    """
    date = dt.datetime.strptime(header['DATE_OBS'], "%Y-%m-%dT%H:%M:%S.000")
    datestr = date.strftime("%Y%m%d_%H%M%S")
    filename = "eis_" + filetype + "_" + datestr + ".fits"
    return filename
