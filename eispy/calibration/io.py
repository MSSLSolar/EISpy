# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101
"""
This module handles file IO for eis_prep.
"""
import numpy as np
import warnings
import re
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
    _remove_unreadable_cards(hdulist[0].header)
    _remove_unreadable_cards(hdulist[1].header)
    header = dict(hdulist[1].header)
    windows = kwargs.get('windows', None)
    waves = [c.name for c in hdulist[1].columns if c.dim is not None]
    data = {wav: np.array(hdulist[1].data[wav], dtype=u.Quantity)
            for wav in waves if windows is None or wav in windows}
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
    win = main_header['NWIN']
    _remove_unreadable_cards(main_header)
    _remove_unreadable_cards(data_header)
    _update_header(main_header, **kwargs)
    _update_header(data_header, **kwargs)
    data_columns = [_get_column(hdulist[1], data, index)
                    for data, _, index in data_and_errors]
    err_columns = [_get_column(hdulist[1], err, index)
                   for _, err, index in data_and_errors]
    data_cdefs = fits.ColDefs(data_columns + list(hdulist[1].columns[win:]))
    err_cdefs = fits.ColDefs(err_columns)
    newdatahdu = fits.BinTableHDU.from_columns(data_cdefs, header=data_header)
    newerrhdu = fits.BinTableHDU.from_columns(err_cdefs, header=data_header)
    _remove_unreadable_cards(newdatahdu.header)
    _remove_unreadable_cards(newerrhdu.header)
    for data, error, index in data_and_errors:
        _clean_data_header(newdatahdu, hdulist[1], data, index)
        _clean_data_header(newerrhdu, hdulist[1], err, index)
    dataoutlist = fits.HDUList([hdulist[0], newdatahdu])
    erroutlist = fits.HDUList([hdulist[0], newerrhdu])
    clobber = kwargs.pop('clobber', True)
    dataoutlist.writeto(outdir + _filename(main_header, 'l1', **kwargs),
                        clobber=clobber, output_verify='ignore')
    erroutlist.writeto(outdir + _filename(main_header, 'er', **kwargs),
                       clobber=clobber, output_verify='ignore')


# =========================    FITS Output utils    ==========================
def _update_header(header, **kwargs):
    """
    Updates the HDU header to include relevant calibration data. This does not
    do data and data-header specific updates.
    """
    header['DATA_LEV'] = 1
    header['NAXIS'] = 3
    if 'windows' in kwargs:
        header['NWIN'] = len(kwargs['windows'])
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


def _clean_data_header(newhdu, oldhdu, data, oldindex):
    """
    Updates the min and max values in the new header.
    """
    newheader = newhdu.header
    oldheader = oldhdu.header
    oldnames = oldhdu.columns.names
    newindex = newhdu.columns.names.index(oldnames[oldindex - 1]) + 1
    newheader['TDMIN' + str(newindex)] = data.min().value
    newheader['TDMAX' + str(newindex)] = data.max().value
    params = ['TDETX', 'TDETXW', 'TDETY', 'TDELT', 'TWAVE', 'TWMIN', 'TWMAX',
              'TWBND', 'WINHDR', 'CCDRON']
    for p in params:
        newcard = p + str(newindex)
        oldcard = p + str(oldindex)
        newheader[newcard] = oldheader[oldcard]


def _filename(header, filetype, **kwargs):
    """
    Returns the filename of the new file
    """
    date = dt.datetime.strptime(header['DATE_OBS'], "%Y-%m-%dT%H:%M:%S.000")
    datestr = date.strftime("%Y%m%d_%H%M%S")
    reduced = "_reduced" if 'windows' in kwargs else ""
    filename = "eis_" + filetype + "_" + datestr + reduced + ".fits"
    return filename


def _remove_unreadable_cards(header):
    try:
        dict(header)
    except fits.VerifyError as err:
        card = re.findall(r'\([^)]*\)', err.message)[0][1:-1]
        warnings.warn("WARNING: Card " + card + " was removed because it had" +
                      " an unreadable ASCII value.", UserWarning)
        header[card] = 0
        _remove_unreadable_cards(header)
