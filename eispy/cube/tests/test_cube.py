import gzip
import urllib.request
import os
import pytest

from eispy.cube import read, EISObservationL2


@pytest.mark.remote_data
@pytest.fixture
def l2_file():
    fname = 'eis_l2_20061103_200206.fits'
    remote_dir = 'http://solar.ads.rl.ac.uk/MSSL-data/eis/level2/2006/11/03/'
    if not os.path.exists(fname) and not os.path.exists(f'{fname}.gz'):
        filename, headers = urllib.request.urlretrieve(f'{remote_dir}{fname}.gz', f'{fname}.gz')

    if not os.path.exists(fname):
        with gzip.open(f'{fname}.gz', 'rb') as f:
            with open(fname, 'wb') as g:
                g.write(f.read())
    return fname


def test_l2_observation(l2_file):
    obs = read(l2_file)
    assert isinstance(obs, EISObservationL2)
