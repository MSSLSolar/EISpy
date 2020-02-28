from sunpy.net import Fido, attrs
import pytest

from eispy.cube import read, EISObservationL2


@pytest.mark.remote_data
@pytest.fixture
def eis_file():
    ###########################################################################
    # Download an EIS level 0 file
    res = Fido.search(attrs.Instrument('EIS'),
                      attrs.Time('2006-11-03 20:02:00', '2006-11-03 20:03:00'))
    downloaded_files = Fido.fetch(res)
    return downloaded_files[0]


def test_observation(eis_file):
    obs = read(eis_file)
    assert isinstance(obs, EISObservation)
