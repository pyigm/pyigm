# Module to run tests on initializing DLASurvey

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from pyigm.surveys.dlasurvey import DLASurvey

remote_data = pytest.mark.remote_data

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


def test_init():
    dlas = DLASurvey(ref='null')
    assert dlas.abs_type == 'DLA'


def test_sdss():
    # All
    sdss = DLASurvey.load_SDSS_DR5(sample='all')
    # Testing
    assert sdss.nsys == 1182
    # Stat
    sdss = DLASurvey.load_SDSS_DR5()
    assert len(sdss.NHI) == 737


def test_read_h100_nosys():
    h100 = DLASurvey.load_H100(load_sys=False)
    assert h100.nsys == 100


@remote_data
def test_read_h100():
    h100 = DLASurvey.load_H100()
    assert h100.nsys == 100

    SiII_clms = h100.ions((14, 2))
    gdSiII = np.where(SiII_clms['flag_N'] > 0)[0]
    assert len(gdSiII) == 98


def test_dat_list():
    """JXP format :: Likely to be Deprecated
    """
    if os.getenv('DLA') is None:
        assert True
        return
    # Load
    dlas = DLASurvey.neeleman13_tree()
    # tests
    assert dlas.nsys == 100



