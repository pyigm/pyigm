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
    sdss_stat = DLASurvey.load_SDSS_DR5()
    assert len(sdss_stat.NHI) == 737
    # Binned
    lX, lX_lo, lX_hi = sdss_stat.calculate_lox([2., 2.5, 3])
    assert np.isclose(lX[0], 0.04625038, atol=1e-5)
    fN, fN_lo, fN_hi = sdss_stat.calculate_fn([20.3, 20.5, 21., 21.5, 22.], [2, 2.5], log=True)
    assert fN.size == 4
    assert np.isclose(fN_lo[0], 0.0682087, atol=1e-5)


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

def test_read_HST():
    """ Neeleman+16
    """
    hst16 = DLASurvey.load_HST16()
    assert hst16.nsys == 4

def test_read_xq100():
    """ XQ-100 """
    xq100 = DLASurvey.load_XQ100(sample='stat')
    assert xq100.nsys == 36

def test_read_p03_g09():
    """ XQ-100 """
    p03 = DLASurvey.load_P03()
    assert p03.nsys == 105

    g09 = DLASurvey.load_G09()
    assert g09.nsys == 38

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



