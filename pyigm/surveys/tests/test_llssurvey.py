# Module to run tests on initializing AbsSurvey

# TEST_UNICODE_LITERALS

import numpy as np
import glob, os, imp, pdb
import pytest

from pyigm.surveys.llssurvey import LLSSurvey
from ..lls_literature import load_lls_lit

remote_data = pytest.mark.remote_data

lt_path = imp.find_module('linetools')[1]

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_read_hdlls_dr1_simple():
    hdlls = LLSSurvey.load_HDLLS(load_sys=False)
    assert hdlls.nsys == 157

@remote_data
def test_read_hdlls_dr1():   # This might actually be local now..
    hdlls = LLSSurvey.load_HDLLS()
    assert hdlls.nsys == 157

    CII_clms = hdlls.ions((6,2))
    gdCII = np.where(CII_clms['flag_N']>0)[0]
    assert len(gdCII) == 103


def test_dat_list():
    """JXP format :: Likely to be Deprecated
    """
    # LLS Survey
    if os.getenv('LLSTREE') is None:
        assert True
        return
    # Load
    lls = LLSSurvey.from_flist('Lists/lls_metals.lst', tree=os.getenv('LLSTREE'))
    # tests
    np.testing.assert_allclose(lls.NHI[0], 19.25)
    assert lls.nsys == 165


def test_sdss():
    """ Test SDSS DR7 -- This is very slow..
    """
    # All
    sdss_dr7_all = LLSSurvey.load_SDSS_DR7(sample='all')
    assert sdss_dr7_all.nsys == 1935
    # Stat
    sdss_dr7_stat = LLSSurvey.load_SDSS_DR7()
    assert len(sdss_dr7_stat.NHI) == 254


def test_hst():
    """ Test HST surveys
    """
    # ACS
    acs = LLSSurvey.load_HST_ACS()
    assert acs.nsys == 34
    assert len(acs.sightlines) == 18
    # WFC3
    wfc3 = LLSSurvey.load_HST_WFC3()
    assert wfc3.nsys == 91
    assert len(wfc3.sightlines) == 53
    # Combined
    HST_LLS = wfc3 + acs
    assert HST_LLS.nsys == 125
    assert len(HST_LLS.sightlines) == 71


def test_z3mage():
    """ Test z~3 MagE
    """
    # All
    z3mage = LLSSurvey.load_mage_z3()
    assert z3mage.nsys == 60
    assert len(z3mage.sightlines) == 105
    # Non-Color
    z3mage_NC = LLSSurvey.load_mage_z3(sample='non-color')
    assert z3mage_NC.nsys == 32
    assert len(z3mage_NC.sightlines) == 61


@remote_data
def test_literature():
    """ Literature list
    """
    lls_lit = load_lls_lit()
    assert lls_lit.nsys == 58
    assert lls_lit.ref == 'Zon04,Jen05,Tri05,Prx06a,Prx06b,Mei06,Mei07,Mei08,Nes08,Mei09,DZ09,Tum11,Kcz12,Bat12'


