# Module to run tests on continuum codes

## # TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest

from pyigm.continuum import quasar as pyicq

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


def test_telfer():
    #
    telfer = pyicq.get_telfer_spec(3.)
    np.testing.assert_allclose(telfer.flux[100].value, 2.297435281318983)


def test_telfer_and_igm():
    telfer = pyicq.get_telfer_spec(zqso=3., igm=True, nproc=4)
    # Test
    np.testing.assert_allclose(telfer.flux[100].value, 0.9383676133341446)


def test_wfc3_conti():
    wfc3, _ = pyicq.wfc3_continuum(wfc3_indx=0, zqso=2.)
    np.testing.assert_allclose(wfc3.flux[100].value, 18.405926715695372)

def test_full_sed():
    lognu, fnu = pyicq.full_nu_spectrum()
    np.testing.assert_allclose(lognu[100], 12.51)
    np.testing.assert_allclose(fnu[100], 89.12509381337513)
