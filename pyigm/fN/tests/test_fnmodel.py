# Module to run tests on fNModel

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
from astropy import units as u

from pyigm.fN.fnmodel import fNModel

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

'''
def test_load_kin():
    # Class
    cos_halos = COSHalos()
    cos_halos.load_mega(skip_ions=True)
    # Load kin
    cos_halos.load_abskin()
'''

def test_init():
    # Class
    fN_P13 = fNModel('Hspline', zmnx=(2.,5.))
    # Mtype
    assert fN_P13.mtype == 'Hspline'
    # Class
    fN_I14 = fNModel('Gamma')
    # Mtype
    assert fN_I14.mtype == 'Gamma'
    assert fN_I14.zmnx == (0., 10.)

def test_default():
    fN_default = fNModel.default_model()
    #
    assert fN_default.mtype == 'Hspline'
    assert fN_default.zmnx == (0.5, 3)

def test_lx():
    fN_default = fNModel.default_model()
    lX = fN_default.calculate_lox(2.4, 17.19+np.log10(2.), 23.)
    np.testing.assert_allclose(lX, 0.36298679339713974)

def test_teff():
    fN_default = fNModel.default_model()
    zval,teff_LL = fN_default.teff_ll(0.5, 2.45)
    #
    np.testing.assert_allclose(zval[0], 0.5)
    np.testing.assert_allclose(teff_LL[0], 1.8176161746504436)
