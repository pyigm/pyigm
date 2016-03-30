# Module to run tests on FNModel

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
from astropy import units as u

from pyigm.fN.fnmodel import FNModel


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_init_I14():
    # Class
    fN_I14 = FNModel('Gamma')
    # Mtype
    assert fN_I14.mtype == 'Gamma'
    assert fN_I14.zmnx == (0., 10.)


def test_default():
    fN_default = FNModel.default_model()
    #
    assert fN_default.mtype == 'Hspline'
    assert fN_default.zmnx == (0.5, 3)


def test_lx():
    fN_default = FNModel.default_model()
    lX = fN_default.calculate_lox(2.4, 17.19+np.log10(2.), 23.)
    np.testing.assert_allclose(lX, 0.36298679339713974)


def test_mfp():
    fN_default = FNModel.default_model()
    z=2.44
    mfp = fN_default.mfp(z)
    # Test
    assert mfp.unit == u.Mpc
    np.testing.assert_allclose(mfp.value, 257.10258545808983)


def test_rhoHI():
    fN_default = FNModel.default_model()
    z=2.44
    rho_HI = fN_default.calculate_rhoHI(z, (20.3, 22.))
    # Test
    assert rho_HI.unit == u.solMass/u.Mpc**3
    np.testing.assert_allclose(rho_HI.value, 83553755.13025708)



