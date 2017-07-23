# Module to run tests on initializing AbslineSystem

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from pyigm.euvb.cuba import CUBA

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_init():
    # Class
    cuba = CUBA()
    # z values
    np.testing.assert_allclose(cuba.z[0:5], 
        [0.,  0.04912,  0.1006 ,  0.1547 ,  0.2114])

def test_zinterp():
    # 
    cuba = CUBA()
    # 
    jnu = cuba.zinterp_jnu(2.1)
    np.testing.assert_allclose(jnu[0:5].value, 
        [1.36534694e-19, 1.33934694e-19, 1.31393878e-19, 1.28893878e-19, 1.26493878e-19])
    assert(jnu.unit == u.erg/u.s/u.cm**2)



def test_phi():
    cuba = CUBA()
    phi = cuba.phi(2.1)
    #
    np.testing.assert_allclose(phi.value, 381837.59540310985)
    assert(phi.unit == 1/u.s/u.cm**2)
    #
    phi = cuba.phi(2.1, min_energy=50.*u.eV)
    np.testing.assert_allclose(phi.value, 4*8171.47189724)


def test_logU():
    cuba = CUBA()
    logU = cuba.logU(0.2)
    #
    np.testing.assert_allclose(logU, -6.13089862882)

