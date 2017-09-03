# Module to run tests on initializing AbslineSystem

# TEST_UNICODE_LITERALS

import numpy as np
import copy
import pytest

import astropy.units as u
from astropy.cosmology import FlatLambdaCDM, LambdaCDM

from pyigm.utils import cosm_xz

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_dxdz():
    # Default cosmology (Vanilla)
    Xz = cosm_xz(3.)
    np.isclose(Xz, 7.68841320742732)
    dxdz = cosm_xz(3., flg_return=1)
    np.isclose(dxdz, 3.5847294011983)
    # Open
    copen = LambdaCDM(H0=70., Om0=0.3, Ode0=0.)
    dxdz = cosm_xz(3., cosmo=copen, flg_return=1)
    np.isclose(dxdz, 2.90092902756)

