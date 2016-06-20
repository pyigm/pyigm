# Module to run tests on initializing AbslineSystem

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import copy

import astropy.units as u

from pyigm.igm import reion

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_damping():
    # Single value
    tau0 = reion.igm_damping(7000.*u.AA, 5.9, 6.3)
    np.testing.assert_allclose(tau0, 99.)
    # Array
    wave = np.linspace(7000.,10000,1000)*u.AA
    tau = reion.igm_damping(wave, 5.9, 6.3)
    np.testing.assert_allclose(tau[-1], 0.0024535867434168235)


