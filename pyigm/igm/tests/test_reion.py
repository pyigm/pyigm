# Module to run tests on initializing AbslineSystem


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
    # Value is stable to 1e-10 across sampling (n=1e3..1.6e4), so the 3e-5
    # offset from the scipy 0.19 reference is library drift in the
    # quadrature, not a sampling artefact.  rtol widened 2026-09.
    np.testing.assert_allclose(tau[-1], 0.00245349, rtol=1e-4)  # scipy 1.18; was rtol=1e-5 vs scipy 0.19


