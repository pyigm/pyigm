# Module to run tests on generating IGMSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from astropy.cosmology import Planck15

from ..analysis import dndx_rvir

def test_dndx_rvir():
    # Simple properties
    Lval, dNdX = dndx_rvir(cosmo=Planck15)
    np.testing.assert_allclose(dNdX[0], 15.515609, rtol=1e-5)

