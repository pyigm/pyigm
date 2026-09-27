# Module to run tests on generating IGMSystem



import numpy as np
import pytest

from astropy.cosmology import Planck15

from ..analysis import dndx_rvir

def test_dndx_rvir():
    # Simple properties
    Lval, dNdX = dndx_rvir(cosmo=Planck15)
    np.testing.assert_allclose(dNdX[0], 15.515609, rtol=1e-5)

