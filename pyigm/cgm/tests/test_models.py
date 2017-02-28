# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
from numpy.random import rand

from astropy import units as u

from ..models import ModifiedNFW


def test_modified_NFW():
    # Init
    mNFW = ModifiedNFW()
    # xyz
    xyz = (1-2*rand(3, 100)) * 100
    # rho
    rho = mNFW.rho(xyz)
    assert rho.size == 100
    assert rho.unit == u.g/u.cm**3
    # nH
    nH = mNFW.nH(xyz)
    assert rho.size == 100
    xyz0 = [100., 0., 0.]
    nH0 = mNFW.nH(xyz0)
    assert np.isclose(nH0, 0.00025117200064614744)
    # ne
    ne = mNFW.ne(xyz)
    assert np.all(ne > nH)


