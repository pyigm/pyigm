# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
from numpy.random import rand

from astropy import units as u
from astropy import constants

from ..models import ModifiedNFW
from ..models import MB04
from ..models import YF17


def test_YF17():
    yf17 = YF17()
    ne = yf17.ne((0.,0.,20.))

def test_MB04():
    mb04 = MB04()
    ne = mb04.ne((0.,0.,20.))


def test_modified_NFW():
    # Init
    mNFW = ModifiedNFW()
    # xyz
    xyz = (1-2*rand(3, 100)) * 100
    # rho
    rho = mNFW.rho_b(xyz)
    assert rho.size == 100
    assert rho.unit == u.g/u.cm**3
    # nH
    nH = mNFW.nH(xyz)
    assert rho.size == 100
    xyz0 = [100., 0., 0.]
    nH0 = mNFW.nH(xyz0)
    assert np.isclose(nH0, 0.000279, rtol=1e-3)
    # ne
    ne = mNFW.ne(xyz)
    assert np.all(ne > nH)


