# Module to run tests on generating IGMSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

from pyigm.abssys.igmsys import IGMSystem, HISystem

import pdb

def test_init():
    # Simple properties
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    igmsys = IGMSystem(radec, 1.244, [-500,500]*u.km/u.s, NHI=16., abs_type='MgII')
    # Test
    np.testing.assert_allclose(igmsys.zabs,1.244)
    #
    hisys = HISystem(radec, 1.244, [-500,500]*u.km/u.s, NHI=15.)
    np.testing.assert_allclose(hisys.NHI, 15.)
    np.testing.assert_allclose(hisys.vlim[0].value, -500.)


