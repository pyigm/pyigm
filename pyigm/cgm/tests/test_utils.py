# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from astropy.cosmology import Planck15 as cosmo
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy

from pyigm.cgm.utils import calc_rho
from pyigm.field.galaxy import Galaxy
from pyigm.abssys.igmsys import IGMSystem


def test_calcrho():
    # Dummy
    galaxy = Galaxy((100.,50.), 0.2)
    igmsys = IGMSystem((100.,50.001), 1.2, None)
    # Calc
    rho, angle = calc_rho(galaxy, igmsys, cosmo)
    # Test
    assert np.isclose(rho.value, 12.2587523534)
    assert rho.unit == astropy.units.kpc
    assert isinstance(angle, astropy.coordinates.Angle)

