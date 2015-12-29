# Module to run tests on generating IGMSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

from ..cgm import CGM, CGMAbsSys
from ..cgmsurvey import CGMAbsSurvey
from pyigm.field.galaxy import Galaxy
from pyigm.abssys.igmsys import IGMSystem

import pdb

def test_init_cgm():
    # Simple properties
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gal = Galaxy(radec,z=0.3)
    cgm = CGM(gal)
    # Test
    np.testing.assert_allclose(cgm.galaxy.z, 0.3)


def test_init_cgmabssys():
    radec = (125*u.deg, 45.2*u.deg)
    gal = Galaxy(radec,z=0.3)
    radec_qso = (125*u.deg, 45.203*u.deg)
    igmsys = IGMSystem('CGM', radec_qso, gal.z, [-500,500]*u.km/u.s)
    # Instantiate
    cgmabs = CGMAbsSys(gal, igmsys)
    # Test
    np.testing.assert_allclose(cgmabs.rho.value, 48.72077748027017)

def test_init_cgmabssurvey():
    cgmsurvey = CGMAbsSurvey(survey='cos-halos', ref='Tumlinson+13, Werk+13')
    # Test
    assert cgmsurvey.survey == 'cos-halos'