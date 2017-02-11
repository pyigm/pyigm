# Module to run tests on generating IGMSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import io
import json

from astropy import units as u
from astropy.coordinates import SkyCoord

from ..cgm import CGM, CGMAbsSys
from ..cgmsurvey import CGMAbsSurvey
from ..utils import cgm_from_galaxy_igmsystems
from pyigm.field.galaxy import Galaxy
from pyigm.abssys.igmsys import IGMSystem
from pyigm.igm.igmsightline import IGMSightline
import pyigm

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
    igmsys = IGMSystem(radec_qso, gal.z, [-500,500]*u.km/u.s, abs_type='CGM')
    # Instantiate
    cgmabs = CGMAbsSys(gal, igmsys)
    # Test
    np.testing.assert_allclose(cgmabs.rho.value, 48.72077748027017)

def test_init_cgmabssurvey():
    cgmsurvey = CGMAbsSurvey(survey='cos-halos', ref='Tumlinson+13, Werk+13')
    # Test
    assert cgmsurvey.survey == 'cos-halos'

def test_to_dict():
    radec = (125*u.deg, 45.2*u.deg)
    gal = Galaxy(radec,z=0.3)
    radec_qso = (125*u.deg, 45.203*u.deg)
    igmsys = IGMSystem(radec_qso, gal.z, [-500,500]*u.km/u.s, abs_type='CGM')
    # Instantiate
    cgmabs = CGMAbsSys(gal, igmsys)
    # Test
    cdict = cgmabs.to_dict()
    with io.open('tmp.json', 'w', encoding='utf-8') as f:
        f.write(unicode(json.dumps(cdict, sort_keys=True, indent=4,
                                   separators=(',', ': '))))

def test_cgm_from_igmsystems():
    # Load sightlines
    sl_file = pyigm.__path__[0]+'/data/sightlines/Blind_CIV/J115120.46+543733.08.json'
    igmsl = IGMSightline.from_json(sl_file)
    igmsys = igmsl.make_igmsystems(vsys=400*u.km/u.s)
    # Galaxy
    galaxy = Galaxy((178.84787, 54.65734), z=0.00283)
    # Go
    cgm_list = cgm_from_galaxy_igmsystems(galaxy, igmsys, correct_lowz=False)
    assert len(cgm_list) == 1
    np.testing.assert_allclose(cgm_list[0].rho.value, 127.8324005876)

def test_cgm_from_igmsystems_lowz():
    # Load sightlines
    sl_file = pyigm.__path__[0]+'/data/sightlines/Blind_CIV/J115120.46+543733.08.json'
    igmsl = IGMSightline.from_json(sl_file)
    igmsys = igmsl.make_igmsystems(vsys=400*u.km/u.s)
    # Galaxy
    galaxy = Galaxy((178.84787, 54.65734), z=0.00283)
    # Go
    cgm_list = cgm_from_galaxy_igmsystems(galaxy, igmsys, correct_lowz=True)
    assert len(cgm_list) == 1
    np.testing.assert_allclose(cgm_list[0].rho.value, 189.01377)
