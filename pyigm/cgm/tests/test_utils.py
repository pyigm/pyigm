# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from astropy.cosmology import Planck15 as cosmo
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import astropy

from pyigm.cgm.utils import calc_cgm_rho, cgmsurvey_from_sightlines_fields
from pyigm.cgm.cgmsurvey import CGMAbsSurvey
from pyigm.field.galaxy import Galaxy
from pyigm.field.igmfield import IgmGalaxyField
from pyigm.abssys.igmsys import IGMSystem
from pyigm.igm.igmsightline import IGMSightline

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.spectra.io import readspec


def test_calcrho():
    # Dummy
    galaxy = Galaxy((100.,50.), 0.2)
    igmsys = IGMSystem((100.,50.001), 1.2, None)
    # Calc
    rho, angle = calc_cgm_rho(galaxy, igmsys, cosmo)
    # Test
    assert np.isclose(rho.value, 12.2587523534)
    assert rho.unit == astropy.units.kpc
    assert isinstance(angle, astropy.coordinates.Angle)


def test_cgmsurvey_from_fields_sightlines():
    # Instantiate fields and add galaxies
    field1 = ('PG1407+265', 212.349634 * u.deg, 26.3058650 * u.deg)
    fieldobj1 = IgmGalaxyField((field1[1], field1[2]), name=field1[0], verbose=False)
    field2 = ('PG1148+549', 177.83526042 * u.deg, 54.625855 * u.deg)
    fieldobj2 = IgmGalaxyField((field2[1], field2[2]), name=field2[0], verbose=False)
    f1ras = [212.33310891,212.329875]*u.deg
    f1decs = [26.3722716934,26.3084391667]*u.deg
    f1zs = [0.974413514137,0.22016787529]
    f1gals = Table([f1ras,f1decs,f1zs],names=['RA','DEC','Z'])
    f2ras = [177.835229336, 178.536958333]*u.deg
    f2decs = [54.625851084, 54.6176108333]*u.deg
    f2zs = [0.0595485754311, 0.00385531364009]
    f2gals = Table([f2ras, f2decs, f2zs], names=['RA', 'DEC', 'Z'])
    fieldobj1.galaxies = f1gals
    fieldobj2.galaxies = f2gals
    fields = [fieldobj1,fieldobj2]

    # Instantiate sightlines
    sl1coord = SkyCoord(field1[1],field1[2],unit = 'deg')
    sl2coord = SkyCoord(field2[1], field2[2],unit = 'deg')
    comp11 = AbsComponent(sl1coord,(8,6),0.9743,[-200,200]*u.km/u.s)
    comp12 = AbsComponent(sl1coord, (1, 1), 0.9743, [-200, 200] * u.km / u.s)
    al11 = AbsLine('OVI 1031')
    al11.attrib['coord'] = sl1coord
    al11.analy['spec'] = 'files/J0042-1037.358_9.fits.gz'
    al12 = AbsLine('HI 1215')
    al12.attrib['coord'] = sl1coord
    al12.analy['spec'] = 'files/J0042-1037.358_9.fits.gz'
    comp11.add_absline(al11,chk_sep=False,chk_vel=False)
    comp12.add_absline(al12, chk_sep=False, chk_vel=False)
    sys1 = IGMSystem.from_components([comp11,comp12])
    sl1 = IGMSightline.from_systems([sys1])
    comp21 = AbsComponent(sl2coord, (6, 4), 0.0037, [-200, 200] * u.km / u.s)
    comp22 = AbsComponent(sl2coord, (1,1), 0.0037, [-200, 200] * u.km / u.s)
    al21 = AbsLine('CIV 1548')
    al21.attrib['coord'] = sl1coord
    al21.analy['spec'] = 'files/J0042-1037.358_9.fits.gz'
    al22 = AbsLine('HI 1215')
    al22.attrib['coord'] = sl1coord
    al22.analy['spec'] = 'files/J0042-1037.358_9.fits.gz'
    comp21.add_absline(al21, chk_sep=False, chk_vel=False)
    comp22.add_absline(al22, chk_sep=False, chk_vel=False)
    sys2 = IGMSystem.from_components([comp21,comp22])
    sl2 = IGMSightline.from_systems([sys2])
    sightlines = [sl1,sl2]

    # Run function
    csurvey = cgmsurvey_from_sightlines_fields(fields,sightlines)
    assert isinstance(csurvey,CGMAbsSurvey)


def test_cgmsurvey_from_fields_sightlines():
    # Instantiate fields and add galaxies
    field1 = ('PG1407+265', 212.349634 * u.deg, 26.3058650 * u.deg)
    fieldobj1 = IgmGalaxyField((field1[1], field1[2]), name=field1[0], verbose=False)
    field2 = ('PG1148+549', 177.83526042 * u.deg, 54.625855 * u.deg)
    fieldobj2 = IgmGalaxyField((field2[1], field2[2]), name=field2[0], verbose=False)
    f1ras = [212.33310891,212.329875]*u.deg
    f1decs = [26.3722716934,26.3084391667]*u.deg
    f1zs = [0.974413514137,0.22016787529]
    f1gals = Table([f1ras,f1decs,f1zs],names=['RA','DEC','Z'])
    f2ras = [177.835229336, 178.536958333]*u.deg
    f2decs = [54.625851084, 54.6176108333]*u.deg
    f2zs = [0.0595485754311, 0.00385531364009]
    f2gals = Table([f2ras, f2decs, f2zs], names=['RA', 'DEC', 'Z'])
    fieldobj1.galaxies = f1gals
    fieldobj2.galaxies = f2gals
    fields = [fieldobj1,fieldobj2]

    # Instantiate sightlines
    sl1coord = SkyCoord(field1[1],field1[2],unit = 'deg')
    sl2coord = SkyCoord(field2[1], field2[2],unit = 'deg')
    comp11 = AbsComponent(sl1coord,(8,6),0.9743,[-200,200]*u.km/u.s)
    comp12 = AbsComponent(sl1coord, (1, 1), 0.9743, [-200, 200] * u.km / u.s)
    sys1 = IGMSystem.from_components([comp11,comp12])
    sl1 = IGMSightline.from_systems([sys1])
    comp21 = AbsComponent(sl2coord, (6, 4), 0.0037, [-200, 200] * u.km / u.s)
    comp22 = AbsComponent(sl2coord, (1,1), 0.0037, [-200, 200] * u.km / u.s)
    sys1 = IGMSystem.from_components([comp21,comp22])
    sl2 = IGMSightline.from_systems([sys1])
    sightlines = [sl1,sl2]

    # Run function
    csurvey = cgmsurvey_from_sightlines_fields(fields,sightlines)
    assert isinstance(csurvey,CGMAbsSurvey)

