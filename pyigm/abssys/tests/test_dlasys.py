# Module to run tests on initializing DLA System

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent
from linetools.spectra import io as lsio
from linetools.abund.relabund import RelAbund
import linetools

from pyigm.abssys.dla import DLASystem


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_simple_dla_init():
	# Init 
    dla = DLASystem((0.*u.deg, 0.*u.deg), 2.5, None, NHI=20.55)
    #
    np.testing.assert_allclose(dla.vlim[0].value,-500.)
    np.testing.assert_allclose(dla.NHI, 20.55)


def test_dat_init():
    # JXP .dat files
    if os.getenv('DLA') is None:
        assert True
        return
    # Read
    datfil = 'Data/PH957.z2309.dat'
    dla = DLASystem.from_datfile(datfil, tree=os.environ.get('DLA'))
    #
    np.testing.assert_allclose(dla.NHI, 21.37)
    np.testing.assert_allclose(dla.zabs, 2.309)


def test_model_abs():
    # Simple system (without an absline)
    dla = DLASystem.from_json(data_path('J010311.38+131616.7_z2.309_ESI.json'))
    spec_fil = linetools.__path__[0]+'/spectra/tests/files/PH957_f.fits'
    spec = lsio.readspec(spec_fil)
    model, lya_lines = dla.model_abs(spec)
    # import pdb; pdb.set_trace()
    # Check core
    ipx = np.argmin(np.abs(spec.wavelength.value-(1+dla.zabs)*1215.67))
    assert model.flux[ipx].value < 1e-4


def test_parse_ion():
    # JXP .ion file
    if os.getenv('DLA') is None:
        assert True
        return
    # Read
    datfil = 'Data/PH957.z2309.dat'
    dla = DLASystem.from_datfile(datfil, tree=os.environ.get('DLA'))
    #
    dla.get_ions(use_Nfile=True)
    assert len(dla._ionN) == 14


def test_dla_from_dict():
    dla = DLASystem.from_json(data_path('J010311.38+131616.7_z2.309_ESI.json'))
    assert len(dla._components) == 16


def test_DLA_from_components():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA, z=2.92939)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['flag_N'] = 1
    lya.attrib['N'] = 3e20 / u.cm**2
    lya.attrib['sig_N'] = 1 / u.cm**2
    lyb = AbsLine(1025.7222*u.AA, z=lya.z)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['N'] = 3e20 / u.cm**2
    lyb.attrib['flag_N'] = 1
    lyb.attrib['sig_N'] = 1 / u.cm**2
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.synthesize_colm()
    abscomp.coord = radec
    # Instantiate
    HIsys = DLASystem.from_components([abscomp])
    # Test
    np.testing.assert_allclose(HIsys.NHI, 20.477121254719663)
    assert len(HIsys._components) == 1
    assert HIsys._components[0].Zion[0] == 1
    assert HIsys._components[0].Zion[1] == 1

"""
def test_default_dla_sample():
    if os.getenv('DLA') is None:
        assert True
        return
    # Load
    dlas = DLASurvey.default_sample()
    assert len(dlas._abs_sys) == 100

def test_default_dla_sample_with_ions():
    if os.getenv('DLA') is None:
        assert True
        return
    # Load
    dlas = DLASurvey.default_sample()
    dlas.fill_ions(use_Nfile=True)
    CIV_clms = dlas.ions((6,4))
    gdCIV = np.where(CIV_clms['flag_N']>0)[0]
    assert len(gdCIV) == 74
"""

def test_dla_XY():
    spec_fil = linetools.__path__[0]+'/spectra/tests/files/PH957_f.fits'
    spec = lsio.readspec(spec_fil)
    dla = DLASystem.from_json(data_path('J010311.38+131616.7_z2.309_ESI.json'))
    #
    dla.measure_aodm(spec=spec)
    dla.update_component_colm()
    dla.fill_ionN()
    dla.XY = RelAbund.from_ionclm_table((1,21.37,0.08), dla._ionN)
    tbl = dla.XY.table()
    assert len(tbl) == 8



