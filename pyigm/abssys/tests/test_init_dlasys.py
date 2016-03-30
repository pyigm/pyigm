# Module to run tests on initializing DLA System

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent

from pyigm.abssys.dla import DLASystem

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


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


def test_DLA_from_components():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lya.attrib['N'] = 3e20 / u.cm**2
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    lyb.attrib['N'] = 3e20 / u.cm**2
    abscomp = AbsComponent.from_abslines([lya,lyb])
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

