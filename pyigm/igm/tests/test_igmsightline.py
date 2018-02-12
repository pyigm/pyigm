# Module to run tests on IGMSightline

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import os

from linetools import utils as ltu

from pyigm.igm.igmsightline import IGMSightline
import pyigm


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_make_igmsystems():
    # Load a sightline
    sl_file = pyigm.__path__[0]+'/data/sightlines/Blind_CIV/J161916.55+334238.41.json'
    sl_dict = ltu.loadjson(sl_file)
    igmsl = IGMSightline.from_dict(sl_dict)
    # Make them
    igm_sys = igmsl.make_igmsystems()
    assert len(igm_sys) == 2


# def test_from_igmguesses():
if 1:
    igms = IGMSightline.from_igmguesses(data_path('IGM_model_reference.json'))
    # Test
    assert len(igms._components) == 2



def test_read_write_igmg():
    # read
    igmg_file = data_path('J1410+2304_model.json')
    comps = ltiio.read_igmg_to_components(igmg_file)
    assert comps[0].name == 'CIV_z-0.00024'
    assert comps[0].reliability == 'a'
    assert comps[8].zcomp == -0.0001
    assert len(comps) == 132
    # write
    # will write a file in directory ./files/
    ltiio.write_igmg_from_components(comps[:2], specfile='test.fits', fwhm=3, outfile=data_path('IGM_model.json'))
    compare_two_files(data_path('IGM_model.json'),
                      resource_filename('linetools', '/data/tests/IGM_model_reference.json'), except_l2_has='2018-Feb-09')

