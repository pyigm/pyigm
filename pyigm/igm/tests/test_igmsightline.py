# Module to run tests on IGMSightline

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import os
from pkg_resources import resource_filename

from linetools import utils as ltu

from pyigm.igm.igmsightline import IGMSightline


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_make_igmsystems():
    # Load a sightline
    sl_file = resource_filename('pyigm', '/data/sightlines/Blind_CIV/J161916.55+334238.41.json')
    sl_dict = ltu.loadjson(sl_file)
    igmsl = IGMSightline.from_dict(sl_dict)
    # Make them
    igmsl._abssystems = igmsl.make_igmsystems()
    assert len(igmsl._abssystems) == 2
    # Write
    igmsl.write_to_json(data_path('tst_sl.json'))
    # Read
    igmsl2 = IGMSightline.from_json(data_path('tst_sl.json'))


def test_from_igmguesses_and_write_igmguesses():
    igms = IGMSightline.from_igmguesses(data_path('J1410+2304_model.json'))
    # Test
    comps = igms._components
    assert comps[0].name == 'CIV_z-0.00024'
    assert comps[0].reliability == 'a'
    assert comps[8].zcomp == -0.0001
    assert len(comps) == 132
    assert len(comps[119]._abslines) == 29

    # write
    # will write a file in directory ./files/
    igms.write_to_igmguesses(outfile=data_path('IGM_model.json'), overwrite=True, date='2018-Feb-12')
    d1 = ltu.loadjson(data_path('IGM_model.json'))
    d2 = ltu.loadjson(data_path('J1410+2304_model.json'))
    added, removed, modified, same = ltu.compare_two_dict(d1,d2)
    pytest.set_trace()
    assert ltu.compare_two_json(data_path('IGM_model.json'), data_path('J1410+2304_model.json'))
