# Module to run tests on IGMSightline


import pytest
import numpy as np
import os
from pyigm.resources import resource_filename

from linetools import utils as ltu

from pyigm.igm.igmsightline import IGMSightline


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_make_igmsystems(tmp_path):
    # Load a sightline
    sl_file = resource_filename('pyigm', '/data/sightlines/Blind_CIV/J161916.55+334238.41.json')
    sl_dict = ltu.loadjson(sl_file)
    igmsl = IGMSightline.from_dict(sl_dict)
    # Make them
    igmsl._abssystems = igmsl.make_igmsystems()
    assert len(igmsl._abssystems) == 2
    # Write
    outfil = str(tmp_path / 'tst_sl.json')
    igmsl.write_to_json(outfil)
    # Read
    igmsl2 = IGMSightline.from_json(outfil)


def test_from_igmguesses_and_write_igmguesses(tmp_path):
    igms = IGMSightline.from_igmguesses(data_path('J1410+2304_model.json'))

    # Test
    comps = igms._components
    assert comps[0].name == 'CIV_z-0.00024'
    assert comps[0].reliability == 'a'
    assert np.isclose(comps[1].zcomp,  -0.00024)
    assert len(comps) == 2
    assert len(comps[1]._abslines) == 2

    # write (to a temporary directory, not into the installed package)
    out_model = str(tmp_path / 'IGM_model.json')
    igms.write_to_igmguesses(outfile=out_model, overwrite=True, date='2018-Feb-12')
    #d1 = ltu.loadjson(data_path('IGM_model.json'))
    #d2 = ltu.loadjson(data_path('J1410+2304_model.json'))

    # assert ltu.compare_two_files(data_path('IGM_model.json'), data_path('J1410+2304_model.json'), verbose=True)
    assert ltu.compare_two_json(out_model, data_path('J1410+2304_model.json'))
