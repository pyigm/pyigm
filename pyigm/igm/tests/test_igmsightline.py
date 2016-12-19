# Module to run tests on IGMSightline

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from linetools import utils as ltu

from pyigm.igm.igmsightline import IGMSightline
import pyigm

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_make_igmsystems():
    # Load a sightline
    sl_file = pyigm.__path__[0]+'/data/sightlines/Blind_CIV/J161916.55+334238.41.json'
    sl_dict = ltu.loadjson(sl_file)
    igmsl = IGMSightline.from_dict(sl_dict)
    # Make them
    igm_sys = igmsl.make_igmsystems()
    assert len(igm_sys) == 2

