# Module to run tests on continuum codes

## # TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest

from pyigm.continuum import core as pyicc

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


def test_dict():
    # Init 
    cdict = pyicc.init_conti_dict(Norm=1.)

    assert isinstance(cdict,dict)
    np.testing.assert_allclose(cdict['Norm'], 1.)

