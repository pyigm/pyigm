# Module to run tests on initializing AbsSurvey

# TEST_UNICODE_LITERALS

import numpy as np
import glob, os, imp, pdb
import pytest

from pyigm.surveys.llssurvey import LLSSurvey
from ..igmsurvey import IGMSurvey

lt_path = imp.find_module('linetools')[1]

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


def test_init_igmsurvey():
    igmsurvey = IGMSurvey('DLA')
    assert igmsurvey.nsys == 0


def test_gz():
    # All
    z3mage = LLSSurvey.load_mage_z3()
    zeval, gz = z3mage.calculate_gz()
    assert gz[4000] == 67
    np.testing.assert_allclose(zeval[4000], 2.8705998897560931)


