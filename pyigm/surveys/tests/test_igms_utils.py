# Module to run tests on initializing AbsSurvey

# TEST_UNICODE_LITERALS

import numpy as np
import pdb
import pytest

from ..igmsurvey import IGMSurvey
from ..llssurvey import LLSSurvey
from ..dlasurvey import DLASurvey
from pyigm.surveys import utils as pyisu

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


def test_init_by_type():
    # LLSSurvey
    llssurvey = pyisu.class_by_type('LLS')()
    assert isinstance(llssurvey, LLSSurvey)
    # DLASurvey
    dlasurvey = pyisu.class_by_type('DLA')()
    assert isinstance(dlasurvey, DLASurvey)

