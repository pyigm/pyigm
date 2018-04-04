# Module to run tests on CGMAbsSurvey
#  Uses COSHalos for testing

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb

import pytest
remote_data = pytest.mark.remote_data

from astropy.table import Table

from pyigm.cgm.cos_halos import COSHalos#, COSDwarfs

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_data_for_survey():
    cos_halos = COSHalos()
    assert isinstance(cos_halos._data, Table)
    # Add an ion
    cos_halos.add_ion_to_data('OVI')
    assert 'flag_N_OVI' in cos_halos._data.keys()



