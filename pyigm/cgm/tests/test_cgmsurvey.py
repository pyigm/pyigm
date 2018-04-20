# Module to run tests on CGMAbsSurvey
#  Uses COSHalos for testing

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb

import pytest
remote_data = pytest.mark.remote_data

from astropy.table import Table

from pyigm.cgm.cos_halos import COSHalos#, COSDwarfs
from pyigm.cgm.cgmsurvey import CGMAbsSurvey

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_write_to_json():
    cos_halos = COSHalos()
    cos_halos.to_json(data_path('cos_halos.json'))

def test_from_json():
    # Without building systems
    cos_halos = CGMAbsSurvey.from_json(data_path('cos_halos.json'))
    assert len(cos_halos.cgm_abs) == 0
    # Build the systems too
    cos_halos2 = CGMAbsSurvey.from_json(data_path('cos_halos.json'), build_sys=True)
    assert len(cos_halos2.cgm_abs) == 44

def test_data_for_survey():
    cos_halos = COSHalos()
    assert isinstance(cos_halos._data, Table)
    # Add an ion
    cos_halos.add_ion_to_data('OVI')
    assert 'flag_N_OVI' in cos_halos._data.keys()



