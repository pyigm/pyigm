# Module to run tests on generating IGMSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

from ..lls import LLSSystem
from ..dla import DLASystem
from pyigm.abssys.igmsys import IGMSystem
from pyigm.abssys import utils as pyasu

import pdb


def test_init_by_type():
    # LLSSurvey
    system = pyasu.class_by_type('LLS')((0.*u.deg, 0.*u.deg),
                                        2.0, None, NHI=17.9)
    assert isinstance(system, LLSSystem)
    # DLASurvey
    system = pyasu.class_by_type('DLA')((0.*u.deg, 0.*u.deg),
                                        2.5, None, NHI=20.55)
    assert isinstance(system, DLASystem)

