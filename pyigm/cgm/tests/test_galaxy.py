# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from ..galaxy import GalaxyCGM


def test_init():
    mwcgm = GalaxyCGM()
    assert len(mwcgm.cgm_abs) > 0
    assert 'Fang+15' in mwcgm.refs

