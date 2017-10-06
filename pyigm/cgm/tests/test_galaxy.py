# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from pyigm.cgm.galaxy import GalaxyCGM


def test_init():
    mwcgm = GalaxyCGM()
    assert len(mwcgm.abs.cgm_abs) > 0
    assert 'Fang+15' in mwcgm.refs
    # OVII table
    ovii_tbl = mwcgm.abs.ion_tbl((8,7))
    assert len(ovii_tbl['sig_logN'][0]) == 2
    # OVI table
    ovi_tbl = mwcgm.abs.ion_tbl((8,6))
    # Usage
    coords = mwcgm.abs.scoord
    #pytest.set_trace()

