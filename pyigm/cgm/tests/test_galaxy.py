# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from pyigm.cgm.galaxy import GalaxyCGM


def test_init_light():
    mwcgm = GalaxyCGM(load=False)

def test_init_full():
    mwcgm = GalaxyCGM()
    assert len(mwcgm.abs.cgm_abs) > 0

    # Cool
    assert 'Richter+17' in mwcgm.refs
    SiII_tbl = mwcgm.abs.ion_tbl((14,2))

    assert (not np.any(np.isnan(SiII_tbl['logN'])))
    assert np.sum(SiII_tbl['flag_N'] > 0) == 188

    # Hot
    assert 'Fang+15' in mwcgm.refs
    # OVII table
    ovii_tbl = mwcgm.abs.ion_tbl((8,7))
    assert len(ovii_tbl['sig_logN'][0]) == 2
    # OVI table
    ovi_tbl = mwcgm.abs.ion_tbl((8,6))
    # Usage
    coords = mwcgm.abs.scoord

