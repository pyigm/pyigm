# Module to run tests on mcmc

# TEST_UNICODE_LITERALS

import numpy as np
import pytest
import pdb

from pyigm.fN import mcmc


def test_z3():
    """ Load FNConstraint from a set of FITS files
    """
    # Load
    all_fN_cs = mcmc.set_fn_data()
    assert len(all_fN_cs) == 17

def test_z4():
    """ Load FNConstraint from a set of FITS files
    """
    # Load
    all_fN_cs = mcmc.set_fn_data(flg=4)
    assert len(all_fN_cs) == 12
