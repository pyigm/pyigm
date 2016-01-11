# Module to run tests on FNConstaint

# TEST_UNICODE_LITERALS

import numpy as np
import imp, pdb

from pyigm.fN.constraints import FNConstraint

def test_init():
    tst = FNConstraint('fN')
    assert tst.fN_dtype == 'fN'

def test_init_from_fits():
    """ Load FNConstraint from a set of FITS files
    """
    # Load
    all_fN_cs = FNConstraint.load_defaults()
    # Test
    assert len(all_fN_cs) == 8
    fN_dtype = [fc.fN_dtype for fc in all_fN_cs]
    itau = fN_dtype.index('teff')
    np.testing.assert_allclose(all_fN_cs[itau].data['TEFF'], 0.19778812)


