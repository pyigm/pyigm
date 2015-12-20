# Module to run tests on FNConstaint

# TEST_UNICODE_LITERALS

import numpy as np
import imp, pdb

from pyigm.fN.constraints import FNConstraint

pyigm_path = imp.find_module('pyigm')[1]


def test_init():
    tst = FNConstraint('fN')
    assert tst.fN_dtype == 'fN'

def test_init_from_fits():
    """ Load FNConstraint from a set of FITS files
    """
    # Load
    fn_file = pyigm_path+'/data/fN/fn_constraints_z2.5_vanilla.fits'
    k13r13_file = pyigm_path+'/data/fN/fn_constraints_K13R13_vanilla.fits'
    n12_file = pyigm_path+'/data/fN/fn_constraints_N12_vanilla.fits'
    all_fN_cs = FNConstraint.from_fitsfile([fn_file,k13r13_file, n12_file])
    # Test
    assert len(all_fN_cs) == 8
    fN_dtype = [fc.fN_dtype for fc in all_fN_cs]
    itau = fN_dtype.index('teff')
    np.testing.assert_allclose(all_fN_cs[itau].data['TEFF'], 0.19778812)


