# Module to run tests on initializing AbslineSystem

# TEST_UNICODE_LITERALS

import numpy as np
import copy

import astropy.units as u

from pyigm.fN.fnmodel import FNModel
from pyigm.fN import tau_eff as pyteff

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_teff():
    fN_default = FNModel.default_model()
    zval,teff_LL = pyteff.lyman_limit(fN_default, 0.5, 2.45)
    #
    np.testing.assert_allclose(zval[0], 0.5)
    np.testing.assert_allclose(teff_LL[0], 1.8176161746504436)

def test_lya():
    # f(N)
    fN_model = FNModel.default_model()

    # tau_eff
    lamb = 1215.6701*(1+2.4)
    teff = pyteff.lyman_ew(lamb, 2.5, fN_model, NHI_MIN=12., NHI_MAX=17.)
    # Test
    np.testing.assert_allclose(teff, 0.19821452949713764)

def test_lyx():
    # f(N)
    fN_model = FNModel.default_model()

    # tau_eff
    lamb = 917.*(1+2.4)
    teff = pyteff.lyman_ew(lamb, 2.5, fN_model, NHI_MIN=12., NHI_MAX=17.)
    # Test
    np.testing.assert_allclose(teff, 0.234694549996041)

def test_parallel():
    import multiprocessing
    from linetools.lists.linelist import LineList
    # f(N)
    fN_model = FNModel.default_model()
    # Lines
    HI = LineList('HI')
    tst_wv = HI._data['wrest']
    #
    adict = []
    for wrest in tst_wv:
        tdict = dict(ilambda=wrest.value*(1+2.4), zem=2.5, fN_model=fN_model, wrest=copy.deepcopy(tst_wv))
        adict.append(tdict)

    pool = multiprocessing.Pool(2) # initialize thread pool N threads
    ateff = pool.map(pyteff.map_lymanew, adict)
    # Test
    np.testing.assert_allclose(ateff[-3:],
                               np.array([0.23440858789182742,
                                0.20263221240650739, 0.21927057420866358]))

def test_DM():
    DM = pyteff.DM(1.)
    assert DM.unit == u.pc/u.cm**3
    np.testing.assert_allclose(DM.value, 1235.8727960316237)
