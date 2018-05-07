# Module to run tests on initializing IGMGalaxyField

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import pdb

from pyigm.clustering.randist import RanDist

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''
seed = 101

def test_randist():
    # check a column density power law distribution
    beta = 1.5

    def ndist(n):
        return n ** -beta

    nvals = np.logspace(12.6, 16, 1000)
    ran = RanDist(nvals, ndist(nvals))
    #ran.self_test(log=1, seed=seed)
    y = ran.random(10000, seed=seed)
    assert np.isclose(y[0], 6.41122601e12)

    bsig = 24.0

    # b param distribution
    def bdist(b):
        b1 = bsig / b
        return b1 ** 5 * np.exp(-b1 ** 4)

    bvals = np.linspace(10, 150, 1000)
    ran2 = RanDist(bvals, bdist(bvals))
    #ran.self_test(seed=seed)
    y2 = ran2.random(10000, seed=seed)
    assert np.isclose(np.median(y2), 26.275860489355459)

    gamma = 2.04

    # z distribution
    def zdist(z):
        return (1 + z) ** gamma

    zp1vals = np.logspace(np.log10(1 + 2.5), np.log10(1 + 4), 1000)
    ran = RanDist(zp1vals, zdist(zp1vals))
    #ran.self_test(log=1, seed=seed, N=1e5, nbins=20)
    #pl.show()
    #'''

