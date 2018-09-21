# Module to run tests on initializing IGMGalaxyField

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import pdb

from astropy.coordinates import SkyCoord

from pyigm.clustering.clustering_field import ClusteringField
from pyigm.clustering.tests import utils

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''
seed = 101

def test_init_cfield():

    rand_gal = utils.make_rand_gal(nrand=1000)

    # Now Clustering
    cf = ClusteringField(SkyCoord(ra=200.5, dec=25.5, unit='deg'))
    cf.galaxies = rand_gal

    # Add em
    gal_idx = rand_gal['ZGAL'] > 0.31
    cf.addGal(gal_idx)

    # Test
    assert len(cf.galrand) == 10*np.sum(gal_idx)

