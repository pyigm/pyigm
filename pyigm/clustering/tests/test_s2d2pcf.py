# Module to run tests on initializing Survey2D2PCF

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import pdb

from astropy.coordinates import SkyCoord

from pyigm.clustering.survey_2d2pcf import Survey2D2PCF
from pyigm.clustering.clustering_field import ClusteringField
from pyigm.clustering.tests import utils

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''
seed = 101

def test_init_survey():


    # Edges (~20 Mpc)
    rbinedges, tbinedges, nbins_T = utils.make_edges()

    # Random galaxies
    rand_gal = utils.make_rand_gal(nrand=1000)

    # Now Clustering
    cf = ClusteringField(SkyCoord(ra=200.5, dec=25.5, unit='deg'))
    cf.galaxies = rand_gal

    # Add em
    gal_idx = rand_gal['ZGAL'] > 0.31
    cf.addGal(gal_idx)

    # Count pairs
    cf.compute_pairs(tbinedges, rbinedges)

    # Now the Survey
    csurvey = Survey2D2PCF(cf)

    # Test
    assert csurvey.DgDg.shape == (40,20)

    # Normalize
    #csurvey.set_normalization(True, Ngal_rand=csurvey.RgRg.size)
    csurvey.set_normalization()

    # Wgg
    Wgg, err_Wgg = csurvey.calc_xi_gg()

    # Transverse
    s=1
    csurvey.calc_xi_gg_transverse(s)

    # Test?
    assert np.all(csurvey.xi_gg_T < 1.)  # A weird random event could lead to a high value and a failed test..
