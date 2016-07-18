# Module to run tests on Mocks of the Forest

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb

from astropy import units as u
from astropy.cosmology import Planck15

from pyigm.fN.mockforest import mk_mock
from pyigm.fN.fnmodel import FNModel



def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_mk_mock():
    # Quasar
    zem = 2.5
    # Spectral properties
    s2n = 10.
    sampling = 2.
    R = 2000.
    # Resultant wavelength array (using constant dwv instead of constant dv)
    disp = 4000/R/sampling # Angstrom
    wave = np.arange(3800., 1300*(1+zem), disp)*u.AA
    # f(N)
    fN_model = FNModel.default_model(cosmo=Planck15)

    #
    mock_spec, HI_comps, _ = mk_mock(wave, zem, fN_model, s2n=s2n,
                                     fwhm=sampling, seed=11223)
    np.testing.assert_allclose(mock_spec.flux[300].value, 21.77, rtol=0.05)  # scipy 0.17
