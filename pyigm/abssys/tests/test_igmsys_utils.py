# Module to run tests on generating IGMSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from astropy import units as u

from linetools.spectra import io as lsio
import linetools

from ..lls import LLSSystem
from ..dla import DLASystem

from pyigm.abssys import utils as pyasu

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_init_by_type():
    # LLSSurvey
    system = pyasu.class_by_type('LLS')((0.*u.deg, 0.*u.deg),
                                        2.0, None, NHI=17.9)
    assert isinstance(system, LLSSystem)
    # DLASurvey
    system = pyasu.class_by_type('DLA')((0.*u.deg, 0.*u.deg),
                                        2.5, None, NHI=20.55)
    assert isinstance(system, DLASystem)


def test_hi_lya():
    # Simple system (without an absline)
    dla1 = DLASystem.from_json(data_path('J010311.38+131616.7_z2.309_ESI.json'))
    dla1.zabs = 4.
    dla1.NHI = 21.
    dla2 = DLASystem.from_json(data_path('J010311.38+131616.7_z2.309_ESI.json'))
    dla2.zabs = 3.5
    dla2.NHI = 20.5
    spec_fil = linetools.__path__[0]+'/spectra/tests/files/PH957_f.fits'
    spec = lsio.readspec(spec_fil)
    # Just lya
    model, lya_lines = pyasu.hi_model([dla1, dla2], spec, lya_only=True)
    ipx = np.argmin(np.abs(spec.wavelength.value-(1+dla1.zabs)*1215.67))
    assert model.flux[ipx].value < 1e-4
    ipx2 = np.argmin(np.abs(spec.wavelength.value-(1+dla1.zabs)*1025.7222))
    assert model.flux[ipx2].value > 0.99
    # Lyman series
    model2, lyman_lines = pyasu.hi_model([dla1, dla2], spec)
    assert len(lyman_lines) > 10
    assert model2.flux[ipx2].value < 1e-4
    # LLS
    model3, lyman_lines = pyasu.hi_model([dla1, dla2], spec, add_lls=True)
    assert model3.flux[0].value < 1e-4
