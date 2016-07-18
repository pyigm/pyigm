# Module to run tests on initializing AbslineSystem

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb

import pytest
remote_data = pytest.mark.remote_data

from ..cos_halos import COSHalos#, COSDwarfs

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

'''
def test_load_kin():
    # Class
    cos_halos = COSHalos()
    cos_halos.load_mega(skip_ions=True)
    # Load kin
    cos_halos.load_abskin()
'''

def test_load_sngl():
    # Class
    cos_halos = COSHalos(fits_path=data_path(''), cdir=data_path(''), load=False)
    # Load
    cos_halos.load_single_fits(('J0950+4831', '177_27'))

def test_write_sngl():
    import io, json
    # Class
    cos_halos = COSHalos(fits_path=data_path(''), cdir=data_path(''), load=False)
    # Load
    cos_halos.load_single_fits(('J0950+4831', '177_27'))
    # Write to JSON
    cdict = cos_halos.cgm_abs[0].to_dict()
    with io.open('tmp.json', 'w', encoding='utf-8') as f:
        f.write(unicode(json.dumps(cdict, sort_keys=True, indent=4,
                                   separators=(',', ': '))))

"""
def test_load_sngl_dwarf():
    # Class
    cos_dwarfs = COSDwarfs(fits_path=data_path(''), kin_init_file='dum', cdir='dum')
    # Load
    cos_dwarfs.load_single( ('J0042-1037', '358_9'))
"""

@remote_data
def test_load_survey():
    # Class
    cos_halos = COSHalos(debug=True)
    assert len(cos_halos.cgm_abs) == 4
    # Load
    #cos_halos.load_mega()  # Only reads one file for the test, actually
    #cos_halos.load_mega(skip_ions=True)

def test_getitem():
    # Class
    cos_halos = COSHalos(fits_path=data_path(''), cdir='dum', load=False)
    # Load
    cos_halos.load_single_fits(('J0950+4831', '177_27'))
    # Get item
    cgm_abs = cos_halos[('J0950+4831','177_27')]
    assert cgm_abs.galaxy.field == 'J0950+4831'

