# Module to run tests on initializing IGMGalaxyField

# TEST_UNICODE_LITERALS

import numpy as np
import pdb

from astropy import units as u

from pyigm.field.igmfield import IgmGalaxyField

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_init():
    field = ('PG1407+265',212.349634*u.deg,26.3058650*u.deg)
    lfield = IgmGalaxyField((field[1],field[2]), name=field[0], verbose=True)
    print(lfield)

