# Module to run tests on usage IGMGalaxyField

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from astropy import units as u
from astropy.table import Table, Column
from pyigm.field.igmfield import IgmGalaxyField

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


def create_fake_target_table():
    targets = Table()
    ra = np.arange(0, 10, 1)
    dec = ra
    # add duplicates
    for i in range(3):
        ra = np.append(ra, ra[:3])
        dec = np.append(dec, dec[:3])

    targets.add_column(Column(data=ra, name='TARG_RA'))
    targets.add_column(Column(data=dec, name='TARG_DEC'))
    return targets


def test_clean_duplicates():
    field = ('PG1407+265', 212.349634*u.deg, 26.3058650*u.deg)
    lfield = IgmGalaxyField((field[1],field[2]), name=field[0], verbose=False)

    targets = create_fake_target_table()
    new_target = lfield.clean_duplicates(targets)
    np.testing.assert_allclose(new_target['TARG_RA'], np.arange(0,10,1))


def test_calc_rho():
    field = ('PG1407+265', 212.349634*u.deg, 26.3058650*u.deg)
    lfield = IgmGalaxyField((field[1],field[2]), name=field[0], verbose=False)
    # Objects
    obj = dict(RA=np.array([212.349634]*3),
               DEC=np.array([26.3057]*3),
               Z=np.array([0.1, 0.2, 0.3]))
    rho = lfield.calc_rhoimpact(obj)
    np.isclose(rho[0].value, 1.244429196296971)
