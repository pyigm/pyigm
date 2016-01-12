"""Tests for field.utils"""

from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
import numpy as np


def create_fake_target_table():
    targets = Table()
    ra = np.arange(0,10,1)
    dec = ra
    # add duplicates
    for i in range(3):
        ra = np.append(ra, ra[:3])
        dec = np.append(dec, dec[:3])

    targets.add_column(Column(data=ra, name='TARG_RA'))
    targets.add_column(Column(data=dec, name='TARG_DEC'))
    return targets


def test_check_dup_coord():
    from pyigm.field.utils import check_dup_coord
    table = create_fake_target_table()

    coords = SkyCoord(table['TARG_RA'], table['TARG_DEC'], unit='deg')
    isdup, idx = check_dup_coord(coords)

    assert np.sum(isdup) == 12
    assert all(isdup[:3])


def test_check_dup_table():
    from pyigm.field.utils import check_dup_table
    table = create_fake_target_table()

    isdup, idx, coords = check_dup_table(table)
    assert np.sum(isdup) == 12
    assert all(isdup[:3])
