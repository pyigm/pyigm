import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from linetools.spectra.io import get_table_column


def check_dup_coord(coords, tol=1*u.arcsec):
    """ Checks for duplicates in a list or array
    of coordinates, within an angular distance
    tolerance `tol`.

    Parameters
    ----------
    coords : SkyCoord array or list
        Coordinates to check for duplicates
    tol : Angle, optional
        Angular tolerance

    Returns
    -------
    isdup : boolean array, shape(len(coords),)
        Whether a given coordinate in `coords` is a duplicate
    idx : array, shape(len(coords),)
        Indices of the first nearest neighbour as
        given by astropy.coord.match_coordinates_sky()

    """

    idx, d2d, d3d = match_coordinates_sky(coords, coords, nthneighbor=2)
    isdup = d2d < tol

    return isdup, idx


def check_dup_table(table, tol=1*u.arcsec, ra_tags=['TARG_RA','ra','RA'], dec_tags=['TARG_DEC','dec','DEC']):
    """
    Checks for duplicates in a Table based on their
    sky coordinates given by RA and DEC.

    Parameters
    ----------
    table : Table
        Table to check duplicates for
    tol : Angle, optional
        Tolerance for the duplicates search
    ra_tags : str, list; optional
        Name(s) of the column corresponding to RA.
        If not given, it guesses it.
    dec_tags : str, list; optional
        Name(s) of the column corresponding to DEC.
        If not given, it guesses it.

    Returns
    -------
    isdup : boolean array, shape(len(table),)
        Whether the row in `table` is duplicate
    idx : array, shape(len(table),)
        Indices of the first nearest neighbour as
        given by astropy.coord.match_coordinates_sky()
    coords : SkyCoord array for the table

    """

    ra, _ = get_table_column(ra_tags, [table], idx=0)
    dec, _ = get_table_column(dec_tags, [table], idx=0)

    coords = SkyCoord(ra, dec, unit='deg')

    isdup, idx = check_dup_coord(coords, tol=tol)

    return isdup, idx, coords
