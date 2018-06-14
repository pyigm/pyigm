from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from astropy.table import Table


def make_rand_gal(nrand=1000):
    # Make a set of random, uniformly distributed galaxies
    RA = 200. + 1.0 * np.random.uniform(size=nrand)  # About 20Mpc at z~0.3
    DEC = 25. + 1.0 * np.random.uniform(size=nrand)
    Z = 0.3 + 0.2*np.random.uniform(size=nrand)
    MAG = 19.0 + 2.*np.random.uniform(size=nrand)

    tbl = Table()
    tbl['RA'] = RA
    tbl['DEC'] = DEC
    tbl['ZGAL'] = Z
    tbl['MAG'] = MAG

    return tbl

def make_edges():
    # lower edge of smallest bin, bin width, and number of bins
    # in comoving Mpc. For radial bins
    rmin = 0.
    rmax = 20
    rwidth = 0.5
    rbinedges = np.arange(rmin, rmax + 0.5 * rwidth, rwidth)
    rcbins = 0.5 * (rbinedges[:-1] + rbinedges[1:])

    # for transverse bins
    tmin = 0.
    tmax = 10
    twidth = 0.5
    tbinedges = np.arange(tmin, tmax + 0.5 * twidth, twidth)
    # tbinedges = np.append(tbinedges,[6**i for i in np.linspace(1.1,2.2,20)])
    tcbins = 0.5 * (tbinedges[:-1] + tbinedges[1:])

    # Bins for integration limits of projected 2D2PCFs
    nbins_T = int(1. / twidth)
    #nbins_LOS = None

    # Return
    return rbinedges, tbinedges, nbins_T
