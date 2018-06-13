""" methods for Maximum Likelihood analyses on IGM surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb

import numpy as np


def powerlaw_loxz(lls, z_gz, gz, guess, zpivot, rnga=(-0.5,0.5), rngl=(-0.2,0.1), ngrid=100):
    """

    Ported from SDSS LLS paper
    Functional form  l = k (1+z/1+zpivot)^a

    Parameters
    ----------
    lls
    z_gz
    gz
    guess
    zpivot
    rnga
    rngl
    ngrid

    Returns
    -------

    """

    # Create vectors of alpha and l_0 values
    lvec = guess[0] * 10.**(rngl[0] + (rngl[1]-rngl[0])*np.arange(ngrid)/(ngrid-1))
    avec = guess[1] * 10.**(rnga[0] + (rnga[1]-rnga[0])*np.arange(ngrid)/(ngrid-1))

    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    #;; Positive term, i.e. h(z) evaluated at z_i
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    #;; Evaluate z_i and sum
    zi_nrm = (1+lls['zabs'])/(1+zpivot)
    sum_lnzi = np.sum(np.log(zi_nrm))

    #;; Weight by alpha
    asum_zi = avec * sum_lnzi

    #;; lvec sum is trivial
    nLLS = len(lls)
    lterm = nLLS * np.log(lvec)

    #;; Make the grids and add
    lgrid = np.outer(lterm, np.ones(ngrid))
    agrid = np.outer(np.ones(ngrid), asum_zi)
    pos_grid = agrid + lgrid

    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    #;; Negative term  :: Quasars
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    #;; Alpha
    nstep = len(gz)
    dzstep = z_gz[1]-z_gz[0]

    agrid = np.zeros((nstep, ngrid))
    f_nrm = (1+z_gz)/(1+zpivot)
    for kk in range(ngrid):
        agrid[:,kk] = f_nrm**avec[kk]

    #;; gofz
    ggrid = np.outer(gz, np.ones(ngrid))
    atot = dzstep * np.sum( agrid * ggrid, axis=0)

    #; New grids
    agrid = np.outer(np.ones(ngrid), atot)
    lgrid = np.outer(lvec, np.ones(ngrid))
    neg_grid = agrid * lgrid

    #;;;;;;;;;;;;;;;;;;;;;;;;;;
    #;;  Likelihood
    lik = pos_grid - neg_grid
    maxL = np.max(lik)
    nrm_lik = lik - maxL

    # Write to disk
    if True:
        from astropy.io import fits
        hdu = fits.PrimaryHDU(nrm_lik)
        hdul = fits.HDUList([hdu])
        hdul.writeto('nrm_lik.fits', overwrite=True)

    # Unravel
    ilmx, iamx = np.unravel_index(nrm_lik.argmax(), nrm_lik.shape)
    print('Max: (l,a) = ', lvec[ilmx], avec[iamx])
    
    # Return
    return lik, lvec, avec


