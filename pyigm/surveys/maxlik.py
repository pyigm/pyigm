""" methods for Maximum Likelihood analyses on IGM surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb

import numpy as np

from scipy.interpolate import interp1d
from scipy.stats import kstest


def powerlaw_loxz(lls_zabs, z_gz, gz, guess, zpivot, rnga=(-0.5,0.5), rngl=(-0.2,0.1), ngrid=100):
    """

    Ported from SDSS LLS paper
    Functional form  l = k (1+z/1+zpivot)^a

    Parameters
    ----------
    lls_zabs
    z_gz
    gz
    guess
    zpivot
    rnga
    rngl
    ngrid

    Returns
    -------
    lik : ndarray
     Normmalized likelihood function
    lvec : ndarray
     l* values for the grid
    avec : ndarray
     alpha values for the grid

    """

    # Create vectors of alpha and l_0 values
    lvec = guess[0] * 10.**(rngl[0] + (rngl[1]-rngl[0])*np.arange(ngrid)/(ngrid-1))
    avec = guess[1] * 10.**(rnga[0] + (rnga[1]-rnga[0])*np.arange(ngrid)/(ngrid-1))

    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    #;; Positive term, i.e. h(z) evaluated at z_i
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    #;; Evaluate z_i and sum
    zi_nrm = (1+lls_zabs)/(1+zpivot)
    sum_lnzi = np.sum(np.log(zi_nrm))

    #;; Weight by alpha
    asum_zi = avec * sum_lnzi

    #;; lvec sum is trivial
    nLLS = len(lls_zabs)
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


def cl_indices(lnL, cl, sigma=False):
    """ Find the indices of a log-Likelihood grid encompassing a
    given confidence interval

    Parameters:
      lnL: np.array
        log-Likelihood image
      sigma: bool, optional
        Return as sigma values [not implemented]

    Returns:
      indices: Tuple of np.where output
    """
    # Max
    mxL = np.max(lnL)

    # Noramlize and flatten
    norm_img = lnL-mxL
    flat_img = norm_img.flatten()

    # Sort
    srt = np.argsort(flat_img)
    norm_lnL = flat_img[srt]

    # Sum
    cumulsum = np.cumsum(np.exp(np.maximum(norm_lnL,-15.)))
    cumul = cumulsum/cumulsum[-1]

    # Interpolation (smoothing a bit)
    fsum = interp1d(norm_lnL, cumul)

    # Finish
    indices = np.where(fsum(norm_img) > (1-cl))

    # Return
    return indices

def powerlaw_ks(lls_zabs, z_gz, gz, lstar, alpha, zpivot):
    """
    Perform a KS test on the LLS zabs distribution vs. that
    predicted from the l(z) power-law model

    Parameters
    ----------
    lls_zabs
    z_gz
    gz
    lstar
    alpha
    zpivot

    Returns
    -------
    D : float
    p_value : float

    """
    # Generate model zabs distribution
    lz = lstar * ((1+z_gz)/(1+zpivot))**alpha

    # Evaluate with g(z)
    dz = z_gz[1]-z_gz[0]
    nLLS_z = lz * gz * dz

    # Build CDF
    cdf = np.cumsum(nLLS_z)
    cdf /= cdf[-1]
    cdf_interp = interp1d(z_gz, cdf)
    def wrap_cdf(z):
        return cdf_interp(z)

    # KS Test
    D, p_value = kstest(lls_zabs, wrap_cdf)

    # Return
    return D, p_value



