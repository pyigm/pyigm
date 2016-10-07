""" module for analysis of CGM outside the CLasses
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import warnings
import pdb

from scipy.special import gamma, gammainc
try:
    import mpmath
except ImportError:
    warnings.warn("I hope your exponent is positive..")
    from scipy.special import gammainc
    gflg=False
else:
    from mpmath import gammainc
    gflg=True

from astropy import units as u
from astropy import constants as const
from astropy.cosmology import Planck15

def dndx_rvir(Lrng=(0.001, 10), nL=1000, beta=0.2, rvir_Lstar=250.*u.kpc,
         phi_str_pref = 1.49, alpha = -1.05, Mstar = -20.44,
         cosmo=None):  #  ; M* - 5 log h
    """ Estimate dN/dX for a set of CGM assuming unit covering fraction
    Following Prochaska+11

    Use beta=0 and rvir_Lstar=300kpc for a constant CGM to 300kpc

    Parameters
    ----------
    Lrng : tuple, optional
      Range to integrate luminosity in terms of L*
    nL : int, optional
      Number of evaluations of L in that interval
    beta : float, optional
      Parameterization of rvir with L
       r_vir = 250 kpc * (L/L*)^beta
    phi_str, alpha, Mstar : float, float, float
        Blanton lum function
          Phi(M) = 0.4 log(10) Phi* 10^(-0.4 [M-M*][alpha+1]) exp(-
                    ; 10^(-0.4[M-M*] ) )
          Phi(L) = Phi* (L/L*)^alpha exp(-L/L*)
          Phi* has the funny units of 1e-2 h^3 Mpc^-3
    cosmo : Cosmology, optional
      Defaults to Planck15

    Returns
    -------
    Lval : ndarray
      Luminosities of evaluation
    dNdX : float
      Cumulative dNdX

    """
    # Cosmology
    if cosmo is None:
        cosmo = Planck15
    hubb = cosmo.H0.value / 100.
    # Constants
    phi_str_cgs = (phi_str_pref * 1e-2 * hubb**3) / u.Mpc**3
    dndx_const = (const.c / cosmo.H0).cgs
    # Cumulative
    Lval = np.linspace(Lrng[0], Lrng[1], nL)
    x = alpha + 1 + beta
    # Integrate
    if gflg:
        igmma = np.zeros_like(Lval)
        i0 = 1.-float(gammainc(x,Lrng[1], regularized=True))
        for kk,iLval in enumerate(Lval):
            igmma[kk] = i0 - (1-float(gammainc(x,iLval, regularized=True)))
    else:
        igmma = gammainc(x,Lrng[1]) - gammainc(x,Lval)
    #from scipy.special import gammainc as gic
    #ig2 = gic(x,Lrng[1]) - gic(x,Lval)
    #pdb.set_trace()
    dNdx_rvir = (dndx_const * phi_str_cgs * (np.pi * rvir_Lstar**2) * (
        gamma(x) * igmma)).decompose().value

    # Return
    return Lval, dNdx_rvir

