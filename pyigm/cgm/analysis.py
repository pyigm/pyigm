""" module for analysis of CGM outside the CLasses
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import warnings
import pdb

from scipy.special import gamma, gammainc

from astropy import units as u
from astropy import constants as const
from astropy.cosmology import Planck15

def dndx_rvir(Lrng=(0.001, 10), beta=0.2, rvir_Lstar=250.*u.kpc,
         phi_str = 1.49/u.Mpc**3, alpha = -1.05, Mstar = -20.44,
         cosmo=None):  #  ; M* - 5 log h
    """ Estimate dN/dX for a set of CGM assuming unit covering fraction
    Use beta=0 and rvir_Lstar=300kpc for a constant CGM to 300kpc

    Parameters
    ----------
    Lrng : tuple, optional
      Range to integrate luminosity in terms of L*
    beta : float, optional
      Parameterization of rvir with L
       r_vir = 250 kpc * (L/L*)^beta
    phi_str, alpha, Mstar : Quantity, float, float
        Blanton lum function
          Phi(M) = 0.4 log(10) Phi* 10^(-0.4 [M-M*][alpha+1]) exp(-
                    ; 10^(-0.4[M-M*] ) )
          Phi(L) = Phi* (L/L*)^alpha exp(-L/L*)
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
    phi_str_cgs = (phi_str * 1e-2 * hubb)
    dndx_const = (const.c / cosmo.H0).cgs
    # CGM extends to rvir with power-law dependence
    Lval = np.linspace(Lrng[0], Lrng[1], 1000)
    x = alpha + 1 + beta
    # Integrate
    dNdx_rvir = (dndx_const * phi_str_cgs * (np.pi * rvir_Lstar**2) * (
        gamma(x) * ( gammainc(x,Lrng[1]) - gammainc(x,Lval)))).decompose()

    # Return
    return Lval, dNdx_rvir

