""" Module for HI IGM calculations related to reionization
"""
import pdb

import numpy as np

from astropy import constants as const
from astropy.cosmology import Planck15

from linetools.lists.linelist import LineList


def igm_damping(wave, zend, zsource, cosmo=Planck15):
    """  Calculates IGM 'damping' wing based on Miralda-Escude 1998 formalism
    However, there were errors with his Appendix..

    Assumes a fully neutral medium at the mean density

    Note: This approximation is only valid for the optical depth away from
    line center

    Parameters
    ----------
    wave : Quantity array
      Wavelengths
    zend : float
      End of reionization epoch
    zsource : float
      Redshift of source (e.g. quasar, GRB)

    Returns
    -------
    tau : np.array
      NaN are replaced by 99

    """
    # Constants
    HIlines= LineList('HI')
    Lya = HIlines['HI 1215']
    Ralpha = (Lya['A'] / 4. / np.pi / (const.c/Lya['wrest'])).decompose()  # 2.02e-8  # A_12 / 4 pi nu_12

    # Mean hydrogen density today
    Y = 0.25  # Helium mass fraction
    rho_b = cosmo.Ob0 * cosmo.critical_density0.cgs
    nH = rho_b * (1-Y) / const.m_p.cgs

    # Offsets
    wave0 = Lya['wrest']*(1+zsource)
    dwave = wave-wave0
    delta = (dwave / wave0).decompose()

    # tau0 -- Note that M98 has a number of typos including the speed of light!
    tau0 = (3 * Lya['wrest']**3 * Lya['A'] * nH * (1+zsource)**1.5 /
            (8 * np.pi * cosmo.H0)).decompose()

    """
    def intx(x):
        # From Wolfram
        integral = np.sqrt(x)*(10*x**4 + 18*x**3 + 42*x**2 + 210*x - 315
                           )/(35*(x-1)) - 4.5*np.log((1+np.sqrt(x))/(1-np.sqrt(x)))
        return integral
    """

    def integral_M98(x):
        inte = x**4.5/(1-x) + (9./7.)* x**3.5 +(9./5.)*x**2.5 + 3*x**1.5
        inte += 9*x**0.5 - 4.5*np.log((1+x**0.5)/(1-x**0.5))

        return inte

    def calc_intx(x1, x2):
        return integral_M98(x2) - integral_M98(x1)
        #return intx(x2) - intx(x1)

    # Calculate x1, x2
    x1 = (1.+zend)/((1.+zsource)*(1.+delta))
    x2 = 1./(1+delta)

    tau = (tau0*Ralpha/np.pi) * (1+delta)**1.5 * calc_intx(x1,x2)

    # Scrub the NaNs
    idx = np.isnan(tau)
    tau[idx] = 99.

    # Return
    return tau.value



