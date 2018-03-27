"""  Module for CGM models
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
from scipy.special import hyp2f1

from astropy import units as u
from astropy import constants as const

try:
    basestring
except NameError:  # For Python 3
    basestring = str

def rad3d2(xyz):
    """ Calculate radius to x,y,z inputted
    Assumes the origin is 0,0,0
    Parameters
    ----------
    xyz : Tuple or ndarray

    Returns
    -------
    rad3d : float or ndarray

    """
    return xyz[0]**2 + xyz[1]**2 + xyz[-1]**2

# Constants
m_p = const.m_p.cgs.value # g

class CGMModel(object):
    """Model of the CGM
    Wraps multiple phases together

    Parameters:
    -----------
    """
    def __init__(self, **kwargs):
        # Groups
        self._pdict = {}   # Dict for the phases

    def __getitem__(self, key):
        """ Access the DB groups

        Parameters
        ----------
        key : str

        Returns
        -------

        """
        # Check
        if not isinstance(key, basestring):
            raise IOError("Item must be str")
        # Try to access the dict
        try:
            return self._pdict[key]
        except KeyError:
            raise IOError("Input phase={:s} is not loaded in the model".format(key))

    def __repr__(self):
        txt = '<{:s}: '.format(self.__class__.__name__)
        # Phases
        txt += '   Phases modeled = {} \n'.format(self._pdict.keys())
        txt += '>'
        return (txt)

class CGMPhase(object):
    """ Model of a single phase of the CGM
    Parameters:
    -----------
    phase : str
      'hot' : T > 10^6 K

    """
    def __init__(self, phase, **kwargs):
        # Check
        assert phase in ['hot']
        #
        self.phase = phase
        self.mass = 0.   # Total mass in Solar masses
        self.r_cgm = 0.  # Radius of the CGM

    def nH(self, xyz):
        """
        Parameters
        ----------
        xyz : ndarray
          Coordinate(s) in kpc

        Returns
        -------
        nH : float or ndarray
          Hydrogen number density in cm**-3
        """
        pass


class ModifiedNFW(CGMPhase):
    """ Generate a modified NFW model, e.g. Mathews & Prochaska 2017
    for the hot gas.  Currently only valid for z=0

    Parameters:
    -----------
    log_Mhalo : float, optional
      log10 of the Halo mass (solar masses)
    c : float, optional
      concentration of the halo
    f_hot : float, optional
      Fraction of the baryons in this hot phase
    alpha : float, optional
      Parameter to modify NFW profile power-law
    y0 : float, optional
      Parameter to modify NFW profile position
    """
    def __init__(self, log_Mhalo=12.2, c=7.67, f_hot=0.9, alpha=0., y0=1., **kwargs):
        # Init
        CGMPhase.__init__(self, 'hot')
        # Param
        self.log_Mhalo = log_Mhalo
        self.M_halo = 10.**self.log_Mhalo * const.M_sun.cgs
        self.c = c
        self.alpha = alpha
        self.y0 = y0
        self.f_hot = f_hot
        # Init more
        self.setup_param()

    def setup_param(self):
        """ Setup key parameters of the model
        """
        # Cosmology
        self.H0 = 70. *u.km/u.s/ u.Mpc
        self.fb = 0.16       # Baryon fraction
        self.rhoc = 9.2e-30 * u.g / u.cm**3
        # Dark Matter
        self.r200 = (((3*self.M_halo) / (4*np.pi*200*self.rhoc))**(1/3)).to('kpc')
        self.rho0 = 200*self.rhoc/3 * self.c**3 / self.fy_DM(self.c)   # Central density
        # Baryons
        self.M_b = self.M_halo * self.fb
        self.rho0_b = (self.M_b / (4*np.pi) * (self.c/self.r200)**3 / self.fy_b(self.c)).cgs
        # Misc
        self.mu = 1.33   # Reduced mass correction for Helium

    def fy_DM(self, y):
        """ Enclosed mass function for the DM
        NFW

        Parameters
        ----------
        y : float or ndarray

        Returns
        -------
        f_y : float or ndarray
        """
        f_y = np.log(1+y) - y/(1+y)
        #
        return f_y

    def fy_b(self, y):
        """ Enclosed mass function for the baryons

        Parameters
        ----------
        y : float or ndarray

        Returns
        -------
        f_y : float or ndarray
        """
        f_y = (y/(self.y0 + y))**(1+self.alpha) * (
            self.y0**(-self.alpha) * (self.y0 + y)**(1+self.alpha) * hyp2f1(
                1+self.alpha, 1+self.alpha, 2+self.alpha, -1*y/self.y0)
            - self.y0) / (1+self.alpha) / self.y0
        return f_y

    def ne(self, xyz):
        """ Calculate n_e from n_H with a correction for Helium
        Assume 25% mass is Helium and both electrons are off

        Parameters
        ----------
        xyz : ndarray
          Coordinate(s) in kpc

        Returns
        -------
        n_e : float or ndarray
          electron density in cm**-3

        """
        ne = self.nH(xyz) * 1.1667
        # Return
        return ne

    def nH(self, xyz):
        """ Calculate the Hydrogen number density
        Includes a correction for Helium

        Parameters
        ----------
        xyz : ndarray
          Coordinate(s) in kpc

        Returns
        -------
        nH : float or ndarray
          Density in cm**-3

        """
        nH = self.rho_b(xyz).value / self.mu / m_p
        # Return
        return nH

    def rho_b(self, xyz):
        """ Mass density in baryons; modified

        Parameters
        ----------
        xyz : ndarray
          Position assumed in kpc

        Returns
        -------
        rho : Quantity
          Density in g / cm**-3

        """
        radius = np.sqrt(rad3d2(xyz))
        y = self.c * (radius/self.r200.to('kpc').value)
        rho = self.rho0_b / y**(1-self.alpha) / (self.y0+y)**(2+self.alpha)
        # Return
        return rho

    def Ne_Rperp(self, Rperp, step_size=0.1*u.kpc, rmax=1., epsrel=1e-4, epsabs=1e-6,
                 *arg, **kwargs):
        """ Calculate N_H at an input impact parameter Rperp

        Parameters
        ----------
        Rperp : Quantity
          Impact parameter, typically in kpc
        step_size : Quantity
          Step size used for numerical integration
        rmax : float
          Maximum radius for integration in units of r200

        Returns
        -------
        Ne : Quantity
          Column density of total electrons
        """
        dz = step_size.to('kpc').value

        # Cut at rmax*rvir
        if Rperp > rmax*self.r200:
            return 0. / u.cm**2
        # Generate a sightline to rvir
        zmax = np.sqrt(self.r200 ** 2 - Rperp ** 2).to('kpc')
        zval = np.arange(-zmax.value, zmax.value+dz, dz)  # kpc
        # Set xyz
        xyz = np.zeros((3,zval.size))
        xyz[0, :] = Rperp.to('kpc').value
        xyz[2, :] = zval

        # Integrate
        ne = self.ne(xyz) # cm**-3
        Ne = np.sum(ne) * dz * 1000  # pc cm**-3

        # Return
        return Ne * u.pc / u.cm**3
