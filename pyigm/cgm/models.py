"""  Module for CGM models
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from scipy.special import hyp2f1
from scipy.interpolate import interp1d

from pkg_resources import resource_filename

from astropy import units as u
from astropy import constants as const
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.utils import isiterable

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
        if not isinstance(key, str):
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
      Will likely use this for all diffuse gas
    alpha : float, optional
      Parameter to modify NFW profile power-law
    y0 : float, optional
      Parameter to modify NFW profile position
    """
    def __init__(self, log_Mhalo=12.2, c=7.67, f_hot=0.75, alpha=0., y0=1., **kwargs):
        # Init
        CGMPhase.__init__(self, 'hot')
        # Param
        self.log_Mhalo = log_Mhalo
        self.M_halo = 10.**self.log_Mhalo * const.M_sun.cgs
        self.c = c
        self.alpha = alpha
        self.y0 = y0
        self.f_hot = f_hot
        self.zero_inner_ne = 0. # kpc

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
        if self.zero_inner_ne > 0.:
            if np.sum(xyz**2) < self.zero_inner_ne**2:
                if isiterable(ne):
                    return np.zeros_like(ne)
                else:
                    return 0.
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
        nH = (self.rho_b(xyz) / self.mu / m_p).cgs.value
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

    def Ne_Rperp(self, Rperp, step_size=0.1*u.kpc, rmax=1.):
        """ Calculate N_H at an input impact parameter Rperp
        Just a simple sum in steps of step_size

        Parameters
        ----------
        Rperp : Quantity
          Impact parameter, typically in kpc
        step_size : Quantity
          Step size used for numerical integration (sum)
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
        zmax = np.sqrt((rmax*self.r200) ** 2 - Rperp ** 2).to('kpc')
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

class MB04(ModifiedNFW):
    def __init__(self, Rc=167*u.kpc, log_Mhalo=12.2, c=7.67, f_hot=0.75, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot, **kwargs)

        # Setup
        self.Rs = self.r200/self.c
        self.Rc = Rc
        self.Cc = (self.Rc/self.Rs).decompose().value
        self.rhoV = 1. * const.m_p/u.cm**3  # Will be renormalized

        # For development
        self.debug=False

        # Normalize
        self.norm_rhoV()

    def norm_rhoV(self):
        # Set rhoV to match expected baryon mass
        r = np.linspace(1., self.r200.to('kpc').value, 1000)  # kpc
        # Set xyz
        xyz = np.zeros((3,r.size))
        xyz[2, :] = r
        #
        dr = r[1] - r[0]
        Mass_unnorm = 4 * np.pi * np.sum(r**2 * self.rho_b(xyz)) * dr * u.kpc**3 # g * kpc**3 / cm**3
        # Ratio
        rtio = (Mass_unnorm/self.M_b).decompose().value
        self.rhoV = self.rhoV.cgs/rtio
        #
        print("rhoV normalized to {} to give M_b={}".format((self.rhoV/const.m_p).cgs,
                                                            self.M_b.to('Msun')))

    def rho_b(self, xyz):
        radius = np.sqrt(rad3d2(xyz))
        x = radius/self.Rs.to('kpc').value
        #
        rho = self.rhoV * (1+ (3.7/x)*np.log(1+x) - (3.7/self.Cc) * np.log(1+self.Cc))**(3/2)
        if self.debug:
            pdb.set_trace()
        #
        return rho


class YF17(ModifiedNFW):
    """
    Y. Faermen (2017) model of the Milky Way

    """
    def __init__(self, log_Mhalo=12.18, c=7.67, f_hot=0.75, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot, **kwargs)

        # Read
        faerman_file = resource_filename('pyigm', '/data/CGM/Models/Faerman_2017_ApJ_835_52-density-full.txt')
        self.yf17 = Table.read(faerman_file, format='ascii.cds')
        self.yf17['nH'] = self.yf17['nHhot'] + self.yf17['nHwarm']

        # For development
        self.debug=False

        # Setup
        self.rhoN = const.m_p/u.cm**3
        self.setup()

    def setup(self):
        # Setup Interpolation
        self.yf17_interp = interp1d(self.yf17['Radius'], self.yf17['nH'], kind='cubic', bounds_error=False, fill_value=0.)

        # Set rhoN to match expected baryon mass
        r = np.linspace(1., self.r200.to('kpc').value, 1000)  # kpc
        # Set xyz
        xyz = np.zeros((3,r.size))
        xyz[2, :] = r
        #
        dr = r[1] - r[0]
        Mass_unnorm = 4 * np.pi * np.sum(r**2 * self.rho_b(xyz)) * dr * u.kpc**3 # g * kpc**3 / cm**3
        # Ratio
        rtio = (Mass_unnorm/self.M_b).decompose().value
        self.rhoN = self.rhoN.cgs/rtio
        #
        print("rhoN normalized to {} to give M_b={}".format((self.rhoN/const.m_p).cgs,
                                                            self.M_b.to('Msun')))

    def rho_b(self, xyz):
        radius = np.sqrt(rad3d2(xyz))
        #
        rho = self.rhoN * self.yf17_interp(radius)
        if self.debug:
            pdb.set_trace()
        #
        return rho


class MilkyWay(ModifiedNFW):
    """
    JXP preferred model for the Galaxy

    """
    def __init__(self, log_Mhalo=12.18, c=7.67, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)


class M31(ModifiedNFW):
    """
    JXP preferred model for M31

    Taking mass from van der Marel 2012

    """
    def __init__(self, log_Mhalo=12.18, c=7.67, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)
        # Position from Sun
        self.distance = 752 * u.kpc # (Riess, A.G., Fliri, J., & Valls - Gabaud, D. 2012, ApJ, 745, 156)
        self.coord = SkyCoord('J004244.3+411609', unit=(u.hourangle, u.deg),
                              distance=self.distance)


class LMC(ModifiedNFW):
    """
    Preferred model for LMC

    Taking data from D'Onghia & Fox ARAA 2016

    """
    def __init__(self, log_Mhalo=np.log10(1.7e10), c=7.67, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)
        # Position from Sun
        self.distance = 50 * u.kpc
        self.coord = SkyCoord('J052334.6-694522', unit=(u.hourangle, u.deg),
                              distance=self.distance)

class SMC(ModifiedNFW):
    """
    Preferred model for SMC

    Taking data from D'Onghia & Fox ARAA 2016

    """
    def __init__(self, log_Mhalo=np.log10(2.4e9), c=7.67, f_hot=0.75, alpha=2, y0=2, **kwargs):

        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot,
                             alpha=alpha, y0=y0, **kwargs)
        # Position from Sun
        self.distance = 61 * u.kpc
        self.coord = SkyCoord('J005238.0-724801', unit=(u.hourangle, u.deg),
                              distance=self.distance)

class ICM(ModifiedNFW):
    """
    ICM model following Vikhilnin et al. 2006

    """
    def __init__(self, log_Mhalo=np.log10(5e14), c=3.5, f_hot=0.75, **kwargs):
        # Init ModifiedNFW
        ModifiedNFW.__init__(self, log_Mhalo=log_Mhalo, c=c, f_hot=f_hot, **kwargs)

        # Using their values for A907
        self.n0 = 6.252e-3 / u.cm**3
        self.rc = 136.9 #* u.kpc
        self.rs = 1887.1 #* u.kpc
        self.alpha = 1.556
        self.beta = 0.594
        self.epsilon = 4.998
        self.n02 = 0.

        # Fixed
        self.gamma = 3

    def ne(self, xyz):
        """

        Parameters
        ----------
        xyz : ndarray
          Coordinate(s) in kpc

        Returns
        -------
        n_e : float or ndarray
          electron density in cm**-3

        """
        radius = np.sqrt(rad3d2(xyz))

        # This ignores the n02 term
        npne = self.n0**2 * (radius/self.rc)**(-self.alpha) / (
                (1+(radius/self.rc)**2)**(3*self.beta - self.alpha/2.)) * (1 /
            (1+(radius/self.rs)**self.gamma)**(self.epsilon/self.gamma))
        if self.n02 > 0:
            pdb.set_trace()  # I didn't code this yet

        ne = self.nH(xyz) * 1.1667
        if self.zero_inner_ne > 0.:
            if np.sum(xyz**2) < self.zero_inner_ne**2:
                if isiterable(ne):
                    return np.zeros_like(ne)
                else:
                    return 0.
        # Return
        return ne
