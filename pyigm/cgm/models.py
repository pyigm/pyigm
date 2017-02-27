"""  Module for CGM models
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from astropy import units as u
from astropy import constants as const

try:
    basestring
except NameError:  # For Python 3
    basestring = str


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
      'hot'

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
    M_halo : float, optional
      log10 of the Halo mass (solar masses)
    """
    def __init__(self, M_halo=12.2, c=7.67, **kwargs):
        # Init
        CGMPhase.__init__(self, 'hot')
        # Param
        self.M_halo = M_halo
        self.setup_param()

    def setup_param(self):
        """ Setup key parameters of the model
        """
        # Cosmology
        self.H0 = 70. *u.km/u.s/ u.Mpc
        self.fb = 0.16       # Baryon fraction
        self.rhoc = 9.2e-30 * u.g / u.cm**3
        # DM
        self.r200 = (((3*10**self.M_halo * const.M_sun.cgs) / (4*np.pi*200*self.rhoc))**(1/3)).to('kpc')
        pdb.set_trace()

    def fy(self, y):
        """ Enclosed mass function
        Parameters
        ----------
        y : float or ndarray

        Returns
        -------
        f_y : float or ndarray

        """
        f_y = 


