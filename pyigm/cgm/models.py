"""  Module for CGM models
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb


from astropy import constants as const


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



