""" Subclass of AbsSystem for IGM Systems
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import pdb

from astropy import units as u

from linetools.isgm.abssystem import AbsSystem


class IGMSystem(AbsSystem):
    """
    Class for an IGM absorption system

    Parameters
    ----------
    abstype : str
      Type of IGM System, e.g. LLS, DLA, MgII
    radec : tuple or coordinate
        RA/Dec of the sightline or astropy.coordinate
    zabs : float
      Absorption redshift
    vlim : Quantity array (2)
      Velocity limits of the system
      Defaulted to +/- 500 km/s if None
    **kwargs : keywords
      passed to AbsSystem.__init__
    """
    def __init__(self, abstype, radec, zabs, vlim, **kwargs):
        """Standard init
        """
        # Generate with type
        AbsSystem.__init__(self, abstype, radec, zabs, vlim, **kwargs)

    # Output
    def __repr__(self):
        return ('<{:s}: {:s} {:s} {:s}, {:g}, NHI={:g}, Z/H={:g}>'.format(
                self.__class__.__name__, self.abs_type,
                 self.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                 self.coord.dec.to_string(sep=':', pad=True),
                 self.zabs, self.NHI, self.ZH))

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'DLA'


