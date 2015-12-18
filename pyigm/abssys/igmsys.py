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
    def __init__(self, abstype, radec, zabs, vlim, ZH=0., **kwargs):
        """Standard init
        """
        # Generate with type
        AbsSystem.__init__(self, abstype, radec, zabs, vlim, **kwargs)
        # Init
        self.ZH = ZH

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


class HISystem(AbsSystem):
    """Class for HI Lyman Absorption Line System
    """
    def __init__(self, radec, zabs, vlim, **kwargs):
        IGMSystem.__init__(self, 'HI', radec, zabs, vlim, **kwargs)

    def chk_component(self,component):
        """Require components are only of HI
        """
        # Require HI
        test = (component.Zion[0] == 1) & (component.Zion[1] == 1)
        if not test:
            warnings.warn('Input component must be HI')
        return test

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'HI'

class AbsSubSystem(object):
    """
    Sub system.  Most frequently used in LLS

    Used to analyze a portion of an AbsSystem
    Must be tied to an AbsSystem

    Parameters
    ----------
    parent : AbsSystem
      Link
    zabs : float
      Absorption redshift
    vlim : Quantity array (2)
      Velocity limits of the system
    lbl : str
      Label for the SubSystem, e.g. 'A', 'B'
    """
    def __init__(self, parent, zabs, vlim, lbl):
        self.parent = parent
        self.zabs = zabs
        self.vlim = vlim
        self.lbl = lbl
        self._ionN = None

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'SubSystem'

        # #############
    def __repr__(self):
        txt = '[{:s}: name={:s}{:s} type={:s}, {:s} {:s}, z={:g}, vlim={:g},{:g}'.format(
            self.__class__.__name__, self.parent.name, self.lbl,
            self.parent.abs_type,
            self.parent.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
            self.parent.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
            self.zabs, self.vlim[0],self.vlim[1])
        # Finish
        txt = txt + ']'
        return (txt)