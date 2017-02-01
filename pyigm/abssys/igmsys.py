""" Subclass of AbsSystem for IGM Systems
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import warnings

from astropy import units as u

from linetools.isgm.abssystem import AbsSystem
from linetools import utils as ltu
from linetools.isgm import utils as ltiu

class IGMSystem(AbsSystem):
    """
    Class for an IGM absorption system

    Parameters
    ----------
    radec : tuple or coordinate
        RA/Dec of the sightline or astropy.coordinate
    zabs : float
      Absorption redshift
    vlim : Quantity array (2)
      Velocity limits of the system
      Defaulted to +/- 500 km/s if None
    abstype : str, optional
      Type of IGM System, e.g. LLS, DLA, MgII
    **kwargs : keywords
      passed to AbsSystem.__init__
    """
    def __init__(self, radec, zabs, vlim, ZH=0., abs_type='IGMSystem', **kwargs):
        """Standard init
        """
        if vlim is None:
            vlim = [-500.,500.]*u.km/u.s
        # Generate with type
        AbsSystem.__init__(self, radec, zabs, vlim, abs_type=abs_type, **kwargs)
        # Init
        self.ZH = ZH

    # Output
    def __repr__(self):
        return ('<{:s}: {:s} {:s} {:s}, {:g}, NHI={:g}, Z/H={:g}>'.format(
                self.__class__.__name__, self.abs_type,
                 self.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                 self.coord.dec.to_string(sep=':', pad=True),
                 self.zabs, self.NHI, self.ZH))


class HISystem(IGMSystem):
    """Class for HI Lyman Absorption Line System
    """
    def __init__(self, radec, zabs, vlim, **kwargs):
        IGMSystem.__init__(self, radec, zabs, vlim, abs_type='HI', **kwargs)

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

    @classmethod
    def from_dict(cls, parent, idict, lbl):
        """ Generate a sub-system from a dict

        Parameters
        ----------
        parent : AbsSystem
          Link
        idict : dict
          Contains the sub-system parameters
        lbl : str
        """
        slf = cls(parent, idict['zabs'], idict['vlim']*u.km/u.s, lbl)
        # Components
        components = ltiu.build_components_from_dict(idict)
        slf._components = components
        # Ion table
        slf._ionN = ltiu.iontable_from_components(components)
        # Return
        return slf

    def __init__(self, parent, zabs, vlim, lbl):
        self.parent = parent
        self.zabs = zabs
        self.vlim = vlim
        self.lbl = lbl
        try:
            self.name = parent.name+'_'+self.lbl
        except AttributeError:
            self.name = ''
        self._ionN = None

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'SubSystem'

    def to_dict(self):
        """ Generate a dict from the sub-system

        Returns
        -------
        outdict : dict
          JSON capatible

        """
        import datetime
        import getpass
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        user = getpass.getuser()
        # Generate the dict
        outdict = dict(abs_type='SubSystem', Name=self.name, zabs=self.zabs,
                       vlim=self.vlim.to('km/s').value,
                       lbl=self.lbl,
                       CreationDate=date,
                       user=user
                       )
        # Components
        outdict['components'] = {}
        for component in self._components:
            outdict['components'][component.name] = component.to_dict()
        # Polish
        outdict = ltu.jsonify(outdict)
        # Return
        return outdict

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