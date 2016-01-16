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

    def update_vlim(self, sub_system=None):
        """ Update vlim in the main or subsystems

        Parameters
        ----------
        sub_system : str, optional
          If provided, apply to given sub-system.  Only used in LLS
        """
        def get_vmnx(components):
            vmin,vmax = 9999., -9999.
            for component in components:
                vmin = min(vmin, component.vlim[0].value)
                vmax = max(vmax, component.vlim[1].value)
            return vmin,vmax

        # Sub-system?
        if sub_system is not None:
            components = self.subsys[sub_system]._components
            vmin, vmax = get_vmnx(components)
            self.subsys[sub_system].vlim = [vmin, vmax]*u.km/u.s
        else:
            components = self._components
            vmin, vmax = get_vmnx(components)
            self.vlim = [vmin, vmax]*u.km/u.s

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
        outdict = ltu.jsonify_dict(outdict)
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