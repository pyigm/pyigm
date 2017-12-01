""" Classes for the Circumgalactic Medium
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import warnings
import pdb

from astropy import units as u
from astropy import cosmology

from linetools import utils as ltu

from pyigm.field.galaxy import Galaxy
from pyigm.abssys.igmsys import IGMSystem

class CGM(object):
    """A CGM Class

    Parameters
    ----------
    galaxy : Galaxy
      Must include redshift

    Attributes
    ----------
    phases : dict
        dict describing the so-called "phases" of the CGM
    rlim : Quantity array (2)
      Radial limits
    cgm_abs : list
      List of CGMAbsSys classes

    """
    # Initialize with a .dat file
    def __init__(self, galaxy):
        # Checks
        if not isinstance(galaxy, Galaxy):
            raise IOError("CGM must be instantiated with a Galaxy")
        self.galaxy = galaxy
        if self.galaxy.z is None:
            raise IOError("Galaxy redshift *must* be specified")
        # Phases
        self.phases = dict(total={},
                           cold=dict(T=[0, 1e3]*u.K),
                           cool=dict(T=[1e3, 1e5]*u.K),
                           warm=dict(T=[1e5, 1e6]*u.K),
                           hot=dict(T=[1e6, 1e10]*u.K))
        # rlim
        self.rlim = None
        # AbsSys
        self.cgm_abs = []

    def __repr__(self):
        return ('<CGM: {:s} {:s}, z={:g}>'.format(
                self.galaxy.coord.icrs.ra.to_string(unit=u.hour, sep=':', pad=True),
                self.galaxy.coord.icrs.dec.to_string(sep=':', pad=True, alwayssign=True),
                self.galaxy.z))


class CGMAbsSys(object):
    """ Class for a CGM Absorption system

    Combines an AbsSystem with a Galaxy

    Parameters
    ----------
    galaxy : Galaxy
      Must include redshift
    igm_sys : IGMSystem
    cosmo : astropy.cosmology, optional
      Defaults to Planck15

    Attributes
    ----------
    rho : Quantity
      Impact parameter (u.kpc)
    """
    @classmethod
    def from_dict(cls, idict, use_angrho=False, **kwargs):
        """ Generate a CGMAbsSys object from a dict

        Parameters
        ----------
        idict : dict
        use_angrho : bool, optional
          Use ang_sep and rho if in idict and cosmo is too
        """
        # Galaxy object
        galaxy = Galaxy.from_dict(idict['galaxy'])
        # IGM system
        igm_sys = IGMSystem.from_dict(idict['igm_sys'], **kwargs)
        # Keywords
        kwargs2 = kwargs.copy()
        kwargs2['name'] = idict['Name']
        if 'cosmo' in idict.keys():
            kwargs2['cosmo'] = getattr(cosmology, idict['cosmo'])
            if use_angrho:
                kwargs2['ang_sep'] = idict['ang_sep']*u.arcsec
                kwargs2['rho'] = idict['rho']*u.kpc
        # Instantiate
        slf = cls(galaxy, igm_sys, **kwargs2)
        # Extras
        for key in ['ebv']:
            if key in idict.keys():
                setattr(slf, key, idict[key])
        # Return
        return slf

    @classmethod
    def from_json(cls, jfile, **kwargs):
        """
        Parameters
        ----------
        jfile : str
        """
        idict = ltu.loadjson(jfile)
        slf = cls.from_dict(idict, **kwargs)
        return slf

    def __init__(self, galaxy, igm_sys, cosmo=None, name=None, rho=None, PA=None,
                 ang_sep=None, correct_lowz=True, debug=False, **kwargs):
        """
        Parameters
        ----------
        galaxy : Galaxy
        igm_sys : IGMSystem
        cosmo : astropy.cosmology, optional
          Defaults to Planck15
        name : str, optional
        rho : Quantity, optional
          Impact parameter; calculated if not input
        PA : Quantity
          Position angle
        ang_sep : Quantity
          Angular separation between galaxy and sightline
        correct_lowz : bool, optional
          If galaxy z < 0.05, correct for peculiar velocies in impact parameter
          calculation

        Returns
        -------

        """
        from pyigm.cgm.utils import calc_cgm_rho
        # Checks
        if not isinstance(galaxy, Galaxy):
            raise IOError('CGMAbsSys instantiated with a Galaxy')
        self.galaxy = galaxy
        if self.galaxy.z is None:
            raise IOError('Galaxy redshift *must* be specified')
        if not isinstance(igm_sys, IGMSystem):
            raise IOError('CGMAbsSys instantiated with an IGMSystem')
        self.igm_sys = igm_sys

        # Calculate rho
        if cosmo is None:
            warnings.warn('cgm.CGMAbsSys: Using Planck15 cosmology')
            self.cosmo = cosmology.Planck15
        else:
            self.cosmo = cosmo

        # Impact parameter and PA
        if rho is None:
            rho, iang = calc_cgm_rho(galaxy, igm_sys, self.cosmo, ang_sep=ang_sep,
                                 correct_lowz=correct_lowz, **kwargs)
            if debug:
                pdb.set_trace()
            if ang_sep is None:
                ang_sep = iang
        self.rho = rho
        self.ang_sep = ang_sep
        self.PA = self.igm_sys.coord.position_angle(self.galaxy.coord)

        # Standard name
        if name is None:
            self.name = 'J{:s}{:s}_{:d}_{:d}'.format(
                    self.igm_sys.coord.icrs.ra.to_string(unit=u.hour,sep='',pad=True)[0:4],
                    self.igm_sys.coord.icrs.dec.to_string(sep='',pad=True,alwayssign=True)[0:5],
                    int(np.round(self.PA.to('deg').value)),
                    int(np.round(self.ang_sep.to('arcsec').value)))
        else:
            self.name = name

    def to_dict(self):
        """ Convert the system to a JSON-ready dict for output

        Returns
        -------
        cdict : dict
        """
        import datetime
        import getpass
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        user = getpass.getuser()
        # Generate the dict
        outdict = dict(Name=self.name, z=self.galaxy.z, rho=self.rho.value,
                       ang_sep=self.ang_sep.value,
                       PA=self.PA.value,
                       RA=self.galaxy.coord.icrs.ra.value,
                       DEC=self.galaxy.coord.icrs.dec.value,
                       cosmo = self.cosmo.name,
                       ebv = self.ebv,
                       CreationDate=date,
                       user=user
                       )
        # IGM_SYS
        outdict['igm_sys'] = self.igm_sys.to_dict()
        # Galaxy
        outdict['galaxy'] = self.galaxy.to_dict()
        # Polish
        outdict = ltu.jsonify(outdict)
        # Return
        return outdict

    def __getattr__(self, k):
        """ Extend Attributes

        Parameters
        ----------
        k : str
          Attribute

        Returns
        -------
          Attribute off main class then galaxy then igm_sys
        """
        # Galaxy?
        try:
            return getattr(self.galaxy, k)
        except AttributeError:
            # Try AbsLine_Sys last
            try:
                return getattr(self.igm_sys, k)
            except AttributeError:
                return None

    def __repr__(self):
        return ('<{:s}: {:s} Galaxy RA/DEC={:s}{:s}, zgal={:g}, rho={:g}>'.format(
                self.__class__.__name__,
                self.name,
                 self.galaxy.coord.icrs.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.galaxy.coord.icrs.dec.to_string(sep=':',pad=True,alwayssign=True),
                 self.galaxy.z, self.rho))

    def write_json(self, outfil=None):
        """ Generate a JSON file from a CGMAbsSys object 

        Parameters
        ----------
        outfil : str
          output file

        Returns
        -------

        """
        # Generate the dict
        odict = self.to_dict()
        # Write
        if outfil is None:
            outfil = self.name+'.json'
        ltu.savejson(outfil, odict, overwrite=True, easy_to_read=True)
        # Finish
        print("Wrote {:s} system to {:s} file".format(self.name, outfil))

