""" Classes for the Circumgalactic Medium
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import warnings

from astropy import units as u

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

    TODO:
      Need to check cosmologies
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
                self.galaxy.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                self.galaxy.coord.dec.to_string(sep=':', pad=True, alwayssign=True),
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
      Defaults to WMAP9

    Attributes
    ----------
    rho : Quantity
      Impact parameter (u.kpc)
    """

    # Initialize
    def __init__(self, galaxy, igm_sys, cosmo=None, name=None):
        # Checks
        if not isinstance(galaxy, Galaxy):
            raise IOError('CGMAbsSys instantiated with a Galaxy')
        self.galaxy = galaxy
        if self.galaxy.z is None:
            raise IOError('Galaxy redshift *must* be specified')
        if not isinstance(igm_sys, IGMSystem):
            raise IOError('CGMAbsSys instantiated with an IGMSystem')
        self.igm_sys = igm_sys

        # Raise error for redshifts not in the Hubble flow
        if galaxy.z < 0.05:
            raise NotImplementedError("Not prepared for such low redshift.  Need to implement corrections.")

        # Calculate rho
        if cosmo is None:
            from astropy.cosmology import WMAP9 as cosmo
            warnings.warn('cgm.CGMAbsSys: Using WMAP9 cosmology')
            self.cosmo = cosmo
        else:
            self.cosmo = cosmo
        ang_sep = self.igm_sys.coord.separation(self.galaxy.coord).to('arcmin')
        kpc_amin = cosmo.kpc_comoving_per_arcmin(self.galaxy.z)  # kpc per arcmin
        self.rho = ang_sep * kpc_amin / (1+self.galaxy.z)  # Physical
        self.ang_sep = ang_sep.to('arcsec')

        # Calculate PA too?
        self.PA = self.igm_sys.coord.position_angle(self.galaxy.coord)

        # Standard name
        if name is None:
            self.name = 'J{:s}{:s}_{:d}_{:d}'.format(
                    self.igm_sys.coord.ra.to_string(unit=u.hour,sep='',pad=True)[0:4],
                    self.igm_sys.coord.dec.to_string(sep='',pad=True,alwayssign=True)[0:5],
                    int(np.round(self.PA.to('deg').value)),
                    int(np.round(self.ang_sep.to('arcsec').value)))
        else:
            self.name = name

    # Output
    def __repr__(self):
        return ('<{:s}: {:s} Galaxy RA/DEC={:s}{:s}, zgal={:g}, rho={:g}]'.format(
                self.__class__.__name__,
                self.name,
                 self.galaxy.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.galaxy.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                 self.galaxy.z, self.rho))