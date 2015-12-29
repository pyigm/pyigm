""" Classes for the Circumgalactic Medium
"""
from __future__ import print_function, absolute_import, division, unicode_literals

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
            raise IOError("CGM instantiated with a Galaxy")
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
    def __init__(self, galaxy, igm_sys, cosmo=None):
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
            from astropy.cosmology import WMAP9 as cosmo
            print('cgm.CGMAbsSys: Using WMAP9 cosmology')
            self.cosmo = cosmo
        else:
            self.cosmo = cosmo
        ang_sep = self.igm_sys.coord.separation(self.galaxy.coord).to('arcmin')
        kpc_amin = cosmo.kpc_comoving_per_arcmin( self.galaxy.z)  # kpc per arcmin
        self.rho = ang_sep * kpc_amin / (1+self.galaxy.z)  # Physical
        # Calculate PA too?

    # Output
    def __repr__(self):
        return ('<{:s}: Galaxy RA/DEC={:s}{:s}, zgal={:g}, rho={:g}]'.format(
                self.__class__.__name__,
                 self.galaxy.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.galaxy.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                 self.galaxy.z, self.rho))
