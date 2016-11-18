""" utility methods for the Circumgalactic Medium
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import warnings
import pdb

from astropy import units as u
from astropy import cosmology
from astropy.coordinates import SkyCoord

from linetools import utils as ltu

from pyigm.field.galaxy import Galaxy
from pyigm.cgm.cgm import CGMAbsSys


def cgm_from_galaxy_igmsystems(galaxy, igmsystems, R_max=300*u.kpc, dv_max=400*u.km/u.s,
                               cosmo=None):
    """ Generate a list of CGMAbsSys objects given an input galaxy and a list of IGMSystems

    Parameters
    ----------
    galaxy : Galaxy
    igmsystems : list
      list of IGMSystems
    R_max : Quantity
      Maximum projected separation from sightline to galaxy
    dv_max
      Maximum velocity offset between system and galaxy

    Returns
    -------
    cgm_list : list
      list of CGM objects generated

    """
    # Cosmology
    if cosmo is None:
        cosmo = cosmology.Planck15
    # Separation
    ra = [igmsystem.coord.ra.value for igmsystem in igmsystems]
    dec = [igmsystem.coord.dec.value for igmsystem in igmsystems]
    igm_coords = SkyCoord(ra=ra, dec=dec, unit='deg')
    ang_sep = galaxy.coord.separation(igm_coords)

    # R
    kpc_amin = cosmo.kpc_comoving_per_arcmin(galaxy.z)  # kpc per arcmin
    rho = ang_sep.to('arcmin') * kpc_amin / (1+galaxy.z)  # Physical

    # dv
    igm_z = [igmsystem.zabs for igmsystem in igmsystems]
    dv = ltu.v_from_z(galaxy.z, igm_z)

    # Rules
    match = np.where((rho<R_max) and (np.abs(dv) < dv_max))[0]
    if len(match) == 0:
        print("No CGM objects match your rules")
        return []
    else:
        # Loop to generate
        cgm_list = []
        for imatch in match:
            cgm = CGMAbsSys(galaxy, igmsystems[imatch], cosmo=cosmo)
            cgm_list.append(cgm)

    # Return
    return cgm_list



