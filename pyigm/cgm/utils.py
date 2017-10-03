""" utility methods for the Circumgalactic Medium
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import warnings
import pdb

from astropy import units as u
from astropy import cosmology
from astropy.coordinates import SkyCoord
from astropy import constants

from linetools import utils as ltu

from pyigm.field.galaxy import Galaxy


def calc_rho(galaxy, igm_sys, cosmo, ang_sep=None, correct_lowz=True,
             Galactic=False, d_Sun=8.0*u.kpc):
    # Loop?
    if isinstance(igm_sys, list):
        rhos = []
        angs = []
        for iigm in igm_sys:
            irho, iang = calc_rho(galaxy, iigm, cosmo)
            rhos.append(irho.value)
            angs.append(iang)
        return np.array(rhos)*u.kpc, angs


def cgm_from_galaxy_igmsystems(galaxy, igmsystems, R_max=300*u.kpc, dv_max=400*u.km/u.s,
                               cosmo=None, **kwargs):
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
    from pyigm.cgm.cgm import CGMAbsSys
    # Cosmology
    if cosmo is None:
        cosmo = cosmology.Planck15

    # R
    rho, angles = calc_rho(galaxy, igmsystems, cosmo)

    # dv
    igm_z = np.array([igmsystem.zabs for igmsystem in igmsystems])
    dv = ltu.dv_from_z(igm_z, galaxy.z)

    # Rules
    match = np.where((rho<R_max) & (np.abs(dv) < dv_max))[0]
    if len(match) == 0:
        print("No CGM objects match your rules")
        return []
    else:
        # Loop to generate
        cgm_list = []
        for imatch in match:
            cgm = CGMAbsSys(galaxy, igmsystems[imatch], cosmo=cosmo, **kwargs)
            cgm_list.append(cgm)

    # Return
    return cgm_list


