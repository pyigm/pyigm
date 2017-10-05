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
from pyigm.abssys.igmsys import IGMSystem


def calc_cgm_rho(galaxy, igm_sys, cosmo, ang_sep=None, correct_lowz=True):
    """ Calculate the impact parameter between the galaxy and IGM sightline
    Mainly a wrapper to pyigm.utils.calc_rho

    Parameters
    ----------
    galaxy : Galaxy
    igm_sys : IGMSystem or list
    cosmo : astropy.cosmology
    ang_sep : Quantity, optional
    correct_lowz : bool, optional

    Returns
    -------
    rho : Quantity
      impact parameter(s) in physical kpc
    ang_sep : Angle
      separation in arcsec

    """
    from pyigm.utils import calc_rho
    # Loop?
    if isinstance(igm_sys, list):
        coords = SkyCoord([iigm.coord for iigm in igm_sys])
        return calc_rho(galaxy.coord, coords, np.array([galaxy.z]*len(coords)), cosmo, correct_lowz=correct_lowz)
    elif isinstance(igm_sys, IGMSystem):
        return calc_rho(igm_sys.coord, galaxy.coord, galaxy.z, cosmo,
                        ang_sep=ang_sep, correct_lowz=correct_lowz)
    else:
        raise IOError("Bad input..  Must be list or IGMSystem")

def get_close_galaxies(field,rholim=300.*u.kpc,minz=0.001,maxz=None):
    '''
    Generate table of galaxies close to sightlines for an IgmGalaxyField

    Parameters
    ----------
    fields : IgmGalaxyField or list of IgmGalaxyField
      Fields to be included
    rholim : quantity, optional
      Maximum impact parameter limit of galaxies
    minz : float, optional
      Minimum redshift of galaxies returned in query
    maxz : float, optional
      Maximum redshift of galaxies returned in query

    Returns
    -------
    closegaltab : astropy Table
      Table of galaxies meeting impact parameter and redshift criteria
    '''
    if maxz is None:
        maxz = np.inf
    field.galaxies['rho'] = field.calc_rhoimpact(field.galaxies,comoving=False)
    closecrit = field.galaxies['rho'].quantity<rholim
    zcrit = (field.galaxies['Z']>minz)&(field.galaxies['Z']<maxz)
    closegaltab=field.galaxies[closecrit&zcrit]
    if len(closegaltab) == 0:
        print('No galaxies found meeting these selection criteria.')
    return closegaltab


def cgmabssys_from_sightline_field(field,sightline,rholim=300.*u.kpc,
                                    minz=0.001,maxz=None,dv_max=400.*u.km/u.s):
    """Instantiate list of CgmAbsSys objects from IgmgGalaxyField and IGMSightline.

    Parameters
    ----------
    field : IgmGalaxyField
    sightline : IGMSightline
    rholim : Quantity, optional
        Maximum impact parameter for associated galaxies
    minz : float, optional
        Minimum redshift for galaxy/absorber search
    maxz : float, optional
        Maximum redshift for galaxy/absorber search
    dv_max : Quantity, optional
        Maximum galaxy-absorber velocity separation

    Returns
    -------
    cgmabslist : list
        List of CgmAbsSys objects
    """

    closegals = get_close_galaxies(field,rholim,minz,maxz)
    cgmabslist = []
    for i,gal in enumerate(closegals):
        galobj = Galaxy((gal['RA'],gal['DEC']),z=gal['Z'])
        cgmobj = cgm_from_galaxy_igmsystems(galobj,sightline._abssystems,
                                                dv_max=dv_max)
        cgmabslist.extend(cgmobj)
    return cgmabslist


def cgmsurvey_from_sightlines_fields(fields,sightlines,name=None,**kwargs):
    """Instantiate CGMAbsSurvey object from lists fo IgmGalaxyFields and IGMSightlines

    Parameters
    ----------
    fields: list of IgmGalaxyFields
    sightlines : list of IGMSightlines
    name : str, optional
        Name for the survey

    Returns
    -------
    cgmsurvey: CGMAbsSurvey
    """

    if ((not isinstance(fields,list))|(not isinstance(sightlines,list))|
        (len(fields) != len(sightlines))):
        raise IOError("Inputs fields and sightlines must lists of the same length")

    from pyigm.cgm.cgmsurvey import CGMAbsSurvey
    cgmsys = []
    for i,ff in enumerate(fields):
        thiscgmlist = cgmabssys_from_sightline_field(ff,sightlines[i],**kwargs)
        cgmsys.extend(thiscgmlist)
    if name is not None:
        cgmsurvey=CGMAbsSurvey.from_cgmabssys(cgmsys,survey=name)
    else:
        cgmsurvey = CGMAbsSurvey.from_cgmabssys(cgmsys)
    return cgmsurvey


def get_close_galaxies(field,rholim=300.*u.kpc,minz=0.001,maxz=None):
    '''
    Generate table of galaxies close to sightlines for an IgmGalaxyField

    Parameters
    ----------
    fields : IgmGalaxyField or list of IgmGalaxyField
      Fields to be included
    rholim : quantity, optional
      Maximum impact parameter limit of galaxies
    minz : float, optional
      Minimum redshift of galaxies returned in query
    maxz : float, optional
      Maximum redshift of galaxies returned in query

    Returns
    -------
    closegaltab : astropy Table
      Table of galaxies meeting impact parameter and redshift criteria
    '''
    if maxz is None:
        maxz = np.inf
    field.galaxies['rho'] = field.calc_rhoimpact(field.galaxies,comoving=False)
    closecrit = field.galaxies['rho'].quantity<rholim
    zcrit = (field.galaxies['Z']>minz)&(field.galaxies['Z']<maxz)
    closegaltab=field.galaxies[closecrit&zcrit]
    if len(closegaltab) == 0:
        print('No galaxies found meeting these selection criteria.')
    return closegaltab


def cgmabssys_from_sightline_field(field,sightline,rholim=300.*u.kpc,
                                    minz=0.001,maxz=None,dv_max=400.*u.km/u.s):
    """Instantiate list of CgmAbsSys objects from IgmgGalaxyField and IGMSightline.

    Parameters
    ----------
    field : IgmGalaxyField
    sightline : IGMSightline
    rholim : Quantity, optional
        Maximum impact parameter for associated galaxies
    minz : float, optional
        Minimum redshift for galaxy/absorber search
    maxz : float, optional
        Maximum redshift for galaxy/absorber search
    dv_max : Quantity, optional
        Maximum galaxy-absorber velocity separation

    Returns
    -------
    cgmabslist : list
        List of CgmAbsSys objects
    """

    closegals = get_close_galaxies(field,rholim,minz,maxz)
    cgmabslist = []
    for i,gal in enumerate(closegals):
        galobj = Galaxy((gal['RA'],gal['DEC']),z=gal['Z'])
        cgmobj = cgm_from_galaxy_igmsystems(galobj,sightline._abssystems,
                                                dv_max=dv_max)
        cgmabslist.extend(cgmobj)
    return cgmabslist


def cgmsurvey_from_sightlines_fields(fields,sightlines,name=None,**kwargs):
    """Instantiate CGMAbsSurvey object from lists fo IgmGalaxyFields and IGMSightlines

    Parameters
    ----------
    fields: list of IgmGalaxyFields
    sightlines : list of IGMSightlines
    name : str, optional
        Name for the survey

    Returns
    -------
    cgmsurvey: CGMAbsSurvey
    """

    if ((not isinstance(fields,list))|(not isinstance(sightlines,list))|
        (len(fields) != len(sightlines))):
        raise IOError("Inputs fields and sightlines must lists of the same length")

    from pyigm.cgm.cgmsurvey import CGMAbsSurvey
    cgmsys = []
    for i,ff in enumerate(fields):
        thiscgmlist = cgmabssys_from_sightline_field(ff,sightlines[i],**kwargs)
        cgmsys.extend(thiscgmlist)
    if name is not None:
        cgmsurvey=CGMAbsSurvey.from_cgmabssys(cgmsys,survey=name)
    else:
        cgmsurvey = CGMAbsSurvey.from_cgmabssys(cgmsys)
    return cgmsurvey


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
    rho, angles = calc_cgm_rho(galaxy, igmsystems, cosmo)

    # dv
    igm_z = np.array([igmsystem.zabs for igmsystem in igmsystems])
    dv = ltu.dv_from_z(igm_z, galaxy.z)

    # Rules
    match = np.where((rho<R_max) & (np.abs(dv) < dv_max))[0]
    if len(match) == 0:
        print("No IGMSystem paired to this galaxy. CGM object not created.")
        return []
    else:
        # Loop to generate
        cgm_list = []
        for imatch in match:
            cgm = CGMAbsSys(galaxy, igmsystems[imatch], cosmo=cosmo, **kwargs)
            cgm_list.append(cgm)

    # Return
    return cgm_list


