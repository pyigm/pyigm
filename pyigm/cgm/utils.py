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
from linetools.lists.linelist import LineList
from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine

from pyigm.field.galaxy import Galaxy
from pyigm.abssys.igmsys import IGMSystem

ismlist = LineList('ISM')

def calc_rho(galaxy, igm_sys, cosmo, ang_sep=None, correct_lowz=True,
             Galactic=False, d_Sun=8.0*u.kpc, **kwargs):
    """ Calculate the impact parameter between the galaxy and IGM sightline

    Parameters
    ----------
    galaxy : Galaxy object
    igm_sys : IGMSystem object or list
    cosmo : astropy.cosmology
    Galactic : bool, optional
      Calculate for our Galaxy!
    correct_lowz : bool, optional
      Apply corrections for the local universe, as desired
      Follows Mould et al. 2000
    ang_sep : Angle or Quantity
      Input angular separation
      May speed up calculation
    d_Sun : Quantity
      Distance to the Sun;  for Galactic calculation only


    Returns
    -------
    rho : Quantity
      impact parameter in kpc
    ang_sep : Angle
      separation in arcsec (deg for Galactic)

    """
    # Loop?
    if isinstance(igm_sys, list):
        rhos = []
        angs = []
        for iigm in igm_sys:
            irho, iang = calc_rho(galaxy, iigm, cosmo)
            rhos.append(irho.value)
            angs.append(iang)
        return np.array(rhos)*u.kpc, angs
    # Galactic?
    if Galactic:
        cosl_cosb = (np.cos(igm_sys.coord.transform_to('galactic').l)*
                          np.cos(igm_sys.coord.transform_to('galactic').b))
        xcomp = d_Sun * (cosl_cosb**2 - 1.)
        ycomp = d_Sun*(cosl_cosb*np.sin(igm_sys.coord.transform_to('galactic').l)*
                       np.cos(igm_sys.coord.transform_to('galactic').b))
        zcomp = d_Sun*(cosl_cosb*np.sin(igm_sys.coord.transform_to('galactic').b))
        # Distance
        d = np.sqrt(xcomp**2 + ycomp**2 + zcomp**2)niti
        # Angle
        ang_sep = np.arccos(cosl_cosb)
        return d, ang_sep.to('deg')

    if ang_sep is None:
        ang_sep = igm_sys.coord.separation(galaxy.coord).to('arcsec')
    # Handle cases where object's distance needs correction from peculiar velocities
    # This is especially important at very low redshifts
    if (galaxy.z < 0.05) and correct_lowz:
        velcorrdict = velcorr_mould(galaxy,cosmo=cosmo)
        kpc_amin = velcorrdict['scale'].to(u.kpc/u.arcmin)
        rho = ang_sep.to('arcmin') * kpc_amin
    else:
        kpc_amin = cosmo.kpc_comoving_per_arcmin(galaxy.z)  # kpc per arcmin
        rho = ang_sep.to('arcmin') * kpc_amin / (1+galaxy.z)
    # Return
    return rho, ang_sep


def get_close_galaxies(field,rho_max=300.*u.kpc,minz=0.001,maxz=None):
    '''
    Generate table of galaxies close to sightlines for an IgmGalaxyField

    Parameters
    ----------
    fields : IgmGalaxyField or list of IgmGalaxyField
      Fields to be included
    rho_max : quantity, optional
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
    closecrit = field.galaxies['rho'].quantity<rho_max
    zcrit = (field.galaxies['Z']>minz)&(field.galaxies['Z']<maxz)
    closegaltab=field.galaxies[closecrit&zcrit]
    if len(closegaltab) == 0:
        print('No galaxies found meeting these selection criteria.')
    return closegaltab


def cgmabssys_from_sightline_field(field,sightline,rho_max=300.*u.kpc,minz=0.001,
                                   maxz=None,dv_max=400.*u.km/u.s,dummysys=True,
                                   dummyspec=None,**kwargs):
    """Instantiate list of CgmAbsSys objects from IgmgGalaxyField and IGMSightline.

    Parameters
    ----------
    field : IgmGalaxyField
    sightline : IGMSightline
    rho_max : Quantity, optional
        Maximum impact parameter for associated galaxies
    minz : float, optional
        Minimum redshift for galaxy/absorber search
    maxz : float, optional
        Maximum redshift for galaxy/absorber search
    dv_max : Quantity, optional
        Maximum galaxy-absorber velocity separation
    dummysys : bool, optional
        Passed on to 'cgm_from_galaxy_igmsystems()'.  If True, create CGMAbsSyS
        even if no matching IGMSystem is found in any sightline for some galaxy.
    dummyspec : XSpectrum1D, optional
        Spectrum object to attach to dummy AbsLine/AbsComponent objects when
        adding IGMSystems if dummysys is True

    Returns
    -------
    cgmabslist : list
        List of CgmAbsSys objects
    """
    if dummyspec is None:
        dummyspec = sightline._abssystems[0]._components[0]._abslines[0].analy['spec']

    closegals = get_close_galaxies(field,rho_max,minz,maxz)
    cgmabslist = []
    for i,gal in enumerate(closegals):
        galobj = Galaxy((gal['RA'],gal['DEC']),z=gal['Z'])
        cgmobj = cgm_from_galaxy_igmsystems(galobj,sightline._abssystems,
                                            dv_max=dv_max, dummysys=dummysys,
                                            dummyspec=dummyspec, rho_max=rho_max,
                                            **kwargs)
        cgmabslist.extend(cgmobj)
    return cgmabslist


def cgmsurvey_from_sightlines_fields(fields, sightlines, rho_max=300*u.kpc,
                                     name=None, dummysys=True,  **kwargs):
    """Instantiate CGMAbsSurvey object from lists fo IgmGalaxyFields and IGMSightlines

    Parameters
    ----------
    fields: list of IgmGalaxyFields
    sightlines : list of IGMSightlines
    name : str, optional
        Name for the survey
    dummysys : bool, optional
        Passed on to 'cgm_from_galaxy_igmsystems()'.  If True, create CGMAbsSyS
        even if no matching IGMSystem is found in any sightline for some galaxy

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
        thiscgmlist = cgmabssys_from_sightline_field(ff,sightlines[i],rho_max=rho_max,
                                                     dummysys=dummysys, **kwargs)
        cgmsys.extend(thiscgmlist)
    if name is not None:
        cgmsurvey=CGMAbsSurvey.from_cgmabssys(cgmsys,survey=name)
    else:
        cgmsurvey = CGMAbsSurvey.from_cgmabssys(cgmsys)
    return cgmsurvey


def cgm_from_galaxy_igmsystems(galaxy, igmsystems, rho_max=300*u.kpc, dv_max=400*u.km/u.s,
                               cosmo=None, dummysys=False, dummyspec=None, **kwargs):
    """ Generate a list of CGMAbsSys objects given an input galaxy and a list of IGMSystems

    Parameters
    ----------
    galaxy : Galaxy
    igmsystems : list
      list of IGMSystems
    rho_max : Quantity
      Maximum projected separation from sightline to galaxy
    dv_max
      Maximum velocity offset between system and galaxy
    dummysys: bool, optional
        If True, instantiate CGMAbsSys even if no match is found in igmsystems
    dummyspec : XSpectrum1D, optional
        Spectrum object to attach to dummy AbsLine/AbsComponent objects when
        adding IGMSystems if dummysys is True.

    Returns
    -------
    cgm_list : list
      list of CGM objects generated

    """
    from pyigm.cgm.cgm import CGMAbsSys
    # Cosmology
    if cosmo is None:
        cosmo = cosmology.Planck15

    if dummysys is True:
        if dummyspec is None:
            dummyspec = igmsystems[0]._components[0]._abslines[0].analy['spec']
        dummycoords = igmsystems[0].coord
        llist = ismlist

    # R
    rho, angles = calc_rho(galaxy, igmsystems, cosmo)

    # dv
    igm_z = np.array([igmsystem.zabs for igmsystem in igmsystems])
    dv = ltu.dv_from_z(igm_z, galaxy.z)

    # Rules
    match = np.where((rho<rho_max) & (np.abs(dv) < dv_max))[0]
    if len(match) == 0:
        if dummysys is False:
            print("No IGMSystem paired to this galaxy. CGM object not created.")
            return []
        else:
            print("No IGMSystem match found. Attaching dummy IGMSystem.")
            dummysystem = IGMSystem(dummycoords,galaxy.z,vlim=None)
            dummycomp = AbsComponent(dummycoords,(1,1),galaxy.z,[-100.,100.]*u.km/u.s)
            dummyline = AbsLine('HI 1215',linelist = llist)  # Need an actual transition for comp check
            dummyline.analy['spec'] = dummyspec
            dummyline.attrib['coord'] = dummycoords
            dummycomp.add_absline(dummyline,chk_vel=False,chk_sep=False)
            dummysystem.add_component(dummycomp,chk_vel=False)
            cgm = CGMAbsSys(galaxy, dummysystem, cosmo=cosmo, **kwargs)
            cgm_list = [cgm]
    else:
        # Loop to generate
        cgm_list = []
        for imatch in match:
            cgm = CGMAbsSys(galaxy, igmsystems[imatch], cosmo=cosmo, **kwargs)
            cgm_list.append(cgm)

    # Return
    return cgm_list

def velcorr_mould(galaxy,cosmo=None):
    '''
    Calculates angular diameter distance, luminosity distance, etc., corrected
    for peculiar motions due to the Shapley Supercluster, the Local Group, etc.,
    using the formalism of Mould et al. (2000..ApJ..529..786).

    Code ported from IDL implementation by John Moustakas:
    github.com/moustakas/impro/blob/master/pro/galaxy/mould_distance.pro

    Parameters
    ----------
    galaxy: Galaxy
    cosmo: astropy.cosmology, optional

    Returns
    -------
    velcorrdict: dictionary
        Resulting quantities of velocity corrected calculation.
        v_hel: Heliocentric velocity of object with no correction, simply z*c
        v_LG: Velocity with respect to Local Group
        v_virgo: Velocity with respect to Virgo Cluster
        v_GA: Velocity with respect to Virgo Cluster
        v_shapley: Velocity with respect to Shapley Supercluster
        v_cosmic: Corrected velocity
        lumdist: Luminosity distance
        angdiamdist: Angular diameter distance
        distmod: Distance modulus
        flag: 1 if object is the direction of Virgo, 2 if GA, 3 if Shapley
        Scale: Proper distance per angular separation on sky
    '''

    # Cosmology
    if cosmo is None:
        cosmo = cosmology.Planck15

    H0 = cosmo.H0
    omega0 = cosmo.Om0
    omega_lambda = cosmo.Ode0

    # Needed parameters
    q0 = omega0 / 2.0 - omega_lambda
    gamma = 2.0

    # Load info from Mould+ 2000 Table 1A
    clnames = ['Virgo', 'GA', 'Shapley']
    clra1950 =  ['12h28m19s','13h20m00s','13h30m00s']
    cldec1950 = ['+12:40:00','-44:00:00','-31:00:00']
    clcoords = SkyCoord(clra1950,cldec1950,unit=(u.hourangle,u.degree),frame='fk4',equinox='B1950.0')
    clhelvel = np.array([1035., 4600., 13800.]) * u.km/u.s
    clLGvel = np.array([957., 4380., 13600.])  * u.km/u.s
    fidvel = np.array([200., 400., 85.]) * u.km/u.s
    clrad = np.array([10., 10., 12.]) * u.degree
    clrangelo = np.array([600., 2600., 10000.]) * u.km/u.s
    clrangehi = np.array([2300., 6600., 16000.]) * u.km/u.s

    # Convert B1950 coordinates (as given) and
    clcoords = clcoords.icrs

   # Convert input coords to Galactic coords
    galcoords_gal = galaxy.coord.transform_to(frame='galactic')

    # Transform to local group frame
    c = constants.c.to(u.km/u.s)
    v_LG = c * galaxy.z - 79.0 * (u.km/u.s) * np.cos(galcoords_gal.l).value \
        * np.cos(galcoords_gal.b).value + 296.0 * (u.km/u.s) * np.sin(galcoords_gal.l).value \
        * np.cos(galcoords_gal.b).value - 36.0 * (u.km/u.s) * np.sin(galcoords_gal.b).value

    # Calculate object-attractor angular and velocity distances (eq. 2 in Mould 2000+)
    theta = galaxy.coord.separation(clcoords)
    costheta = np.cos(theta).value
    r0a = np.sqrt(v_LG ** 2 + clLGvel ** 2 - 2. * v_LG * clLGvel * costheta)

    # Determine if object is in direction of one of the attractors.  If not, calculate velocity!
    virgo = False
    GA = False
    shapley = False

    if (theta[0] < clrad[0]) & (
        ((clLGvel[0] - r0a[0]) > clrangelo[0]) | ((clLGvel[0] + r0a[0]) < clrangehi[0])):
        virgo = True
    if (theta[1] < clrad[1]) & (
        (clLGvel[1] - r0a[1] > clrangelo[1]) | (clLGvel[1] + r0a[1] < clrangehi[1])):
        GA = True
    if (theta[2] < clrad[2]) & (
        (clLGvel[2] - r0a[2] > clrangelo[2]) | (clLGvel[2] + r0a[2] < clrangehi[2])):
        shapley = True

    if virgo or GA or shapley:
        v_infall = np.zeros(3) * u.km/u.s
        if virgo:
            v_cosmic = clLGvel[0]
            flag = 1
        if GA:
            v_cosmic = clLGvel[1]
            flag = 2
        if shapley:
            v_cosmic = clLGvel[2]
            flag = 3
    else:
        v_infall = fidvel * (costheta + (v_LG - clLGvel * costheta) / r0a * (r0a / clLGvel) ** (1. - gamma))
        v_cosmic = v_LG + np.sum(v_infall)
        flag = 0

    # Derive remaining parameters to report
    z_cosmic = v_cosmic / c
    propdist = c / H0 * (z_cosmic - z_cosmic ** 2 / 2. * (1. + q0))
    lumdist = propdist * (1. + z_cosmic)
    angdiamdist = lumdist / (1. + z_cosmic) ** 2
    distmod = 5. * np.log10(lumdist.to(u.pc).value) - 5.
    kpcarcsec = angdiamdist * np.pi * 2. / 360. / 3600. * 1000.
    scale = angdiamdist * np.pi * 2./ (360. * u.degree)

    # Create dictionary to return
    velcorrdict = dict(v_hel=c * galaxy.z, v_LG=v_LG, v_virgo=v_infall[0], v_GA=v_infall[1],
                  v_shapley=v_infall[2], v_cosmic=v_cosmic, lumdist=lumdist, angdiamdist=angdiamdist,
                  distmod=distmod, flag=flag, scale=scale.to(u.kpc/u.arcsec))

    return velcorrdict


