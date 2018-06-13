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
from linetools.analysis.zlimits import zLimits

from pyigm.field.galaxy import Galaxy
from pyigm.abssys.igmsys import IGMSystem


def calc_cgm_rho(galaxy, igm_sys, cosmo, **kwargs):
    """ Calculate the impact parameter between the galaxy and IGM sightline
    Mainly a wrapper to pyigm.utils.calc_rho

    Parameters
    ----------
    galaxy : Galaxy
    igm_sys : IGMSystem or list
    cosmo : astropy.cosmology
    **kwargs :
      Passed to pyigm.utils.calc_rho

    Returns
    -------
    rho : Quantity
      impact parameter(s) in comoving or physical kpc
    ang_sep : Angle
      separation in arcsec

    """
    from pyigm.utils import calc_rho
    # Loop?
    if isinstance(igm_sys, list):
        coords = SkyCoord([iigm.coord for iigm in igm_sys])
        return calc_rho(galaxy.coord, coords, np.array([galaxy.z]*len(coords)), cosmo, **kwargs)
    elif isinstance(igm_sys, IGMSystem):
        return calc_rho(igm_sys.coord, galaxy.coord, galaxy.z, cosmo, **kwargs)
    else:
        raise IOError("Bad input.  Must be list or IGMSystem")


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
                                   maxz=None,dv_max=400.*u.km/u.s,embuffer=None,
                                   dummysys=True,dummyspec=None,linelist=None,
                                   debug=False,
                                   **kwargs):
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
    embuffer : Quantity, optional
        Velocity buffer between background source (e.g., QSO) and CGMAbsSys
    dummysys : bool, optional
        Passed on to 'cgm_from_galaxy_igmsystems()'.  If True, create CGMAbsSyS
        even if no matching IGMSystem is found in any sightline for some galaxy.
    dummyspec : XSpectrum1D, optional
        Spectrum object to attach to dummy AbsLine/AbsComponent objects when
        adding IGMSystems if dummysys is True
    linelist : LineList, optional
        ListList from which to add dummy line in case of no IGMSystem match

    Returns
    -------
    cgmabslist : list
        List of CgmAbsSys objects
    """
    if dummyspec is None:
        dummyspec = sightline._abssystems[0]._components[0]._abslines[0].analy['spec']

    if linelist is None:
        from linetools.spectralline import AbsLine
        linelist = LineList('ISM')

    if embuffer is not None:
        try:
            bufmax = ltu.z_from_dv(-embuffer,field.zem)
            if maxz is not None:
                zmax = np.max(bufmax,maxz)
            else:
                zmax = bufmax
        except:
            zmax = maxz
    else:
        zmax = maxz


    closegals = get_close_galaxies(field,rho_max,minz,zmax)
    if debug:
        closegals = closegals[0:10]
    cgmabslist = []
    for i,gal in enumerate(closegals):
        #print('i={:d}'.format(i))
        galobj = Galaxy((gal['RA'],gal['DEC']),z=gal['Z'])
        cgmobj = cgm_from_galaxy_igmsystems(galobj,sightline._abssystems,
                                            dv_max=dv_max, dummysys=dummysys,
                                            dummyspec=dummyspec, rho_max=rho_max,
                                            linelist=linelist,**kwargs)
        cgmabslist.extend(cgmobj)
    return cgmabslist


def cgmsurvey_from_sightlines_fields(fields, sightlines, rho_max=300*u.kpc,
                                     name=None, dummysys=True, embuffer=None,
                                     **kwargs):
    """Instantiate CGMAbsSurvey object from lists fo IgmGalaxyFields and IGMSightlines

    Parameters
    ----------
    fields: list of IgmGalaxyFields
    sightlines : list of IGMSightlines
    name : str, optional
        Name for the survey
    dummysys : bool, optional
        Passed on to 'cgm_from_galaxy_igmsystems()'.  If True, create CGMAbsSys
        even if no matching IGMSystem is found in any sightline for some galaxy
    embuffer : Quantity, optional
        Velocity buffer between background source (e.g., QSO) and CGMAbsSys

     Returns
    -------
    cgmsurvey: CGMAbsSurvey
    """
    if ((not isinstance(fields,list))|(not isinstance(sightlines,list))|
        (len(fields) != len(sightlines))):
        raise IOError("Inputs fields and sightlines must lists of the same length")

    if dummysys is True:
        from linetools.spectralline import AbsLine
        ismlist = LineList('ISM')

    from pyigm.cgm.cgmsurvey import CGMAbsSurvey
    cgmsys = []
    for i,ff in enumerate(fields):
        print(ff.name)
        thiscgmlist = cgmabssys_from_sightline_field(ff,sightlines[i],rho_max=rho_max,
                                                     dummysys=dummysys,linelist=ismlist,
                                                     embuffer=embuffer,**kwargs)
        cgmsys.extend(thiscgmlist)
    if name is not None:
        cgmsurvey=CGMAbsSurvey.from_cgmabssys(cgmsys,survey=name)
    else:
        cgmsurvey = CGMAbsSurvey.from_cgmabssys(cgmsys)
    return cgmsurvey


def cgm_from_galaxy_igmsystems(galaxy, igmsystems, rho_max=300*u.kpc, dv_max=400*u.km/u.s,
                               cosmo=None, dummysys=False, dummyspec=None, verbose=True, **kwargs):
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

    # R -- speed things up
    rho, angles = calc_cgm_rho(galaxy, igmsystems, cosmo, **kwargs)
    if len(igmsystems) == 1:  # Kludge
        rho = u.Quantity([rho])
        angles = u.Quantity([angles])


    # dv
    igm_z = np.array([igmsystem.zabs for igmsystem in igmsystems])
    dv = ltu.dv_from_z(igm_z, galaxy.z)

    # Rules
    match = np.where((rho<rho_max) & (np.abs(dv) < dv_max))[0]

    ### If none, see if some system has a component that's actually within dv_max
    if (len(match) == 0)&(rho[0] < rho_max):
        zcomps = []
        sysidxs = []
        for i, csys in enumerate(igmsystems):
            thesezs = [comp.zcomp for comp in csys._components]
            sysidxs.extend([i] * len(thesezs))
            zcomps.extend(thesezs)
        zcomps = np.array(zcomps)
        sysidxs = np.array(sysidxs)
        dv_comps = ltu.dv_from_z(zcomps,galaxy.z)
        match = np.unique(sysidxs[np.where(np.abs(dv_comps)<dv_max)[0]])

    if len(match) == 0:
        if dummysys is False:
            print("No IGMSystem paired to this galaxy. CGM object not created.")
            return []
        else:
            if verbose:
                print("No IGMSystem match found. Attaching dummy IGMSystem.")
            dummysystem = IGMSystem(dummycoords,galaxy.z,vlim=None)
            dummycomp = AbsComponent(dummycoords,(1,1),galaxy.z,[-100.,100.]*u.km/u.s)
            dummycomp.flag_N = 3
            dummyline = AbsLine('HI 1215',**kwargs)  # Need an actual transition for comp check
            dummyline.analy['spec'] = dummyspec
            dummyline.attrib['coord'] = dummycoords
            dummycomp.add_absline(dummyline,chk_vel=False,chk_sep=False)
            dummysystem.add_component(dummycomp,chk_vel=False,chk_sep=False)
            cgm = CGMAbsSys(galaxy, dummysystem, cosmo=cosmo, **kwargs)
            cgm_list = [cgm]
    else:
        # Loop to generate
        cgm_list = []
        for imatch in match:
            # Instantiate new IGMSystem
            # Otherwise, updates to the IGMSystem cross-pollinate other CGMs
            sysmatch = igmsystems[imatch]
            newisys = sysmatch.copy()
            # Handle z limits
            zlim = ltu.z_from_dv((-dv_max.value,dv_max.value)*u.km/u.s,galaxy.z)
            newlims = zLimits(galaxy.z,zlim.tolist())
            newisys.limits = newlims
            # Allow for components extending beyond dv_max
            newisys.update_vlim()
            newisys.update_component_vel()
            # Finish
            cgm = CGMAbsSys(galaxy, newisys, cosmo=cosmo, rho=rho[imatch], ang_sep=angles[imatch], **kwargs)
            cgm_list.append(cgm)

    # Return
    return cgm_list

def covering_fraction(iontable,colthresh,sys_attr='rho_impact',bins=None,
                      returncounts=False,**kwargs):
    """Given some detection threshold, calculate the ion covering fraction
    within bins of some quantity (impact parameter by default)

    Parameters
    ----------
    iontable : Table
        Output of CGMAbsSurvey.ion_tbl()
    colthresh : float or None
        Detection threshold in log column density; if None, impose no threshold
        Nondetections with upper limits above this threshold are ignored
    sys_attr : str, optional
        Column name in iontable corresponding to given bins 
    bins,optional : list of ints or floats
        Bins in impact parameter (units of kpc) or other quantity
    returncounts : bool, optional
        If True, return numbers of hits and total systems within each bin

    Returns
    -------
    fracs : list of floats
        Covering fraction for each bin
    lolims : list of floats
        Lower limits on covering fractions
    uplims : list of floats
        Upper limits on covering fractions
    dethist : list of ints, optional
        Number of detections in each bin
    tothist : list of ints, optional
        Total number of systems in each bin

    """
    from pyigm import utils as pu

    # Clean up table to remove nulls if they exist
    itab = iontable.copy()
    goodones = np.where(itab['flag_N']>0)[0]
    itab = itab[goodones]

    flag = itab['flag_N']
    col = itab['logN']
    rho = itab[sys_attr]

    if bins is None:
        bins = [0,np.max(rho)]

    # Impose threshold if it exists
    if colthresh != None:
        dets = np.where(((flag == 1) | (flag == 2))&(col > colthresh))[0]
        nondets = np.where(((flag == 3) & (col <= colthresh)) |
                            ((flag == 1) & (col <= colthresh)) )[0]
    else:
        dets = np.where((flag == 1) | (flag == 2))[0]
        nondets = np.where((flag == 3))[0]

    # Set up and sort bins
    dethist, histbins = np.histogram(rho[dets], bins=bins)
    nondethist, histbins = np.histogram(rho[nondets], bins=bins)

    # Get the stats
    tothist = nondethist + dethist
    fracs, lolims, uplims = pu.confintervals(dethist, tothist,**kwargs)
    upbars = uplims - fracs
    lobars = fracs - lolims

    # Return
    if returncounts == False:
        return fracs, lobars, upbars
    else:
        return fracs, lobars, upbars, dethist, tothist