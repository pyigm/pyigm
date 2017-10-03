""" Utilities for IGM calculations
"""
import numpy as np
import pdb

from astropy import constants as const
from astropy import units as u
from astropy.units.quantity import Quantity
from astropy import cosmology
from astropy.coordinates import SkyCoord

from pyigm.field.galaxy import Galaxy

def calc_rho(coords1, coords2, z1, cosmo=None, ang_sep=None, correct_lowz=True,
             z_low=0.05):
    """ Calculate the impact parameter between the galaxy and IGM sightline

    Parameters
    ----------
    coords1 : SkyCoord (one or more)
    coords2 : SkyCoord (one or more)
    z1 : float or ndarray
      Redshifts of sources at coords1
    cosmo : astropy.cosmology, optional
    z_low : float, optional
      Redshift below which corrections for the local universe would be applied
    correct_lowz : bool, optional
      Apply corrections for the local universe, as desired
      Only those with z < z_low
      Follows Mould et al. 2000
    ang_sep : Angle or Quantity
      Input angular separation
      May speed up calculation

    Returns
    -------
    rho : Quantity
      impact parameter(s) in physical kpc
    ang_sep : Angle
      separation in arcsec
    """
    # Handle inputs
    if isinstance(z1,float):
        z1 = np.array([z1])
    if not (coords1.size == z1.size):
        raise IOError("Length of z1 must match coords1")
    if ang_sep is None:
        ang_sep = coords1.separation(coords2).to('arcsec')
    # Init rho
    rho = np.zeros_like(z1) * u.kpc
    # Handle cases where object's distance needs correction from peculiar velocities
    # This is especially important at very low redshifts
    lowz = z1 < 0.05
    if correct_lowz and np.any(lowz):
        # Ugly for loop for now
        for idx in np.where(lowz)[0]:
            # Deal with scalar vs. array
            if coords1.size == 1:
                icoord = coords1
            else:
                icoord = coords1[idx]
            # Call Mould correction
            pdb.set_trace()
            velcorrdict = velcorr_mould(Galaxy(icoord, z=z1[idx]), cosmo=cosmo)
            kpc_amin = velcorrdict['scale'].to(u.kpc/u.arcmin)
            rho[idx] = ang_sep[idx].to('arcmin') * kpc_amin
    else:
        kpc_amin = cosmo.kpc_comoving_per_arcmin(z1)  # kpc per arcmin
        rho = ang_sep.to('arcmin') * kpc_amin / (1+z1)
    # Return
    return rho, ang_sep


def calc_Galactic_rho(coords, d_Sun=8.0*u.kpc):
    """  Calculate impact parameter to Galactic center

    Parameters
    ----------
    coords
    d_Sun

    Returns
    -------
    d : Quantity
      impact parameter(s) from Galctic center in physical kpc
    ang_sep : Angle
      separation in deg
    """
    # Transform to Galactic
    gcoords = coords.to('galactic')
    # Calculate
    cosl_cosb = (np.cos(gcoords.l)*
                 np.cos(gcoords.b))
    xcomp = d_Sun * (cosl_cosb**2 - 1.)
    ycomp = d_Sun * (cosl_cosb*np.sin(gcoords.l)*
                   np.cos(coords.b))
    zcomp = d_Sun * (cosl_cosb*np.sin(gcoords.b))
    # Distance
    d = np.sqrt(xcomp**2 + ycomp**2 + zcomp**2)
    # Angle
    ang_sep = np.arccos(cosl_cosb)
    return d, ang_sep.to('deg')

def cosm_xz(z, cosmo=None, flg_return=0):
    """ Calculates X(z) -- absorption path length or dXdz

    Parameters
    ----------
    z : float or ndarray
       Redshift to evaluate at.  May be an array
    cosmo : astropy.cosmology
       Cosmological model to adopt
    flg_return : int, optional (0)
       Flag controlling the output
         0 = X(z)
         1 = dX/dz at z

    Returns
    -------
    Xz or dXdz :
      Absorption path or its derivative
      Set by flag
    """

    # Cosmology
    if cosmo is None:
        from astropy.cosmology import core as acc
        cosmo = acc.FlatLambdaCDM(70., 0.3)

    # Flat?
    if cosmo.Ok(0.) == 0:
        if flg_return == 0:  # X(z)
            rslt = cosmo.absorption_distance(z)
        elif flg_return == 1:  # dX/dz
            try:
                rslt = cosmo.abs_distance_integrand(z)
            except AttributeError:
                rslt = cosmo._xfunc(z)
        else:
            raise ValueError('pyigm.utils.cosm_xz: Bad flg {:d}'.format(flg_return))
    else:
        raise ValueError('pyigm.utils.cosm_xz: Not prepared for non-flat cosmology')

    #
    return rslt

def mk_ew_lyman_spline(bval, ew_fil=None, chk=False):
    """ Generate a pickle file of a Spline of EW vs NHI for the Lyman series

    Parameters
    ----------
    bval : float or Quantity
      Doppler parameter (km/s)
    ew_fil : string ('EW_SPLINE_b##.p')
      Name of output pickle file
    chk : bool, optional
        Check the calculation

    """
    import pickle
    import yaml
    from linetools.lists.linelist import LineList
    from linetools.analysis import voigt as lav
    from scipy import interpolate

    # Outfil
    if ew_fil == None:
        ew_fil = 'EW_SPLINE_b'+str(int(bval))+'.p'

    # Units
    if not isinstance(bval,u.quantity.Quantity):
        bval = bval * u.km/u.s  # km/s

    # NHI
    nspl = 100
    log_NHI = 11.0 + 11*np.arange(nspl)/(nspl-1.)

    # Lines
    HI = LineList('HI')
    wrest = HI._data['wrest']

    # Output
    owrest = [row['wrest'].value for row in HI._data]
    outp = {'wrest': owrest, 'tck': []}

    # Setup
    nvel = 60001
    velo = (-30000. + np.arange(nvel,dtype='float64'))*u.km/u.s # km/s
    dvel = 1. * u.km/u.s # km/s
    uval = velo / bval

    # Loop
    for cnt, line in enumerate(HI._data):

        # Wave array
        dwv = dvel.to(u.cm/u.s) * line['wrest'] / const.c.cgs  # Ang

        # Voigt
        vd = (bval/line['wrest']).to(u.Hz)  # Frequency
        a = line['gamma'].value / (12.56637 * vd.value)
        vgt = lav.voigt_wofz(uval.value, a)

        # tau
        tau = 0.014971475*line['f']*vgt/vd  # Normalized to N_HI = 1 cm^-2

        # Flux
        tau_array = np.outer(tau, 10.**log_NHI)
        fx = np.exp(-1.*tau_array)

        # EW
        EW = np.sum(1.-fx, 0) * dwv

        # Spline
        tck = interpolate.splrep(log_NHI, EW)

        # Check?
        if chk:
            from matplotlib import pyplot as plt
            plt.clf()
            plt.plot(log_NHI, EW, 'o')
            # Spline
            xnew = np.linspace(np.amin(log_NHI),np.amax(log_NHI), nspl*10)
            ynew = interpolate.splev(xnew, tck, der=0)
            plt.plot(xnew, ynew, '-')
            plt.show()

        # Output
        print('line = {:g}'.format(line['wrest']))
        outp['tck'].append(tck)

    # Write
    print('Writing {:s}'.format(ew_fil))
    with open(ew_fil, 'w') as yamlf:
        yamlf.write(yaml.dump(outp))

def lst_to_array(lst, mask=None):
    """ Simple method to convert a list to an array

    Allows for a list of Quantity objects

    Parameters
    ----------
    lst : list
      Should be number or Quantities
    mask : boolean array, optional

    Returns
    -------
    array or Quantity array

    """
    if mask is None:
        mask = np.array([True]*len(lst))
    if isinstance(lst[0], Quantity):
        return Quantity(lst)[mask]
    else:
        return np.array(lst)[mask]


def velcorr_mould(galaxy,cosmo=None):
    """
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
    """

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
    c = const.c.to(u.km/u.s)
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
