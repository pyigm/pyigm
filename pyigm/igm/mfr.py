""" Module for mean-flux regulation
Port of K-G Lee IDL code
"""
import pdb

import numpy as np

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import FittableModel, Parameter

from pyigm.fN import tau_eff as pyteff

try:
    basestring
except NameError:  # For Python 3
    basestring = str


def fit_forest(wave, flux, sigma, zqso, wavemin=3200., normfac=None,
               user_mask=None, mask_dlas=None, coord=None,
               for_lohi=(1041.,1185)):
    """  Perform mean-flux regulation on an input spectrum

    Parameters
    ----------
    wave : ndarray
      Observed wavelength values in Ang
    flux : ndarray
      flux array
    sigma : ndarray
      sigma array.  0 values are masked
    zqso : float
      Quasar emission redshift
    wavemin : float, optional
      Minimum wavelength to include in analysis
    normfac : float, optional
      optional scaling factor
    user_mask : list (of lists)
      User can mask out regions in the anaylsis `by-hand'
    mask_dlas : Table or str, optional
      Table of DLAs to mask based on coordinates and redshift
        Required columns are RA, DEC, z
      Can be the filename of the table to input
    coord : SkyCoord, optional
      Coordinate of the sightline.  Used to match against DLAs
    for_lohi : tuple
      wavemin, wavemax that define the forest for analysis

    Returns
    -------
    new_ff : ndarray
      Mean-flux regulated flux
    parm : FittableModel
      astropy model of the continuum
    """
    from astropy.modeling import fitting

    # invvar
    ivar = (sigma != 0.) / (sigma**2 + (sigma == 0))

    # Scale?
    if normfac is not None:
        ivar = ivar * normfac**2
        flux = flux / normfac

    # User Mask?
    if user_mask is not None:
        for imask in user_mask:
            msk = (wave >= imask[0]) & (wave <= imask[1])
            ivar[msk] = 0.

    # DLA mask?
    if mask_dlas is not None:
        if coord is None:
            raise IOError("Must input coord to mask DLAs")
        # Load DLAs, as needed
        if isinstance(mask_dlas, basestring):
            dlas = Table.read(mask_dlas)
        elif isinstance(mask_dlas, Table):
            dlas = mask_dlas
            # Check for coord column
            if 'coord' not in dlas.keys():
                dla_coord = SkyCoord(ra=dlas['RA'], dec=dlas['DEC'], unit='deg')
            else:
                dla_coord = dlas['coord']
        else:
            raise IOError("Not ready for this type of DLA input")

        # Search
        sep = coord.separation(dla_coord)
        match = np.where(sep < 2*u.arcsec)[0]
        for imatch in match:
            zdla = dlas['z'][imatch]
            skip_fg = np.abs(wave-1215.6701*(1+zdla)) < 100.
            ivar[skip_fg] = 0.
            #print, 'Removed a DLA'

    # Isolate forest for analysis
    lambda_r = wave / (1+zqso)
    forestrange = (lambda_r >= for_lohi[0]) & (lambda_r <= for_lohi[1]) & (
                    wave > wavemin)
    lamb_forest = lambda_r[forestrange]
    z_for = (lamb_forest/1215.67)*(1.+zqso) -1.  # Redshift at each pixel
    fforest = flux[forestrange]
    ivarforest = ivar[forestrange]

    # Estimate weights for each pixel
    var_F = forestvar(z_for) * (np.exp(-pyteff.lyman_alpha_obs(z_for)))**2
    #var_F = forestvar(z_for) * (np.exp(-old_taueff_evo(z_for)))**2
    var_noise = (ivarforest != 0) / (ivarforest + (ivarforest == 0))
    var_total = var_F + var_noise
    weights_forest = (var_total != 0) / (var_total + (var_total == 0))

    # But need to make sure that masked pixels remain masked...
    maskedpix = ivarforest == 0
    weights_forest[maskedpix] = 0

    # Astropy modeling
    model = mflux_tauevo(p0=0., p1=0., zqso=zqso, lamb_piv=1113.)
    fitter = fitting.LevMarLSQFitter()
    parm = fitter(model, lamb_forest, fforest, weights=weights_forest)

    # Apply
    full_forest = lambda_r < 1220.
    new_ff = flux
    new_ff[full_forest] = flux[full_forest] / mfluxcorr(lambda_r[full_forest],
                                                        parm.p0.value, parm.p1.value,
                                                        lamb_piv=parm.lamb_piv.value)
    # Return
    return new_ff, parm


def forestvar(z_in):
    """  Return intrinsic variance of LyaF variance for weighting. This
    estimate is roughly from McDonald et al 2006

    Parameters
    ----------
    z_in : float or ndarray

    Returns
    -------
    fvar : float or ndarray
      Variance
    """
    fvar = 0.065 * ((1.+z_in)/(1.+2.25))**3.8
    # Return
    return fvar


def mfluxcorr(lambda_r, p0, p1, lamb_piv=1113.):
    """ Correction factor to power-law fit
    Parameters
    ----------
    lambda_r
    p
    lamb_piv

    Returns
    -------

    """
    lamb_piv = 1113. # This is the pivot point in the restframe spectrum
    return p0 + p1*(lambda_r/lamb_piv - 1.)


class mflux_tauevo(FittableModel):
    """ Mean flux evolution * exp(delta*(lambda/1280-1)),
    Meant for use with astropy.modeling to correct the fitted Lya forest continuum
    Abscissa parameter x is the restframe wavelength, and free parameters p0, p1
    set the power law.
    This function NEEDS the quasar redshift zqso to be set

    input: lambda_r :: Rest wavelength Assumed in Angstroms
    output: absorbed, normalized flux
    Parameters: logN,b,z,wrest,f,gamma,fwhm
    """
    inputs = ('lambda_r',)
    outputs = ('flux',)

    # Free parameters (generally)
    p0 = Parameter()
    p1 = Parameter()

    # Fixed parameters
    zqso = Parameter(fixed=True)
    lamb_piv = Parameter(fixed=True)

    @staticmethod
    def evaluate(lambda_r, p0, p1, zqso, lamb_piv): #logN,b,z,wrest,f,gamma,fwhm):
        zfor = (lambda_r/1216.) * (1. + zqso) - 1.
        tau = pyteff.lyman_alpha_obs(zfor)
        #tau = old_taueff_evo(zfor)
        fmean = np.exp(-1*tau)

        mfluxtauevo = fmean * mfluxcorr(lambda_r, p0, p1, lamb_piv=lamb_piv)
        return mfluxtauevo


def old_taueff_evo(z):
    """ F-G taueff. Mainly for testing
    Parameters
    ----------
    z

    Returns
    -------

    """
    tauevo  = 0.001845 * (1+z)**3.924
    return tauevo