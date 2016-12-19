""" Module for tau_eff calculations
"""
import pdb
import yaml
import imp

import numpy as np
from scipy import interpolate

from astropy import constants as const
from astropy import units as u
from astropy import cosmology
from astropy.table import QTable

from linetools.analysis import absline as ltaa
from linetools.lists import linelist as lll
from linetools.analysis import voigt as lav
from linetools.analysis.absline import photo_cross
from linetools import spectralline
from linetools.spectra.xspectrum1d import XSpectrum1D

from pyigm.fN.fnmodel import FNModel
from pyigm import utils as pyigmu
from pyigm.continuum import quasar as pycq

pyigm_path = imp.find_module('pyigm')[1]


def dopp_val(x,bsig=24*u.km/u.s,bmnx=(15.,80)*u.km/u.s):
    """ Generate random distribution of b-values
    Follows Hui&Rutledge 1999 distribution

    Parameters
    ----------
    x: float or array  [between 0 and 1]
    bsig: Quantity, optional  [parameterization of the distribution]
    bmnx: Quantity array, optional
      Set limits on the allowed values

    Returns
    -------
    bx: Quantity
      b-values
    """
    bmin,bmax = bmnx
    # Normalization constant
    Const = (4*bsig**4) / (np.exp(-bsig**4/bmax**4)-np.exp(-bsig**4/bmin**4))
    # Now the evaluation
    bx4 = bsig**4 / (-1*np.log(4*x*bsig**4/Const + np.exp(-bsig**4/bmin**4)))
    #Lastly
    bx = bx4**(0.25)
    return bx


def monte_HIcomp( zmnx, fN_model, NHI_mnx=None, dz=0.001, cosmo=None,
    rstate=None, seed=None):
    """ Generate a Monte Carlo draw of HI components (z,N,b)

    Parameters
    ----------
    zmnx : tuple (float,float)
      Redshift range for analysis.
      Should correspond to Lya
    fN_model : fN_Model class
    NHI_mnx : tuple, optional (float,float)
      Range of logNHI for linelist
    dz : float, optional
      Step size for z discretization
    cosmo : astropy Cosmology, optional
    rstate : RandomState, optional
      For random number generation
    seed : int, optional
      Seed for random numbers

    Returns:
    -----------
    HI_comps : list
      List of HI components drawn for the given sightline
    """
    # Init
    # NHI range
    if NHI_mnx is None:
        NHI_mnx = (12., 22.)
    # seed
    if rstate is None:
        if seed is None:
            seed = 12345
        rstate = np.random.RandomState(seed)

    # Check fN_model type
    if fN_model.fN_mtype != 'Hspline':
        raise ValueError('monte_HIlines: Can only handle Hspline so far.')

    # Calculate lX at pivot
    lX, cum_lX, lX_NHI = fN_model.calculate_lox(fN_model.zpivot,
        NHI_mnx[0],NHI_max=NHI_mnx[1], cumul=True)

    # Interpolator for NHI distribution (assumed independent of redshift)
    #   Uses lowest NHI value for the first bit (kludgy but ok)
    invfN = interpolate.interp1d(cum_lX/lX,lX_NHI,bounds_error=False,fill_value=lX_NHI[0])#, kind='quadratic')

    # z evaluation
    zeval = np.arange(zmnx[0], zmnx[1]+dz, dz)

    # Cosmology
    if cosmo is None:
        print('Using a Flat LCDM cosmology: h=0.7, Om=0.3')
        cosmo = cosmology.core.FlatLambdaCDM(70., 0.3)

    # dXdz
    dXdz = pyigmu.cosm_xz(zeval, cosmo=cosmo, flg_return=1)

    # Generate loz
    loz = lX * dXdz * ( (1+zeval)/(1+fN_model.zpivot) )**fN_model.gamma

    # Calculate average number of lines for analysis
    sum_loz = np.cumsum(loz*dz)

    # Interpolator
    #   Uses lowest NHI value for the first bit (kludgy but ok)
    invz = interpolate.interp1d(sum_loz/sum_loz[-1],zeval, bounds_error=False, fill_value=zeval[0])

    # Assume Gaussian stats for number of lines
    nlines = int(np.round(sum_loz[-1] + np.sqrt(sum_loz[-1])*rstate.randn(1)))

    # NHI values
    randNHI = rstate.random_sample(nlines)
    lgNHI = invfN(randNHI)

    # z values
    randz = rstate.random_sample(nlines)
    zval = invz(randz)

    # b values
    randb = rstate.random_sample(nlines)
    bval = dopp_val(randb)

    # Pack em up as a QTable
    HI_comps = QTable([zval, lgNHI, bval], names=('z','lgNHI','bval'))
    return HI_comps


def generate_tau(iwave, HIlines, HI_comps, kludge=True):
    """Generate optical depth array given lines and components

    Parameters:
    -----------

    iwave : Quantity array
      Input Spectrum wavelengths
    HIlines : list of AbsLines
    HI_comps : QTable of components
    kludge : bool, optional
      Kludge the opacity

    Returns:
    --------
    tau : ndarray
      Optical depth at subgrid wavelengths.  Will need to rebin back
    """
    # Rebin to subgrid

    wmin = np.min(iwave.to('AA').value)
    wmax = np.max(iwave.to('AA').value)
    nsub = int(np.round( (np.log10(wmax)- np.log10(wmin)) / 1.449E-6)) + 1
    wave = 10.**(np.log10(wmin) + np.arange(nsub)*1.449E-6) * u.AA

    # Voigt for Lyman series
    tau_Lyman = lav.voigt_from_abslines(wave,HIlines,fwhm=0.,ret='tau')

    # Continuum opacity
    LL_comps = HI_comps['lgNHI'] > 15.0
    tau_LL = np.zeros(wave.size)
    for row in HI_comps[LL_comps]:
        # Energies in LLS rest-frame
        wv_rest = wave / (row['z']+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())
        # Get photo_cross and calculate tau
        itau_LL = (10.**row['lgNHI'] / u.cm**2) * photo_cross(1,1,energy)
        # Kludge around the limit
        if kludge:
            pix_LL = np.argmin(np.fabs(wv_rest- 911.3*u.AA))
            pix_kludge = np.where((wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA))[0]
            itau_LL[pix_kludge] = itau_LL[pix_LL]
        # Sum
        tau_LL += itau_LL
    # Total
    return wave, tau_LL + tau_Lyman


def mock_HIlines(HI_comps, wvmnx, tau0_min=5e-3):
    """Generate list of absorption lines for a sightline
    Parameters:
    ------------
    HI_comps : QTable
      HI components ('z','lgNHI','bval')
    wvmnx : tuple (Quantity,Quantity)
      min/max wavelength
    tau0_min : float, optional
      Minimum optical depth to include

    Returns:
    --------
    abs_lines: list
      List of HI Lyman lines
    """
    # Linelist
    HIlines = lll.LineList('HI')
    wrest = HIlines._data['wrest']

    # All wvobs
    grid_zp1 = np.outer(HI_comps['z']+1, np.ones(len(wrest)))
    grid_wr = np.outer(np.ones(len(HI_comps)), wrest)
    all_wvobs = grid_zp1 * grid_wr

    # All tau0
    # Grids
    grid_NHI = np.outer(10.**HI_comps['lgNHI'], np.ones(len(wrest)))
    grid_b = np.outer(HI_comps['bval'].to('cm/s'), np.ones(len(wrest)))
    grid_f = np.outer(np.ones(len(HI_comps)), HIlines._data['f'])
    # tau0 (Spitzer 1979)
    grid_tau0 = 1.497e-2 * grid_wr * 1e-8 * grid_f * grid_NHI / grid_b

    # Good ones
    gdlin = np.where( (all_wvobs>wvmnx[0].value) &
        (all_wvobs<wvmnx[1].value) & (grid_tau0>tau0_min))

    # Generate the line list for Lyman series
    grid_lgNHI = np.log10(grid_NHI)
    grid_z = grid_zp1 - 1.
    grid_bkms = grid_b / 1e5 * u.km/u.s

    # Parameters
    parm = []
    for ii,jj in zip(gdlin[0],gdlin[1]): # lambda, z, N, b
        parm.append([grid_wr[ii,jj]*u.AA, grid_z[ii,jj],
            grid_lgNHI[ii,jj], grid_bkms[ii,jj]])

    # Generate a large batch of AbsLines
    all_wrest = [iparm[0] for iparm in parm]
    abs_lines = spectralline.many_abslines(all_wrest, HIlines)
    # Fill in the rest
    for jj,iparm in enumerate(parm):
        abs_lines[jj].setz(iparm[1])
        abs_lines[jj].attrib['logN'] = iparm[2]
        abs_lines[jj].attrib['N'] = 10**iparm[2] / u.cm**2
        abs_lines[jj].attrib['b'] = iparm[3]

    # Retrun
    return abs_lines


def mk_mock(wave, zem, fN_model, out_spec=None, add_conti=True,
            out_tbl=None, s2n=15., fwhm=3., seed=None):
    """ Generate a mock
    Parameters
    ----------
    wave : Quantity array
    zem : float
    fN_model
    out_spec : str, optional
    out_tbl : str, optional
    s2n : float, optional
    fwhm : float, optional
    seed : int, optional

    Returns
    -------
    full_mock : XSpectrum1D of the mock
    HI_comps : Table of the components
    misc : tuple
      Other bits and pieces [may deprecate]
    """
    # Init
    rstate=np.random.RandomState(seed)

    # Wavelength Range
    wvmin = np.min(wave) #(1+zem)*912.*u.AA - 1000*u.AA
    wvmax = np.max(wave)
    dwv = np.median(wave - np.roll(wave,1))

    # zrange for Lya
    zmin = (wvmin/(1215.6700*u.AA))-1.

    # Components
    HI_comps = monte_HIcomp((zmin,zem), fN_model, rstate=rstate)

    # Main call
    HIlines = mock_HIlines(HI_comps, (wvmin,wvmax))

    # Optical depth
    sub_wave, tot_tau = generate_tau(wave, HIlines, HI_comps)
    dwv_sub = np.median(sub_wave-np.roll(sub_wave,1))

    # Normalized mock spectrum (sub-grid of wavelengths)
    sub_flux = np.exp(-1.*tot_tau)
    sub_mock = XSpectrum1D.from_tuple( (sub_wave,sub_flux) )

    # Smooth
    smooth_mock = sub_mock.gauss_smooth(fwhm=fwhm*(dwv/dwv_sub))

    # Rebin
    mock = smooth_mock.rebin(wave)

    # Add noise
    mock.sig = np.ones(len(mock.flux))/s2n
    noisy_mock = mock.add_noise(rstate=rstate)

    # Continuum
    if add_conti:
        conti, wfc3_idx = pycq.wfc3_continuum(zqso=zem,wave=wave,rstate=rstate)
        cflux = conti.flux
    else:
        cflux = np.ones_like(noisy_mock.flux)
        wfc3_idx = -1

    # Full
    full_mock = XSpectrum1D.from_tuple((wave,noisy_mock.flux*cflux,
                                 cflux*np.ones(len(noisy_mock.flux))/s2n))

    # Write spectrum (as desired)
    full_mock.meta.update(dict(zQSO=zem, S2N=s2n, seed=seed,
            fN_gamma=fN_model.gamma, fN_type=fN_model.fN_mtype, conti=wfc3_idx))
    if out_spec is not None:
        full_mock.write_to_fits(out_spec, clobber=True)

    # Save HI data (as desired)
    if out_tbl is not None:
        HI_comps.write(out_tbl, overwrite=True)
    # Return
    return full_mock, HI_comps, (sub_mock, tot_tau, mock)
