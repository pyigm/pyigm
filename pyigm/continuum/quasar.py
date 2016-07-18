""" Module for quasar continuum code
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import imp
import pdb

from scipy.interpolate import interp1d

from astropy import units as u
from astropy.io import fits, ascii
from astropy.table import Table
from astropy import constants as const

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.lists.linelist import LineList

from pyigm.fN import tau_eff as pyift
from pyigm.fN.fnmodel import FNModel


pyigm_path = imp.find_module('pyigm')[1]


def get_telfer_spec(zqso=0., igm=False, fN_gamma=None,
                    LL_flatten=True, nproc=6):
    """Generate a Telfer QSO composite spectrum

    Parameters
    ----------
    zqso : float, optional
      Redshift of the QSO
    igm : bool, optional
      Include IGM opacity? [False]
    fN_gamma : float, optional
      Power-law evolution in f(N,X)
    LL_flatten : bool, optional
      Set Telfer to a constant below the LL?
    nproc : int, optional
      For multi-processing with IGM

    Returns
    -------
    telfer_spec : XSpectrum1D
      Spectrum
    """
    # Read
    telfer = ascii.read(
        pyigm_path+'/data/quasar/telfer_hst_comp01_rq.ascii', comment='#')
    scale = telfer['flux'][(telfer['wrest'] == 1450.)]
    telfer_spec = XSpectrum1D.from_tuple((np.array(telfer['wrest'])*(1+zqso),
        np.array(telfer['flux'])/scale[0]))  # Observer frame

    # IGM?
    if igm is True:
        """The following concept is rather experimental.
        Use at your own risk.
        """
        import multiprocessing
        fN_model = FNModel.default_model()
        # Expanding range of zmnx (risky)
        fN_model.zmnx = (0.,5.5)
        if fN_gamma is not None:
            fN_model.gamma = fN_gamma
        # Setup inputs
        #EW_FIL = pyigm_path+'/data/fN/EW_SPLINE_b24.yml'
        #with open(EW_FIL, 'r') as infile:
        #    EW_spline = yaml.load(infile)  # dict from mk_ew_lyman_spline
        HI = LineList('HI')
        twrest = HI._data['wrest']
        # Parallel
        if LL_flatten:
            igm_wv = np.where((telfer['wrest'] > 900.) & (telfer['wrest'] < 1220.))[0]
        else:
            igm_wv = np.where(telfer['wrest'] < 1220.)[0]
        adict = []
        for wrest in telfer_spec.wavelength[igm_wv].value:
            tdict = dict(ilambda=wrest, zem=zqso, fN_model=fN_model,
                         wrest=twrest.copy())
            adict.append(tdict)
        # Run
        if nproc > 1:
            pool = multiprocessing.Pool(nproc) # initialize thread pool N threads
            ateff = pool.map(pyift.map_lymanew, adict)
        else:
            ateff = map(pyift.map_lymanew, adict)
        # Apply
        new_flux = telfer_spec.flux.value
        new_flux[igm_wv] *= np.exp(-1.*np.array(ateff))
        # Flatten?
        if LL_flatten:
            wv_LL = np.where(np.abs(telfer_spec.wavelength/(1+zqso)-914.*u.AA)<3.*u.AA)[0]
            f_LL = np.median(new_flux[wv_LL])
            wv_low = np.where(telfer_spec.wavelength/(1+zqso)<911.7*u.AA)[0]
            new_flux[wv_low] = f_LL
        # Regenerate spectrum
        telfer_spec = XSpectrum1D.from_tuple(
                (np.array(telfer['wrest'])*(1+zqso), new_flux))

    # Return
    return telfer_spec

def full_nu_spectrum(alpha_euv=1.7):
    """Modified version of the IDL code by JFH

    Parameters
    ----------
    alpha_euv: float, optional
      Power-law in EUV

    Returns:
    ----------
    lognu : ndarray
      log10 frequency of evaluated SED
    fnu_qso : ndarray
      log10 fnu of the QSO
    """

    clight = const.c.cgs
    # Beta spliced to vanden Berk template with host galaxy  removed
    van_file = pyigm_path+'/data/quasar/VanDmeetBeta_fullResolution.txt'
    van_tbl = Table.read(van_file,format='ascii')
    isort = np.argsort(van_tbl['nu'])
    nu_van = van_tbl['nu'][isort]
    fnu_van = van_tbl['f_nu'][isort]/1e-17
    #
    lam_van = (clight/(nu_van*u.Hz)).to('AA')
    lognu_van = np.log10(nu_van)
    logfnu_van = np.log10(fnu_van)
    # Beta spliced to vanden Berk template with host galaxy  remove
    gtr_file = pyigm_path+'/data/quasar/richards_2006_sed.txt'
    lognu_gtr, logL_gtr  = np.loadtxt(gtr_file, skiprows=29, usecols=(0,1), unpack=True)
    # Sort
    isort = np.argsort(lognu_gtr)
    lognu_gtr = lognu_gtr[isort]
    logL_gtr = logL_gtr[isort]
    #
    logLnu_gtr = logL_gtr - lognu_gtr  # Lnu units erg/s/Hz
    # create a vector of frequencies in log units
    lognu_min = np.min(lognu_gtr)  # This is about 1e-3 Rydberg and limit of richards template
    lognu_max = (np.log10(1e5*const.Ryd*clight/u.Hz)).value  # 1.36 Mev
    # SDSS pixel scale
    dlognu = 0.0001
    n_nu = int(np.round((lognu_max - lognu_min)/dlognu))
    lognu = np.arange(n_nu)*dlognu + lognu_min
    nu = 10.0**lognu
    logfnu = np.zeros(n_nu)
    #;; Some relevant frequencies
    lognu_2500 = np.log10((clight/2500./u.AA).to('Hz').value) #;; 2500A for alpha_OX
    lognu_8000 = np.log10((clight/8000./u.AA).to('Hz').value) #;; IR splice is at 8000A
    lognu_1Ryd = np.log10((const.Ryd*clight).to('Hz').value) #
    lognu_30Ryd = np.log10(30.0) + lognu_1Ryd
    lognu_2keV = np.log10((2000.0*u.eV/const.h).to('Hz').value)
    lognu_100keV = np.log10(50.) + lognu_2keV
    i8000 = (lam_van >= 8000.0*u.AA) & (lam_van <= 8100.0*u.AA)
    #;; median as template is noisy here
    logfnu_van_8000 = np.median(logfnu_van[i8000])
    logLnu_gtr_8000 = float(interp1d(lognu_gtr, logLnu_gtr)(lognu_8000))
    #logLnu_gtr_2500 = float(interp1d(lognu_gtr, logLnu_gtr)(lognu_2500))
    #;; IR part: nu < 8000A/c;  use the Richards al. 2006 template
    i_8000 = lognu <= lognu_8000
    logfnu[i_8000] = interp1d(lognu_gtr, logLnu_gtr)(lognu[i_8000])
    #;; UV part: c/8000A < nu < 1 Ryd/h; use the template itself
    i_UV = (lognu > lognu_8000) & (lognu <= lognu_1Ryd)
    logfnu[i_UV] = (logLnu_gtr_8000 - logfnu_van_8000 +
                interp1d(lognu_van, logfnu_van)(lognu[i_UV]))
    logfnu_van_1Ryd =  (logLnu_gtr_8000 - logfnu_van_8000 +
                    interp1d(lognu_van, logfnu_van)(lognu_1Ryd))
    logfnu_2500 = interp1d(lognu, logfnu)(lognu_2500)
    '''
    ;; This equation is from Strateva et al. 2005 alpha_OX paper.  I am
    ;; evaluating the alpha_OX at the L_nu of the template, which is
    ;; based on the normalization of the Richards template. A more
    ;; careful treatment would actually use the 2500A luminosity of the
    ;; quasar template, after it is appropriately normalized
    '''
    alpha_OX = -0.136*logfnu_2500 + 2.630
    logfnu_2keV = logfnu_2500 + alpha_OX*(lognu_2keV - lognu_2500)
    #;; FUV par 1 Ryd/h < nu < 30 Ryd/h;
    #;; use the alpha_EUV power law
    i_FUV = (lognu > lognu_1Ryd) &  (lognu <= lognu_30Ryd)
    logfnu[i_FUV] = logfnu_van_1Ryd - alpha_euv*(lognu[i_FUV] - lognu_1Ryd)
    logfnu_30Ryd  =  logfnu_van_1Ryd - alpha_euv*(lognu_30Ryd - lognu_1Ryd)
    '''
    ;; soft X-ray part 30 Ryd/h < nu < 2kev/h;
    ;; use a power law with a slope alpha_soft chosen to match the
    ;; fnu_2Kev implied by the alpha_OX
    '''
    i_soft = (lognu > lognu_30Ryd) & (lognu <= lognu_2keV)
    alpha_soft = (logfnu_2keV - logfnu_30Ryd)/(lognu_2keV - lognu_30Ryd)
    logfnu[i_soft] = logfnu_30Ryd + alpha_soft*(lognu[i_soft] - lognu_30Ryd)
    #;; X-ray part 2 kev/h < nu < 100 keV/h
    i_X = (lognu > lognu_2keV) & (lognu <= lognu_100keV)
    alpha_X = -1.0 #;; adopt this canonical 'flat X-ray slope'
    logfnu[i_X] = logfnu_2keV + alpha_X*(lognu[i_X] - lognu_2keV)
    logfnu_100keV = logfnu_2keV + alpha_X*(lognu_100keV - lognu_2keV)
    #;; hard X-ray part nu > 100 keV/h
    i_HX = lognu > lognu_100keV
    alpha_HX = -2.0 #;; adopt this canonical 'flat X-ray slope'
    logfnu[i_HX] = logfnu_100keV + alpha_HX*(lognu[i_HX] - lognu_100keV)

    #;; subtract off 10^30 to make units manageable
    fnu_qso = 10.0**(logfnu - 30.0)

    # Return
    return lognu, fnu_qso

def wfc3_continuum(wfc3_indx=None, zqso=0., wave=None,
                   smooth=3., NHI_max=17.5, rstate=None):
    """Use the WFC3 data + models from O'Meara+13 to generate a continuum

    Parameters
    ----------
    wfc3_indx : int, optional
      Index of WFC3 data to use
    zqso : float, optional
      Redshift of the QSO
    wave : Quantity array, optional
      Wavelengths to rebin on
    smooth : float, optional
      Number of pixels to smooth on
    NHI_max : float, optional
      Maximum NHI for the sightline
    rstate : RandomState, optional
      Useful for generating random mocks that you can re-generate

    Returns
    -------
    wfc3_continuum : XSpectrum1D 
       of the continuum
    idx : int
      Index of the WFC3 spectrum used    
    """
    # Random number
    if rstate is None:
        rstate = np.random.RandomState()
    # Open
    wfc3_models_hdu = fits.open(pyigm_path+'/data/quasar/wfc3_conti_models.fits.gz')
    nwfc3 = len(wfc3_models_hdu)-1
    # Load up models
    wfc_models = []
    for ii in range(1,nwfc3):
        wfc_models.append( Table(wfc3_models_hdu[ii].data) )
    # Grab a random one
    if wfc3_indx is None:
        need_c = True
        while need_c:
            idx = rstate.randint(0,nwfc3-2)
            if wfc_models[idx]['TOTNHI'] > NHI_max:
                continue
            if wfc_models[idx]['QSO'] in ['J122836.05+510746.2',
                                          'J122015.50+460802.4']:
                continue  # These QSOs are NG
            need_c = False
    else:
        idx = wfc3_indx

    # Generate spectrum
    wfc_spec = XSpectrum1D.from_tuple((
        wfc_models[idx]['WREST'].data.flatten()*(1+zqso),
        wfc_models[idx]['FLUX'].data.flatten()))
    # Smooth
    wfc_smooth = wfc_spec.gauss_smooth(fwhm=smooth)

    # Rebin?
    if wave is not None:
        wfc_rebin = wfc_smooth.rebin(wave)
        return wfc_rebin, idx
    else:
        return wfc_smooth, idx
