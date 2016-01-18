""" Module for quasar continuum code
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import imp
import pdb
import yaml
import copy

from astropy import units as u
from astropy.io import fits, ascii
from astropy.table import Table

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
    telfer_spec = XSpectrum1D.from_tuple((telfer['wrest']*(1+zqso),
        telfer['flux']/scale[0]))  # Observer frame

    # IGM?
    if igm is True:
        """The following concept is rather experimental.
        Use at your own risk.
        """
        import multiprocessing
        fN_model = FNModel.default_model()
        # Expanding range of zmnx (risky)
        fN_model.zmnx = (0.,5.)
        if fN_gamma is not None:
            fN_model.gamma = fN_gamma
        # Setup inputs
        #EW_FIL = pyigm_path+'/data/fN/EW_SPLINE_b24.yml'
        #with open(EW_FIL, 'r') as infile:
        #    EW_spline = yaml.load(infile)  # dict from mk_ew_lyman_spline
        HI = LineList('HI')
        twrest = HI._data['wrest']
        # Parallel
        igm_wv = np.where(telfer['wrest'] < 1220.)[0]
        adict = []
        for wrest in telfer_spec.dispersion[igm_wv].value:
            tdict = dict(ilambda=wrest, zem=zqso, fN_model=fN_model,
                         wrest=copy.deepcopy(twrest))
            adict.append(tdict)
        # Run
        if nproc > 1:
            pool = multiprocessing.Pool(nproc) # initialize thread pool N threads
            ateff = pool.map(pyift.map_lymanew, adict)
        else:
            ateff = map(pyift.map_lymanew, adict)
        # Apply
        telfer_spec.flux[igm_wv] *= np.exp(-1.*np.array(ateff))
        # Flatten?
        if LL_flatten:
            wv_LL = np.where(np.abs(telfer_spec.dispersion/(1+zqso)-914.*u.AA)<3.*u.AA)[0]
            f_LL = np.median(telfer_spec.flux[wv_LL])
            wv_low = np.where(telfer_spec.dispersion/(1+zqso)<911.7*u.AA)[0]
            telfer_spec.flux[wv_low] = f_LL

    # Return
    return telfer_spec


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
    for ii in range(1,nwfc3-1):
        wfc_models.append( Table(wfc3_models_hdu[ii].data) )
    # Grab a random one
    if wfc3_indx is None:
        need_c = True
        while need_c:
            idx = rstate.randint(0,nwfc3-1)
            if wfc_models[idx]['TOTNHI'] > NHI_max:
                continue
            if wfc_models[idx]['QSO'] in ['J122836.05+510746.2',
                                          'J122015.50+460802.4']:
                continue  # These QSOs are NG
            need_c = False
    else:
        idx = wfc3_indx

    # Generate spectrum
    wfc_spec = XSpectrum1D.from_tuple((wfc_models[idx]['WREST'].flatten()*(1+zqso),
        wfc_models[idx]['FLUX'].flatten()))
    # Smooth
    wfc_smooth = wfc_spec.gauss_smooth(fwhm=smooth)

    # Rebin?
    if wave is not None:
        wfc_rebin = wfc_smooth.rebin(wave)
        return wfc_rebin, idx
    else:
        return wfc_smooth, idx
