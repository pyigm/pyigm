""" Utilities for IGMSystem and IGMSurvey
Best to keep these separate from the Class modules
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import warnings
import pdb
import numpy as np
from scipy.interpolate import interp1d

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from astropy import units as u
from astropy.table import Table, Column

from linetools.analysis import absline as ltaa
from linetools.spectra.xspectrum1d import XSpectrum1D


def plot_contours(lnL3D, xtup):
    """
    Parameters
    ----------
    lnL3D
    C0_val
    C1_val
    NHI_eval

    Returns
    -------

    """
    from scipy import stats as scistats
    if len(xtup) == 3:
        NHI_eval, C0_val, C1_val = xtup

    # Best values
    cl2D = cl_interval(lnL3D)
    # CL image
    cl_img_C = cl_image(lnL3D[cl2D[0][0],:,:])
    cl_img_C = cl_img_C.T

    # Begin plot
    plt.figure(figsize=(8, 5))
    plt.clf()
    gs = gridspec.GridSpec(1, 2)

    # Continuum
    ax = plt.subplot(gs[0])
    ax.set_xlabel(r'$C_0 \; \rm (erg/s/cm^2/\AA)$ ', size=15.)
    ax.set_ylabel(r'$C_1 \; \rm (erg/s/cm^2/\AA^2)$', size=15.)

    extent=[C0_val[0], C0_val[-1], C1_val[0], C1_val[-1], ]
    aspect=np.abs((C0_val[-1]-C0_val[0])/(C1_val[-1]-C1_val[0]))
    cax = ax.imshow(cl_img_C, cmap=plt.cm.gist_heat, interpolation='none',
        extent=extent, origin='lower', aspect=aspect)
    # Contour levels
    lvls = [1.-2*(1.-scistats.norm.cdf(ii+1)) for ii in range(3)]
    CS = ax.contour(cl_img_C, lvls, colors='k', linewidths=1.4, hold='on', extent=extent)
    plt.clabel(CS, CS.levels, inline=True, fmt='%0.3f', fontsize=10)

    # NHI
    # Confidence level image
    cl_img_NC = cl_image(lnL3D[:,:,cl2D[0][2]])
    ax = plt.subplot(gs[1])
    ax.set_xlabel(r'$C_0 \; \rm (erg/s/cm^2/\AA)$ ', size=15.)
    ax.set_ylabel(r'$\log \, N_{\rm HI}$', size=15.)
    #
    extent=[C0_val[0], C0_val[-1], NHI_eval[0], NHI_eval[-1], ]
    aspect=np.abs((C0_val[-1]-C0_val[0])/(NHI_eval[-1]-NHI_eval[0]))
    cax = ax.imshow(cl_img_NC, cmap=plt.cm.gist_heat, interpolation='none',
        extent=extent, origin='lower', aspect=aspect)
    # Contour levels
    lvls = [1.-2*(1.-scistats.norm.cdf(ii+1)) for ii in range(3)]
    CS2 = ax.contour(cl_img_NC, lvls, colors='k', linewidths=1.4, hold='on', extent=extent)
    ax.clabel(CS2, CS2.levels, inline=True, fmt='%0.3f', fontsize=10)


    # End
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    plt.show()


def plot_NHI_model(lls_dict, ax, lsz=12., touch=False, scl=1., csz=10.):
    """  Plot the NHI model on the spectrum

    Parameters
    ----------
    lls_dict : dict
      See maxlik_linearc for a description
    ax : matplotlib.axis
    lsz: float, optional
      Font size for axes labels
    touch : bool, optional
    scl : float, optional
    csz : float, optional

    Returns
    -------

    """
    from linetools.spectra.plotting import get_flux_plotrange

    spec, xspec, gdp, NHI, tau0 = setup_lls_fit_analy(lls_dict['spec_fil'], lls_dict['z'], lls_dict['windows'], lls_dict['NHI_mnx'])
    # Scale
    xspec.data['flux'] *= scl
    # Limits
    xmnx = [lls_dict['windows'][0][0], 940.*(1+lls_dict['z'])]
    if lls_dict['cdict']['type'] == 'Gaussian':
        ymnx = [-1*lls_dict['cdict']['sig'], lls_dict['cdict']['best']+4*lls_dict['cdict']['sig']]
    elif lls_dict['cdict']['type'] == 'Fixed':
        ymnx = [-0.1*lls_dict['cdict']['value'], 1.5*lls_dict['cdict']['value']]
    elif lls_dict['cdict']['type'] == 'Fit_const':
        ymnx = [-1*(lls_dict['cdict']['fit_val'][0]-lls_dict['cdict']['fit_val'][1]),
                3*(lls_dict['cdict']['fit_val'][2]-lls_dict['cdict']['fit_val'][0])+
                lls_dict['cdict']['fit_val'][0]]
    elif lls_dict['cdict']['type'] == 'Fit_line':
        if gdp is None:
            gdp = (xspec.wavelength>xmnx[0]*u.AA) & (xspec.wavelength<xmnx[1]*u.AA)
        conti = lls_dict['cdict']['best'][0] + lls_dict['cdict']['best'][1]*(
            xspec.wavelength[gdp].value-lls_dict['cdict']['slope_pivot']*(1+lls_dict['z']))
        mx = np.max(conti)
        ymnx = [-0.1*mx, mx*1.3]
    else:
        raise ValueError("Need to setup this continuum model")
    # Extend xmnx
    if lls_dict['cdict']['type'] in ['Fit_line', 'Fit_const']:
        xmx = 0.
        for rng in lls_dict['cdict']['analy']:
            xmx = max(xmx, rng[1])
        xmnx[1] = xmx+3.
    # Scale
    ymnx = np.array(ymnx)*scl
    # Finally
    #idx = (xspec.wavelength > xmnx[0]*u.AA) & (xspec.wavelength < xmnx[1]*u.AA)
    idx = gdp
    f_ymnx = get_flux_plotrange(xspec.flux[idx].value)
    ymnx[1] = max(ymnx[1],f_ymnx[1])


    # Axes
    #ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_major_locator(plt.MultipleLocator(20.))
    #ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    #ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.set_xlim(xmnx)
    ax.set_ylim(ymnx)
    if scl == 1.:
        ax.set_ylabel(r'$f_\lambda$ (cgs)', size=lsz)
    else:
        ax.set_ylabel(r'$f_\lambda$ ($10^{-15}$ cgs)', size=lsz)


    # Plot data
    ax.plot(xspec.wavelength, xspec.flux, color='black', drawstyle='steps-mid',
            zorder=2)
    try:
        ax.plot(xspec.wavelength, scl*xspec.sig, ':', color='red', zorder=1)
    except ValueError:
        pdb.set_trace()

    # Binned
    if False:
        binsz = 5.
        binwv = np.arange(1040., 1200., binsz)*u.AA
        binspec = xspec.rebin(binwv)
        gdp = binspec.wavelength.value < 910.*(1+lls_dict['z'])
        ax.scatter(binspec.wavelength.value[gdp]+binsz/2., binspec.flux[gdp],
                   color='yellow', zorder=300)
                   #edgecolor='none')#, alpha=0.5)

    # Best continuum
    if lls_dict['cdict']['type'] == 'Gaussian':
        conti = lls_dict['cdict']['best']*np.ones_like(xspec.flux.value)
    elif lls_dict['cdict']['type'] == 'Fixed':
        conti = lls_dict['cdict']['value']*np.ones_like(xspec.flux.value)
    elif lls_dict['cdict']['type'] == 'Fit_const':
        conti = lls_dict['cdict']['fit_val'][0]*np.ones_like(xspec.flux.value)
    elif lls_dict['cdict']['type'] == 'Fit_line':
        conti = lls_dict['cdict']['best'][0] + lls_dict['cdict']['best'][1]*(
            xspec.wavelength.value-lls_dict['cdict']['slope_pivot']*(1+lls_dict['z']))
    ax.plot(xspec.wavelength, conti*scl, '--', color='green', zorder=3)
    ax.minorticks_on()
    if touch is True:
        pass
        #ax.get_xaxis().set_ticks([])
    else:
        ax.set_xlabel('Wavelength (Ang)', size=lsz)

    # Best Model
    mclr = 'lightblue'
    wv_rest = xspec.wavelength / (lls_dict['z']+1)
    energy = wv_rest.to(u.eV, equivalencies=u.spectral())
    tau0 = (10.**lls_dict['fit_NHI'][0] / u.cm**2) * ltaa.photo_cross(1, 1, energy)
    if lls_dict['analy_type'] in ['Fit_Conti', 'Vary_Conti']:
        if lls_dict['fit_NHI'][0] != lls_dict['fit_NHI'][2]:
            best_model = scl*conti * np.exp(-1*tau0)
            abs = tau0 > 0.
            ax.plot(xspec.wavelength[abs], best_model[abs], color=mclr, zorder=100)

    # Continuum Error
    clr_ce = 'lightgreen'
    alpha_ce = 0.4
    if lls_dict['cdict']['type'] == 'Gaussian':
        cwv = tau0 == 0.
        npix = np.sum(cwv)
        ax.fill_between(xspec.wavelength.value[cwv],
                        [scl*lls_dict['cdict']['best']+lls_dict['cdict']['sig']]*npix,
                        [scl*lls_dict['cdict']['best']-lls_dict['cdict']['sig']]*npix,
                        color=clr_ce, alpha=alpha_ce, zorder=50)
    elif lls_dict['cdict']['type'] == 'Fit_const':
        for rng in lls_dict['cdict']['analy']:
            idx = ((xspec.wavelength > rng[0]*u.AA) &
                   (xspec.wavelength < rng[1]*u.AA) &
                   (xspec.sig > 0))
            gdC = np.where(idx)[0]
            ax.fill_between(xspec.wavelength.value[gdC],
                        [scl*lls_dict['cdict']['fit_val'][1]]*gdC.size,
                        [scl*lls_dict['cdict']['fit_val'][2]]*gdC.size,
                        color=clr_ce, alpha=alpha_ce, zorder=50)
    elif lls_dict['cdict']['type'] == 'Fit_line':
            #xdb.set_trace()
            if 'fit_val' in lls_dict['cdict']:
                for rng in lls_dict['cdict']['analy']:
                    idx = ((xspec.wavelength > rng[0]*u.AA) &
                           (xspec.wavelength < rng[1]*u.AA) &
                           (xspec.sig > 0))
                    gdC = np.where(idx)[0]
                #
                    sig0 = (lls_dict['cdict']['fit_val'][0][2]-lls_dict['cdict']['fit_val'][0][1])/2.
                    sig1 = (lls_dict['cdict']['fit_val'][1][2]-lls_dict['cdict']['fit_val'][1][1])/2.
                    sigl = np.sqrt(sig0**2 +
                               sig1**2*(lls_dict['cdict']['slope_pivot']*(1+lls_dict['z'])-
                                        xspec.wavelength.value[gdC])**2)
                    ax.fill_between(xspec.wavelength.value[gdC],
                                scl*(conti[gdC] + sigl),
                                scl*(conti[gdC] - sigl),
                                color=clr_ce, alpha=alpha_ce, zorder=50)

    # Model with error (limits too)
    taulow = (10.**lls_dict['fit_NHI'][1] / u.cm**2) * ltaa.photo_cross(1, 1, energy)
    try:
        low_model = scl*conti * np.exp(-1*taulow)
    except ValueError:
        pdb.set_trace()
    tauhi = (10.**lls_dict['fit_NHI'][2] / u.cm**2) * ltaa.photo_cross(1, 1, energy)
    hi_model = scl*conti * np.exp(-1*tauhi)
    mwv = tau0 > 0.
    ax.fill_between(xspec.wavelength.value[mwv], low_model[mwv],
                    hi_model[mwv], color=mclr, alpha=0.3, zorder=100)

    # Finish
    ax.plot(xmnx, [0.,0.], '--', color='gray')


def setup_lls_fit_analy(spec_fil, zlls, lls_windows, NHI_mnx, nNHI=100, spec_keys=None):
    """ Setup for LLS fit analysis.  This is pretty
    specific to COS spectra

    Parameters
    ----------
    spec_fil : str
      Path+name of spectrum file
    zlls : float
      Redshift of the LLS
    lls_windows : list
      List of lists of float specifying windows to perform LLS analysis e.g.  [ [1323.3, 1355.3], [1360.3, 1366.1]]
      Allows user to skip over unwanted features or data (e.g. chip gap)
    NHI_mnx : list or tuple
      log10 NHI min and max
    nNHI : int, optional
      Number of evaluations
    spec_keys : dict, optional
      specifies tags for sig, wave, flux (defaults are 'ERROR', 'WAVE', 'FLUX')

    Returns
    -------
    spec : Table
    xspec : XSpectrum1D
    NHI : ndarray
      NHI values for evaluation
    tau0 : ndarray
      Optical depth in analysis range

    """
    # Init
    if spec_keys is None:
        spec_keys = dict(sig='ERROR', flux='FLUX', wave='WAVE')
    # Load up spectrum (Table and xspec)
    spec = Table.read(spec_fil)
    # Deal with NANs
    sig = spec[spec_keys['sig']].data.flatten()
    sig[np.isnan(sig)] = 0.
    xspec = XSpectrum1D.from_tuple((np.array(spec[spec_keys['wave']].data.flatten()),
                                    np.array(spec[spec_keys['flux']].data.flatten()),
                                    sig), masking='none')

    # Analysis pixels
    pixels = []
    for window in lls_windows:
        gdwv = np.where((xspec.wavelength >= window[0]*u.AA) &
                    (xspec.wavelength <= window[1]*u.AA))[0]
        pixels.append(gdwv)
    gdwv = np.concatenate(pixels)

    # NHI
    NHI = np.linspace(NHI_mnx[0], NHI_mnx[1], num=nNHI)
    wv_rest = xspec.wavelength[gdwv] / (zlls+1)
    energy = wv_rest.to(u.eV, equivalencies=u.spectral())
    # Get photo_cross and calculate tau
    tau0 = (10.**NHI[0] / u.cm**2) * ltaa.photo_cross(1, 1, energy)

    # Return
    return spec, xspec, gdwv, NHI, tau0


def maxlik_fitlinearc(lls_dict, neval=100, nevalC=50, slope_pivot=911.*u.AA, **kwargs):
    """ Max Likelihood analysis of a Lyman limit with fit to a linear continuum over a given range
    Uses ~25Gb with neval = 100

    Parameters
    ----------
    lls_dict : dict
      Contains all relevant input values for the analysis
      Here are the required keys:
      'spec_fil' : str -- Path+filename of the spectrum (a COS spectrum in binary FITS table)
      'z': float -- Redshift of the LLS
      'windows' : list
        List of lists of float specifying windows to perform LLS analysis e.g.  [ [1323.3, 1355.3], [1360.3, 1366.1]]
        Allows user to skip over unwanted features or data (e.g. chip gap)
       'NHI_mnx' : tuple or list -- lower/upper bound on NHI values for analysis
       'cdict' : dict
         dict describing the continuum with keys:
           'analy_type' : str  -- must be 'Fit_line'
           'C0_range' : list -- low/high range for fitting continuum normalization
           'C1_range' : list -- low/high range for fitting continuum slope
           'analy' : list of lists -- list of windows for fitting continuum, e.g.
               [[1303.3, 1303.9], [1298.0, 1299.2], [1296.2,1296.9]]
    neval : int, optional
      Size of grid
    slope_pivot: Quantity, optional
      Defines the line continuum model:   conti = C0 + C1 * (lambda - slope_pivot)

    Returns
    -------
    NHI : ndarray
      Array of NHI values considered
    C0_val : ndarray
      Array of C0 values considered
    C1_val : ndarray
      Array of C1 values considered
    lnL : ndarray
      3D array of likelihood grid corresponding to NHI, C0, C1 values

    """
    from scipy.special import gammaln
    # Setup
    spec, xspec, gdwv, NHI, tau0 = setup_lls_fit_analy(lls_dict['spec_fil'], lls_dict['z'], lls_dict['windows'],
                                                       lls_dict['NHI_mnx'], nNHI=neval, **kwargs)
    cdict = lls_dict['cdict']
    npix = gdwv.size

    # Check continuum
    if cdict['type'] == 'Fit_line':
        gdC = []
        for rng in cdict['analy']:
            idx = ((xspec.wavelength > rng[0]*u.AA) &
                   (xspec.wavelength < rng[1]*u.AA) &
                   (xspec.sig > 0))
            gdC = gdC + list(np.where(idx)[0])
        gdC = np.array(gdC)
        # Grab continuum
        conti = xspec.flux[gdC]
        sig_conti = xspec.sig[gdC]
        # Set continuum range
        medC = np.median(conti.value)
        C0_val = np.linspace(cdict['C0_range'][0], cdict['C0_range'][1], num=nevalC)
        C1_val = np.linspace(cdict['C1_range'][0], cdict['C1_range'][1], num=nevalC)
        cdict['slope_pivot'] = slope_pivot.value
    else:
        raise ValueError('Not ready for this type of continuum')

    # Generate arrays (value x NHI x Cval)
    # Arrays -- The following is for HSLA outputs
    count_array = np.outer(np.round(spec['GROSSCOUNTS']*spec['EXP_PIX']).data.flatten()[gdwv], np.ones(neval))
    dark_dumb = np.median(((spec['GROSSCOUNTS']-spec['NETCOUNTS'])*spec['EXP_PIX']).data.flatten()[gdwv])
    dark_array = np.outer(dark_dumb*np.ones(len(gdwv)), np.ones(neval))
    calib_array = np.outer(1./spec['FLUXFACTOR'].data.flatten()[gdwv], np.ones(neval))
    expt_array = np.outer(spec['EXP_PIX'].data.flatten()[gdwv], np.ones(neval))
    #pdb.set_trace()
    # Continuum
    obs_conti_array = np.outer(conti.value, np.ones(nevalC))
    sig_conti_array = np.outer(sig_conti, np.ones(nevalC))
    # Grids
    count_grid = np.zeros((npix, neval, nevalC, nevalC))
    dark_grid = np.zeros((npix, neval, nevalC, nevalC))
    calib_grid = np.zeros((npix, neval, nevalC, nevalC))
    expt_grid = np.zeros((npix, neval, nevalC, nevalC))
    # Fill
    for ii in range(nevalC):
        for jj in range(nevalC):
            count_grid[:, :, ii, jj] = count_array
            dark_grid[:, :, ii, jj] = dark_array
            calib_grid[:, :, ii, jj] = calib_array
            expt_grid[:, :, ii, jj] = expt_array
    # Model Continuum in LL
    LL_conti_sub_grid = np.zeros((npix, nevalC, nevalC))
    LL_conti_grid = np.zeros((npix, neval, nevalC, nevalC))
    LL_conti_array = np.outer((xspec.wavelength[gdwv]-
                               slope_pivot*(1+lls_dict['z'])), C1_val)
    for jj,C0 in enumerate(C0_val):
        LL_conti_sub_grid[:,jj,:] = LL_conti_array + C0
    for jj in range(neval):
        LL_conti_grid[:,jj,:,:] = LL_conti_sub_grid

    # Model Continuum redward
    conti_array = np.outer((xspec.wavelength[gdC]-slope_pivot*(1+lls_dict['z'])), C1_val)
    model_conti_grid = np.zeros((gdC.size, nevalC, nevalC))
    for ii,C0 in enumerate(C0_val):
        model_conti_grid[:, ii, :] = conti_array + C0

    # Observed + error (Continuum points only)
    obs_conti_grid = np.zeros((gdC.size, nevalC, nevalC))
    sig_conti_grid = np.zeros((gdC.size, nevalC, nevalC))
    for ii in range(nevalC):
        obs_conti_grid[:, ii, :] = obs_conti_array
        sig_conti_grid[:, ii, :] = sig_conti_array

    # tau
    tau_array = np.zeros( (npix, neval) )
    for kk, iNHI in enumerate(NHI):
        tau_array[:, kk] = tau0 * 10.**(iNHI-NHI[0])
    tau_grid = np.zeros((npix, neval, nevalC, nevalC))
    for ii in range(nevalC):
        for jj in range(nevalC):
            tau_grid[:, :, ii, jj] = tau_array

    # Giddy up -- Likelihood
    model_flux = LL_conti_grid * np.exp(-1 * tau_grid)
    model_flux = np.maximum(model_flux, 0.)
    model_counts = model_flux*calib_grid*expt_grid + dark_grid
    lnP_LL = -1*model_counts + count_grid * np.log(model_counts) - gammaln(count_grid+1)
    lnP_C = -1*(obs_conti_grid-model_conti_grid)**2 / 2 / (sig_conti_grid**2)

    #pdb.set_trace()
    # Sum
    sum_LL = np.sum(lnP_LL, axis=0)
    sum_C = np.sum(lnP_C, axis=0)
    lnL = sum_LL
    for ii in range(neval):
        lnL[ii,:,:] += sum_C

    # Free up memory
    del model_counts, model_flux, LL_conti_grid, tau_grid, lnP_LL, expt_grid, calib_grid

    # Return
    return NHI, C0_val, C1_val, lnL


def analyze_lnl(lls_dict, xtup, lnL, CL=0.68):
    """ Analyze an input likelihood grid and fill best answers
    into input dict

    Parameters
    ----------
    lls_dict : dict
    xtup : tuple
      Arrays of values for the likelihood grid
    lnL : ndarray

    Returns
    -------
    lls_dict is modified in place
      'fit_NHI' -- best, lower bound, upper bound on NHI
      'cdict/fit_val' -- best, lower bound, upper bound on C0, and C1
    """
    if len(xtup) == 3:
        NHI_eval, C0_val, C1_val = xtup
    cl3D = cl_interval(lnL, CL=CL)  # 1-sigma!
    lls_dict['fit_NHI'] = [NHI_eval[cl3D[0][0]], NHI_eval[cl3D[1][0][0]], NHI_eval[cl3D[1][0][1]]]
    lls_dict['cdict']['best'] = (float(C0_val[cl3D[0][1]]), float(C1_val[cl3D[0][2]]))
    lls_dict['cdict']['fit_val'] = [[float(C0_val[cl3D[0][1]]),
                            float(C0_val[cl3D[1][1][0]]),
                            float(C0_val[cl3D[1][1][1]])],
                           [float(C1_val[cl3D[0][2]]),
                            float(C1_val[cl3D[1][2][0]]),
                            float(C1_val[cl3D[1][2][1]])]]
    # Return
    return


def cl_interval(lnL, sigma=None, CL=0.68):
    """ Calculate a confidence level interval from a log-likelihood image
    Simple area under the curve with the image collapsed along each dimension
    Taken from xastropy.stats.likelihood

    Parameters:
      lnL: np.array
        log-Likelihood image
      CL: float, optional
      sigma: float, optional
        Use to calculate confindence interval

    Returns:
      best_idx, all_error: Lists
        [best] [-, +] indices for each dimension
    """
    import copy
    # Confidence limits
    if sigma is None:
        c0 = (1. - CL)/2.
        c1 = 1.-c0
    # Image dimensions
    shape = lnL.shape
    ndim = len(shape)
    slc = [slice(None)]*ndim
    # Find max
    norm_L = np.exp(np.maximum(lnL - np.max(lnL),-15.))
    # Find best indices
    indices = np.where(lnL == np.max(lnL))
    best_idx = [bi[0] for bi in indices]

    # Error intervals
    all_error = []
    for kk in range(ndim):
        # Collapse on this dimension
        slc = copy.deepcopy(best_idx)
        slc[kk] = slice(None)
        Lslice = norm_L[slc].flatten()
        # Interpolate and go
        cumul_area = np.cumsum(Lslice)
        f_area = interp1d(cumul_area/cumul_area[-1], np.arange(len(Lslice)))
        # Here we go
        idx0 = int(np.round(f_area(c0)))
        idx1 = int(np.round(f_area(c1)))
        all_error.append([idx0,idx1])

    # Return
    return best_idx, all_error

def cl_image(lnL, sigma=False):
    """ Calculate a confidence level image from a lnL image
    Simple area under the curve with cubic spline interpolation

    Parameters:
      lnL: np.array
        log-Likelihood image
        Should probably be 2D
      sigma: bool, optional
        Return as sigma values [not implemented]

    Returns:
      cl_img: np.array
        Image with the same dimensions with confidence levels
    """
    # Max
    mxL = np.max(lnL)

    # Noramlize and flatten
    norm_img = lnL-mxL
    flat_img = norm_img.flatten()

    # Sort
    srt = np.argsort(flat_img)
    norm_lnL = flat_img[srt]

    # Sum
    cumul_area = np.cumsum(np.exp(np.maximum(norm_lnL,-15.)))
    tot_area = np.max(cumul_area)
    cumul_area = cumul_area/tot_area

    # Interpolation (smoothing a bit)
    f_area = interp1d(norm_lnL, cumul_area)

    # Finish
    area_img = f_area(norm_img)
    cl_img = 1.-area_img

    # Return
    return cl_img