""" Module for core continuum codes
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import numpy.polynomial.legendre as L
from numpy.linalg import svd
from astropy import units as u
from scipy import stats



def init_conti_dict(Norm=0., tilt=0., tilt2=0.,
                    piv_wv=0., piv_wv2=None, igm='None',
                    fN_gamma=-1., LL_flatten='True'):
    """Initialize a continuum conti_dict

    Parameters
    ----------
    Norm : float, optional
      Normaliztion
    tilt : float, optional
      Power-law tilt to continuum
    piv_wv : float, optional
      Pivot wave for tilt.  Best kept *without* units
    piv_wv2 : float, optional
      Pivot wave for a second tilt. Better be at wavelength < piv_wv
    igm : str, optional
      Adopt average IGM model? ['None']
    LL_flatten : bool, optional
      Set Telfer to a constant below the LL?

    Returns
    -------
    conti_dict : dict
      Useful for simple modeling.  Keep as a dict for JSON writing
    """
    conti_dict = dict(Norm=Norm, tilt=tilt, piv_wv=piv_wv, piv_wv2=piv_wv2,
                      tilt2=tilt2, igm=igm, fN_gamma=fN_gamma,
                      LL_flatten=LL_flatten)
    # Checks
    if piv_wv2 is None:
        conti_dict.pop('piv_wv2')
    else:
        if piv_wv is None:
            raise IOError("piv_wv required if piv_wv2 set")
        else:
            if piv_wv2 > piv_wv:
                raise ValueError("piv_wv2 < piv_wv required!")
    #
    return conti_dict


def quick_cont(spec,fitrange,initcontpts = 400):
    '''Generates an automatic continuum fit over a chunk of spectrum.

    Parameters
    ----------
    spec : XSpectrum1d
        Spectrum to be fitted
    fitrange : 2-element tuple of Quantity
        Tuple of wavelength range over which to snip spectrum and fit continuum


    Returns
    -------
    specchunk_cont : XSpectrum1d
        Snipped spectrum with continuum added as specchunk_cont.co
    '''

    wave = spec.wavelength.to(u.Angstrom).value
    flux = spec.flux
    sig = spec.sig
    wave1 = fitrange[0].to(u.Angstrom).value
    wave2 = fitrange[1].to(u.Angstrom).value
    #global fitidx,cont,wrange
    waveidx1=closest(wave,wave1)
    waveidx2=closest(wave,wave2)
    wrange=range(waveidx1,waveidx2+1)
    ### Generate automatic first-pass fit
    # Initialize points for fit
    initcontpts=400
    #ptspacing=(wave2-wave1)/initcontpts
    initx = np.linspace(wave1,wave2,initcontpts)
    #initx=wave1+np.arange(initcontpts)*ptspacing
    initidx=[closest(wave,ix) for ix in initx]
    fitidx=initidx
    #initidx=[]
    #for i in range(initcontpts):
    #	initidx.append(closest(wave,initx[i]))


    initwave = wave[fitidx]
    initflux = flux[initidx]
    initsig=sig[initidx]

    # Fit initial points and reject fit points far from continuum
    coeff,covmtx=fitcont(spec,fitidx,4)
    cont=L.legval(initwave,coeff)
    resarr=(cont-initflux)
    toremove=[]
    for i in range(len(fitidx)):
        if resarr[i]>(1.*initsig[i]):
            toremove.append(i)
    fitidx=np.delete(fitidx,toremove)
    contwave=wave[fitidx]
    contflux=flux[fitidx]
    contsig=sig[fitidx]

    # Fit again and reject again
    coeff,covmtx=fitcont(spec,fitidx,5)
    cont=L.legval(contwave,coeff)
    resarr=(cont-contflux)
    toremove=[]
    for i in range(len(cont)):
        if resarr[i]>(1.*contsig[i]):
            toremove.append(i)
    fitidx=np.delete(fitidx,toremove)
    contwave=wave[fitidx]
    contflux=flux[fitidx]
    contsig=sig[fitidx]

    # Now fit again and using all initial x-values as domain for continuum, reject intial pts from new cont.
    coeff,covmtx=fitcont(spec,fitidx,5)
    cont=L.legval(initwave,coeff)
    resarr=abs(cont-initflux)
    toremove=[]
    for i in range(len(initsig)):
        if resarr[i]>(1.*initsig[i]):
            toremove.append(i)
    fitidx=np.delete(initidx,toremove)
    contwave=wave[fitidx]
    contflux=flux[fitidx]
    contsig=sig[fitidx]

    # Final first-pass fit!
    coeff,covmtx=fitcont(spec,fitidx,6)
    cont=evalcont(wave[wrange],coeff)
    spec.wavelength = wave[wrange]*u.Angstrom
    spec.flux = flux[wrange]
    spec.sig = sig[wrange]
    return coeff,covmtx

def fitcont(spec,fitidx,maxorder=6):
    #global fitcoeff,fitcovmtx,cont,wrange
    wave = spec.wavelength.to(u.Angstrom).value
    flux = spec.flux
    sig = spec.sig
    wavepts=wave[fitidx]; fluxpts=flux[fitidx]; sigpts=sig[fitidx]
    ### Cycle through fits of several orders and decide the appropriate order by F-test
    coeffs=[]; covmtxs=[]; fits=[]; chisqs=[]; dfs=[]
    fprob=0.
    i=0; order=1
    while (fprob<=0.95) & (order<=maxorder):
        order=i+1
        coeff,covmtx=fitlegendre(wavepts,fluxpts,sigpts,order)
        coeffs.append(coeff) ; covmtxs.append(covmtx)
        fit=L.legval(wavepts,coeff,order)
        fits.append(fit)
        chisq,df=redchisq(fluxpts,fit,sigpts,order)
        chisqs.append(chisq),dfs.append(df)
        if i>0:
            fval=chisqs[i]/chisqs[i-1]
            fprob=stats.f.cdf(fval,dfs[i],dfs[i-1])
        i+=1
    ### Choose fit of order just prior to where F-test indicates no improvement
    fitcoeff=coeffs[i-2]; fitcovmtx=covmtxs[i-2]
    wrange=range(fitidx[0],fitidx[-1])
    #cont=L.legval(wave[wrange],fitcoeff)
    return fitcoeff,fitcovmtx

def fitlegendre(wavepts,fluxpts,sigpts,order):
	vander=L.legvander(wavepts,order)
	design=np.zeros(vander.shape)
	for i in range(len(wavepts)):
		design[i]=vander[i]/sigpts[i]
	U,s,v=svd(design,full_matrices=False,compute_uv=True)
	V=v.transpose()
	solvec=np.zeros(order+1)
	for i in range(order+1):
		solvec+=np.dot(U[:,i],fluxpts/sigpts)/s[i]*V[:,i]
	### Build covariance matrix
	covmtx=np.zeros([order+1,order+1])
	for j in range(order+1):
		for k in range(order+1):
			covmtx[j,k]=np.sum(V[:,j]*V[:,k]/s)
	return solvec,covmtx

def closest(array,value):
    '''Return the index of the closest element in array to value

    Parameters
    ----------
    array : 1d array
        Array from which to find element
    value : scalar
        Value to search for within array

    Returns
    -------
    idx : int
        Index of closest element to value in array

    '''
    idx = (np.abs(array-value)).argmin()
    return idx

def redchisq(fluxpts1,fitpts1,sigpts,order):
    diff=sum(fluxpts1-fitpts1)
    df=len(fluxpts1)-order-1
    Xsq=1./df*sum(diff**2/sigpts**2)
    return Xsq,df

def evalcont(wave,coeff):
	cont=L.legval(wave,coeff)
	return cont