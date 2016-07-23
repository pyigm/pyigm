""" This module will fit f(N) data
using a Markov chain Monte Carlo (MCMC) approach.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, imp
import numpy as np
import pdb
import pymc
#import MCMC_errors

from pyigm.fN import plots
from pyigm.fN.fnmodel import FNModel
from pyigm.fN.constraints import FNConstraint
from pyigm.fN import tau_eff

from time import gmtime, strftime

pyigm_path = imp.find_module('pyigm')[1]


def set_fn_model(flg=0):
    """ Load up f(N) model

    Parameters
    ----------
    flg : int, optional
      * 0 = Default model
      * 1 = Gamma

    Returns
    -------
    sfN_model : FNModel
    """
    if flg==0: # I may choose to pickle a few of these
        sfN_model = FNModel.default_model(use_mcmc=True) # Hermite Spline
    elif flg==1:
        sfN_model = FNModel('Gamma')
    else: 
        raise ValueError('mcmc.set_model: Not ready for this type of fN model {:d}'.format(flg))
    #
    return sfN_model


def set_fn_data(sources=None, extra_fNc=[]):
    """ Load up f(N) data

    Parameters
    ----------
    sources : list, optional
      References for constraints
    extra_fNc : list, optional

    Returns
    -------
    fN_data :: List of fN_Constraint Classes
    """
    if sources is None:
        sources = ['OPB07', 'OPW12', 'OPW13', 'K05', 'K13R13', 'N12']

    fn_file = pyigm_path+'/data/fN/fn_constraints_z2.5_vanilla.fits'
    k13r13_file = pyigm_path+'/data/fN/fn_constraints_K13R13_vanilla.fits'
    n12_file = pyigm_path+'/data/fN/fn_constraints_N12_vanilla.fits'
    all_fN_cs = FNConstraint.from_fitsfile([fn_file,k13r13_file,n12_file])

    # Add on, e.g. user-supplied
    if len(extra_fNc) > 0:
        raise IOError("NOT READY FOR THIS YET")
        #for src in extra_fNc:
        #	all_fN_cs.append(xifd.fN_data_from_ascii_file(os.path.abspath(src)))

    # Include good data sources
    fN_cs = []
    for fN_c in all_fN_cs:
        # In list?
        if fN_c.ref in sources:
            print('Using {:s} as a constraint'.format(fN_c.ref))
            # Append
            fN_cs.append(fN_c)
            # Pop
            idx = sources.index(fN_c.ref)
            sources.pop(idx)
    
    # Check that all the desired sources were used
    if len(sources) > 0:
        pdb.set_trace()

    return fN_cs

##########################################
#   Prepare the variables and their limits
##########################################
def set_pymc_var(fN_model,lim=2.):
    """ Generate pymc variables

    Parameters
    ----------
    fN_model : FNModel
    lim: float, optional
      Range limits for Uniform stochastic value

    Returns
    -------
    Array of pymc Stochastic variables
    """
    iparm=np.array([])

    # Deal with model Type
    if fN_model.fN_mtype == 'Hspline': 
        # Loop on parameters to create an array of pymc Stochatsic variable objects
        nparm = len(fN_model.param['sply'])
        rand = 0.2*(np.random.rand(nparm)-0.5)
        #pdb.set_trace()
        for ii in range(nparm):
            nm = str('p')+str(ii)
            doc = str('SplinePointNHI_')+str(fN_model.pivots[ii])
            #iparm = np.append(iparm, pymc.Uniform(nm, lower=fN_model.param[ii]-lim,
            #                                    upper=fN_model.param[ii]+lim, doc=doc))
            iparm = np.append(iparm, pymc.Normal(nm, mu=fN_model.param['sply'][ii]*(1+rand[ii]), tau=1./0.025, doc=doc))
    elif fN_model.fN_mtype == 'Gamma':  # Inoue+14
        raise ValueError("NOT UPDATED FOR THIS")
        rand = 0.2*(np.random.rand(4)-0.5) 
        fN_model.param[2][0] = 10.
        #fN_model.param[2][0] = fN_model.param[2][0]*(1.+rand[0])
        fN_model.param[2][1] = fN_model.param[2][1]*(1.+rand[1])
        fN_model.param[3][0] = fN_model.param[3][0]*(1.+rand[2])
        fN_model.param[3][1] = fN_model.param[3][1]*(1.+rand[3])
        # LAF: Only vary A and beta as a first go
        iparm = np.append(iparm, pymc.Normal(str('p0'), mu=fN_model.param[2][0], tau=1./50., doc=str('ALAF')))
        iparm = np.append(iparm, pymc.Normal(str('p1'), mu=fN_model.param[2][1], tau=1./0.05, doc=str('bLAF')))
        # DLA: Only vary A and beta as a first go
        iparm = np.append(iparm, pymc.Normal(str('p2'), mu=fN_model.param[3][0], tau=1./0.25, doc=str('ADLA')))
        iparm = np.append(iparm, pymc.Normal(str('p3'), mu=fN_model.param[3][1], tau=1./0.05, doc=str('bDLA')))
    else:
        raise ValueError('mcmc: Not ready for this type of fN model {:s}'.format(fN_model.fN_mtype))
    # Return
    return iparm


def run(fN_cs, fN_model, parm, debug=False):
    """ Run the MCMC

    Parameters
    ----------
    fN_cs
    fN_model
    parm
    debug : bool, optional

    Returns
    -------

    """
    #
    pymc_list = [parm]

    # Parse data and combine as warranted
    all_NHI = []
    all_fN = []
    all_sigfN = []
    all_z = []
    flg_teff = 0
    flg_LLS = 0
    for fN_c in fN_cs: 
        # Standard f(N)
        if fN_c.fN_dtype == 'fN':
            ip = range(fN_c.data['NPT'])
            val = np.where(fN_c.data['FN'][ip] > -90)[0] # Deal with limits later
            ipv = np.array(ip)[val]
            # Append the NHI
            NHI = np.median(fN_c.data['BINS'][:,ipv],0)
            all_NHI += list(NHI)
            # Append the f(N)
            all_fN += list(fN_c.data['FN'][ipv])
            # Append the Error
            fNerror = np.median(fN_c.data['SIG_FN'][:,ipv],0)
            all_sigfN += list(fNerror)
            # Append zeval
            for ii in range(len(ipv)):
                all_z.append(fN_c.zeval)
        elif fN_c.fN_dtype == 'teff': # teff_Lya
            if flg_teff:
                raise ValueError('Only one teff allowed for now!')
            else:
                flg_teff = 1
            teff=float(fN_c.data['TEFF'])
            D_A = 1. - np.exp(-1. * teff)
            SIGDA_LIMIT = 0.1  # Allows for systemtics and b-value uncertainty
            sig_teff = np.max([fN_c.data['SIG_TEFF'], (SIGDA_LIMIT*teff)])
            teff_zeval = float(fN_c.data['Z_TEFF'])

            # Save input for later usage
            teff_input = (teff_zeval, fN_c.data['NHI_MNX'][0], fN_c.data['NHI_MNX'][1])
        elif fN_c.fN_dtype == 'l(X)': # teff_Lya
            if flg_LLS:
                raise ValueError('Only one teff allowed for now!')
            else:
                flg_LLS = 1
            LLS_lx = fN_c.data['LX']
            LLS_siglx = fN_c.data['SIG_LX']
            LLS_input = (fN_c.zeval, fN_c.data['TAU_LIM'])
            
    # 
    fN_input = (np.array(all_NHI), np.array(all_z))
    #flg_teff = 0

    """
    Generate the Models
    """
    @pymc.deterministic(plot=False)
    def pymc_fn_model(parm=parm):
        # Define f(N) model for PyMC
        # Set parameters
        aparm = np.array([parm[i] for i in range(parm.size)])
        fN_model.update_parameters(aparm)
        #
        log_fNX = fN_model.evaluate( fN_input, 0. )
        #
        return log_fNX
    pymc_list.append(pymc_fn_model)

    if flg_teff:
        # Define teff model for PyMC
        @pymc.deterministic(plot=False)
        def pymc_teff_model(parm=parm):
            # Set parameters
            aparm = np.array([parm[i] for i in range(parm.size)])
            fN_model.update_parameters(aparm)
            # Calculate teff
            model_teff = tau_eff.lyman_ew(1215.6701*(1+teff_input[0]), teff_input[0]+0.1,
                                               fN_model, NHI_MIN=teff_input[1], NHI_MAX=teff_input[2])
            return model_teff
        pymc_list.append(pymc_teff_model)

    if flg_LLS:
        # Define l(X)_LLS model for PyMC
        @pymc.deterministic(plot=False)
        def pymc_lls_model(parm=parm): 
            # Set parameters 
            aparm = np.array([parm[i] for i in range(parm.size)])
            fN_model.update_parameters(aparm)
            # Calculate l(X)
            lX = fN_model.calculate_lox(LLS_input[0], 17.19+np.log10(LLS_input[1]), 22.)
            return lX
        pymc_list.append(pymc_lls_model)

    """
    Generate the Data
    """

    # Define f(N) data for PyMC
    fNvalue=np.array(all_fN)
    #xdb.set_trace()
    pymc_fN_data = pymc.Normal(str('fNdata'), mu=pymc_fn_model, tau=1.0/np.array(all_sigfN)**2,
                               value=fNvalue, observed=True)
    pymc_list.append(pymc_fN_data)

    # Define teff data for PyMC
    if flg_teff:
        pymc_teff_data = pymc.Normal(str('teffdata'), mu=pymc_teff_model, tau=1.0/np.array(sig_teff)**2,
                                value=teff, observed=True)
        pymc_list.append(pymc_teff_data)

    # Define l(X)_LLS model for PyMC
    if flg_LLS:
        pymc_lls_data = pymc.Normal(str('LLSdata'), mu=pymc_lls_model, tau=1.0/np.array(LLS_siglx)**2,
                                value=LLS_lx, observed=True)
        pymc_list.append(pymc_lls_data)


    """
    RUN THE MCMC
    """

    MC = pymc.MCMC(pymc_list)#,verbose=2)
    # Force step method to be Metropolis!
    for ss in MC.stochastics-MC.observed_stochastics:
        MC.use_step_method(pymc.Metropolis, ss, proposal_sd=0.025, proposal_distribution='Normal')
    #xdb.set_trace()

    # Run a total of 40000 samples, but ignore the first 10000.
    # Verbose just prints some details to screen.
    #xdb.set_trace()
    #MC.sample(20000, 3000, verbose=2, tune_interval=500)
    #MC.sample(5000, 500, verbose=2, tune_interval=200)
    #MC.sample(20000, 5000, verbose=2, tune_interval=500)
    MC.sample(2000, 400, verbose=2, tune_interval=200)
    #MC.isample(10000, 1000, verbose=2)
    
    if debug:
        pdb.set_trace()
        #xifd.tst_fn_data(fN_model=fN_model)
        #xdb.xhist(MC.trace(str('p0'))[:])
        #xdb.set_trace()

    
    # Draw a contour plot with 1 & 2 sigma errors
    #MCMC_errors.draw_contours(MC, 'p0', 'p1')

    # Save the individual distributions to a file to check convergence
    #pymc.Matplot.plot(MC)
    #xdb.set_trace()
    return MC


def geterrors(array):
    """
    Parameters
    ----------
    array

    Returns
    -------

    """
    arrsort = np.sort(array)
    arrsize = np.size(array)
    value = arrsort[int(round(0.5*arrsize))]
    err1 = np.array([value-arrsort[int(round(0.15866*arrsize))],arrsort[int(round(0.84134*arrsize))]-value])
    err2 = np.array([value-arrsort[int(round(0.02275*arrsize))],arrsort[int(round(0.97725*arrsize))]-value])
    return value, err1, err2

def print_errors(MC):
    """
    Parameters
    ----------
    MC

    Returns
    -------

    """
    keys = MC.stats().keys()
    keys_size = len(keys)

    all_pval = []
    for ival in range(keys_size):
        teststr = str('p'+str(ival))
        #xdb.set_trace()
        if not teststr in keys:
            continue
        try: pval, perr1, perr2 = geterrors(MC.trace(teststr)[:])
        except: continue
        print('{:s} {:5.4f} +{:5.4f}-{:5.4f} (1sig) +{:5.4f}-{:5.4f} (2sig)'.format(
            teststr, pval, perr1[0], perr1[1], perr2[0], perr2[1]))
        all_pval.append(pval)
        #ival += 1
    return all_pval
    
def save_figures(MC, fN_model, email=None):
    """
    Parameters
    ----------
    MC
    fN_model
    email

    Returns
    -------

    """
    #xdb.set_trace()
    #######################################
    #   SAVE THE RESULTS
    #######################################
    #Creates new directory for output
    t = strftime("%Y-%m-%d%H%M%S", gmtime())
    #newpath = 'C:/Xastropy Output Files/' + email
    #if not os.path.exists(newpath): os.makedirs(newpath)
    
    #creates ascii file with the best values and their errors 
    #and saves to correct directory
    if email is not None:
        asciifilename = email + t
    else:
        asciifilename = 'mcmc_' + t
    completeAsciiName = os.path.join(asciifilename+".ascii")
    f = open(completeAsciiName, 'w+')
    best_pval = print_errors(MC)
    #db.set_trace()
    f.write(str(best_pval))
    f.close()
    
    #creates PNG file with test plot from data
    #png1filename= email + t + 'png1'
    #completepng1name= os.path.join(newpath, png1filename + ".png")
    #g = open(completepng1name, 'w+')
    #g.write(xifd.tst_fn_data(fN_model=fN_model))
    #g.close()
    
    #creates PNG file with bottom plot (individual distributions?)
    if email is not None:
        png2filename = email + t
    else:
        png2filename = 'mcmc_' + t
    completepng2name= os.path.join(png2filename + ".png")
    pymc.Matplot.plot(MC)
    pymc.Matplot.savefig(completepng2name)


def mcmc_main(email, datasources, extrasources, flg_model=0, flg_plot=0):
    """
    Parameters
    ----------
    email
    datasources
    extrasources
    flg_model : int
     Flag controlling the f(N) model fitted
       0: JXP spline
       1: Inoue+14 functional form
    flg_plot

    Returns
    -------

    """
    import argparse

    # PARSE 
    #parser = argparse.ArgumentParser(description='MCMC for f(N)')
    #parser.add_argument("model", type=int, help="Model flag (0=JXP, 1=Inoue+14)")
    #parser.add_argument("--plot", help="Plot the final model", action="store_true")
    #parser.add_argument("-id_dir", type=str,
    #                    help="Directory for ID files (ID_LINES is default)")
    
    #args = parser.parse_args()

    #flg_model = args.model
    #if args.plot:
        #flg_plot = 1

    # ##########################
    # Set Data
    fN_data = set_fn_data(datasources, extrasources)
    
    # Set f(N) functional form 
    fN_model = set_fn_model(flg=flg_model)
    
    # Set variables
    parm = set_pymc_var(fN_model)
    
    # Check plot
    if flg_plot:
        plots.tst_fn_data(fN_model=fN_model)

    # Run
    MC = run(fN_data, fN_model, parm, email)

    # Save files
    save_figures(MC, email, fN_model)
    
    # Plot?
    #if flg_plot:
        #xifd.tst_fn_data(fN_model=fN_model)


def chain_stats(chain_file, burn_frac=0.3, cl=0.683):
    """ Turn an MCMC chain into stats
    Port of x_mcmc_chain_stats from XIDL

    Parameters:
      chain_file: string
          Name of MCMC file
      burn_frac: float (0.3)
          Fraction of chain to burn
      cl: float (0.683)
          Confidence interval

    Returns:
      A dictionary with the key outputs

    JXP 07 Nov 2014
    """
    from astropy.io import fits

    # Read
    hdu = fits.open(chain_file)
    chain = hdu[0].data
    like = hdu[1].data

    # Param
    nparm = chain.shape[-1]
    if len(chain.shape) < 3: nchain = 1
    else: nchain = chain.shape[0]

    # Output
    outp = {}

    # Burn
    burn = int( np.round(chain.shape[1] * burn_frac ) )
    chain = chain[:,burn:,:]    # Burn
    like = like[:,burn:]    # Burn
    sz = chain.shape

    # Reshape
    if len(sz) == 3:
        chain = chain.reshape(sz[0]*sz[1],sz[2])
        chain = chain.transpose()
        sz = chain.shape
        like = like.flatten()

    # Maximize
    imx = np.argmax(like)
    outp['best_p'] = chain[:,imx]
    outp['sig'] = np.zeros((sz[0],2))

    # Confidence limits (68%)
    cnt = 0
    for row in chain: #qq=0L,nparm-1 do begin
        #; Sort
        srt = np.sort(row)
        # Simple CDF
        lowv = srt[np.round(sz[1]*(1-cl)/2.)]
        outp['sig'][cnt,0] = np.fabs(outp['best_p'][cnt] - lowv)
        #
        hiv = srt[np.round(sz[1]*(1.-((1-cl)/2.)))]
        outp['sig'][cnt,1] = np.fabs(hiv-outp['best_p'][cnt])
        cnt += 1

    return outp


