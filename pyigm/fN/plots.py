""" Class for f(N) constraints
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import imp
import pdb

from pyigm.fN.fnmodel import FNModel
from pyigm.fN.constraints import FNConstraint
from pyigm.fN import tau_eff

# Path for pyigm
pyigm_path = imp.find_module('pyigm')[1]


def tst_fn_data(fN_model=None, model_two=None, data_list=None, outfil=None):
    """ Make a matplotlib plot like the final figure from P14
    See the Notebook for Bokeh plots
    Taken from xastropy without actually trying to run

    Parameters
    ----------
    fN_model
    model_two
    data_list
    outfil
    """
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'STIXGeneral-Regular' # Not for PDF
    mpl.rcParams['lines.linewidth'] = 2

    from matplotlib import pyplot as plt
    #from matplotlib.backends.backend_pdf import PdfPages

    # fN data
    all_fN_cs = FNConstraint.load_defaults()

    # Remove K12
    if data_list is None:
        fN_cs = [fN_c for fN_c in all_fN_cs if ((fN_c.ref != 'K02') & (fN_c.ref != 'PW09'))]
    else:
        fN_cs = [fN_c for fN_c in all_fN_cs if fN_c.ref in data_list]
    fN_dtype = [fc.fN_dtype for fc in fN_cs]

    fig = plt.figure(figsize=(8, 5))
    fig.clf()
    main = fig.add_axes( [0.1, 0.1, 0.8, 0.8] ) # xypos, xy-size

    # f(N) data
    main.set_ylabel(r'$\log f(N_{\rm HI})$')
    main.set_xlabel(r'$\log N_{\rm HI}$')
    main.set_ylim(-25., -9)

    for fN_c in fN_cs:
        if fN_c.fN_dtype == 'fN':
            # Length
            ip = range(fN_c.data['NPT'])
            val = np.where(fN_c.data['FN'][ip] > -90)[0]
            if len(val) > 0:
                ipv = np.array(ip)[val]
                xval = np.median(fN_c.data['BINS'][:,ipv],0)
                xerror = [ fN_c.data['BINS'][1,ipv]-xval, xval-fN_c.data['BINS'][0,ipv] ]
                yerror = [ fN_c.data['SIG_FN'][1,ipv], fN_c.data['SIG_FN'][0,ipv] ] # Inverted!
                main.errorbar(xval, fN_c.data['FN'][ipv], xerr=xerror, yerr=yerror, fmt='o',
                              label=fN_c.ref,capthick=2)
    main.legend(loc='lower left', numpoints=1)

    # Model?
    if fN_model is not None:
        xplt = 12.01 + 0.01*np.arange(1100)
        yplt = fN_model.evaluate(xplt, 2.4)
        main.plot(xplt,yplt,'-',color='black')
        print(xplt[0],yplt[0])
    if model_two is not None:
        xplt = 12.01 + 0.01*np.arange(1100)
        yplt = model_two.evaluate(xplt, 2.4)
        main.plot(xplt,yplt,'-',color='gray')

    # Extras
    #mpl.rcParams['lines.capthick'] = 2

    inset = fig.add_axes( [0.55, 0.6, 0.25, 0.25] ) # xypos, xy-size
    inset.set_ylabel('Value') # LHS
    inset.xaxis.set_major_locator(plt.FixedLocator(range(5)))
    inset.xaxis.set_major_formatter(plt.FixedFormatter(['',r'$\tau_{\rm eff}^{\rm Ly\alpha}$',
                                                        r'$\ell(X)_{\rm LLS}$',
                                                        r'$\lambda_{\rm mfp}^{912}$', '']))
    inset.set_ylim(0., 0.6)

    # tau_eff
    flg_teff = 1
    try:
        itau = fN_dtype.index('teff') # Passes back the first one
    except:
        flg_teff = 0

    if flg_teff:
        teff=float(fN_cs[itau].data['TEFF'])
        D_A = 1. - np.exp(-1. * teff)
        SIGDA_LIMIT = 0.1  # Allows for systemtics and b-value uncertainty
        sig_teff = np.max([fN_cs[itau].data['SIG_TEFF'], (SIGDA_LIMIT*teff)])
        # Plot
        inset.errorbar(1, teff, sig_teff, fmt='_', capthick=2)
        # Model
        if fN_model is not None:
            model_teff = tau_eff.lyman_ew(1215.6701*(1+fN_cs[itau].zeval), fN_cs[itau].zeval+0.1,
                                          fN_model, NHI_MIN=fN_cs[itau].data['NHI_MNX'][0],
                                          NHI_MAX=fN_cs[itau].data['NHI_MNX'][1])
            inset.plot(1, model_teff, 'ko')
            #xdb.set_trace()

    ## #######
    # LLS constraint
    flg_LLS = 1
    try:
        iLLS = fN_dtype.index('LLS') # Passes back the first one
    except:
        #raise ValueError('fN.data: Missing LLS type')
        flg_LLS = 0
    if flg_LLS:
        inset.errorbar(2, fN_cs[iLLS].data['LX'], yerr=fN_cs[iLLS].data['SIG_LX'],
                       fmt='_', capthick=2)
        # Model
        if fN_model != None:
            lX = fN_model.calculate_lox(fN_cs[iLLS].zeval,
                                        17.19+np.log10(fN_cs[iLLS].data['TAU_LIM']), 22.)
            inset.plot(2, lX, 'ko')

    ## #######
    # MFP constraint
    flg_MFP = 1
    try:
        iMFP = fN_dtype.index('MFP') # Passes back the first one
    except:
        #raise ValueError('fN.data: Missing MFP type')
        flg_MFP = 0

    if flg_MFP:
        inset2 = inset.twinx()
        inset2.errorbar(3, fN_cs[iMFP].data['MFP'], yerr=fN_cs[iMFP].data['SIG_MFP'],
                        fmt='_', capthick=2)
        inset2.set_xlim(0,4)
        inset2.set_ylim(0,350)
        inset2.set_ylabel('(Mpc)')

        # Model
        if fN_model is not None:
            #fN_model.zmnx = (0.1, 20.) # Reset for MFP calculation
            mfp = fN_model.mfp(fN_cs[iMFP].zeval)
            inset2.plot(3, mfp, 'ko')

    # Show
    if outfil != None:
        plt.savefig(outfil,bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':

    flg_fig = 0
    #flg_fig += 2**0   # Without model
    flg_fig += 2**1   # With model

    # Without model
    if (flg_fig % 2**1) >= 2**0:
        tst_fn_data()

    # With model
    if (flg_fig % 2**2) >= 2**1:
        fN_default = FNModel.default_model()
        tst_fn_data(fN_model=fN_default)
