""" Analysis code for IGMSurvey objects
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob
import json
import pdb

from astropy.table import Table

from linetools import utils as ltu



def calc_slgrid_atan(surveys, Agrid, Bgrid, Cgrid, C2grid):
    """ Calculate the sightline grid for a Atan l(z) fit
    Breaking this off for bootstrap speed-up

    Parameters
    ----------
    surveys : list of DLASurvey objects
    Agrid
    Bgrid
    Cgrid
    C2grid

    Returns
    -------
    slgrid : ndarray
      Sightline term in likelihood function
    """
    # Integrating over the sightlines
    slgrid = np.zeros_like(Agrid)
    # Int(atan[x-a]) = (a-x) atan(a-x) - 0.5 * ln(a**2 - 2ax + x**2 + 1)
    for isurvey in surveys:
        slines = isurvey.sightlines
        gds = slines['Z_START'] < slines['Z_END']
        zstart = slines['Z_START'][gds]
        zend = slines['Z_END'][gds]
        # Integrate constant term
        AAgrid = Agrid * np.sum(zend-zstart)
        slgrid += AAgrid
        # Integrate second term
        for iz in zend:
            CCgrid = (Cgrid-iz) * np.arctan(Cgrid-iz) - 0.5 * np.log(
                C2grid - 2*Cgrid*iz + iz**2 + 1)
            slgrid += Bgrid * CCgrid
            if np.min(CCgrid) < -0.1:
                pdb.set_trace()
        for iz in zstart:
            CCgrid = (Cgrid-iz) * np.arctan(Cgrid-iz) - 0.5 * np.log(
                C2grid - 2*Cgrid*iz + iz**2 + 1)
            slgrid -= Bgrid * CCgrid
    # Return
    return slgrid


def fit_atan_dla_lz(surveys, nstep=20, bootstrap=True,
                    nboot=10, nproc=2,
                    fit_out=None, boot_out=None,
                    verbose=True):
    """ Fit a A + B * atan(z-C)  l(z) model to AbsSys data
    Writes bootstrap analysis to hard-drive

    Code used in Prochaska & Neeleman 2017 for DLAs

    Parameters
    ----------
    surveys : list of IGMSurvey objects
      If None, a default list is loaded
    nstep : int, optional
      Steps in each dimension of the grid
    bootstrap : bool, optional
      Perform bootstrap analysis
    nboot : int, optional
      Number of bootstrap iterations
    nproc : int, optional
      Number of processors to use
    fit_out : str, optional
      Output filename for best fit (JSON)
    boot_out : str, optional
      Output filename for bootstrap analysis
    verbose : bool, optional

    Returns
    -------
    dfits : dict
      Best fit parameters
    boot_tbl : Table
      Returned if bootstrap=True
      else return None


    """
    # Name and date
    # Init
    if boot_out is None:
        boot_out = './lz_boot.fits.gz'
    if fit_out is None:
        fit_out = './lz_fit.json'
    # Synthesize
    all_z = np.concatenate([isurvey.zabs for isurvey in surveys])
    ndla = len(all_z)

    # Model :  l(z) = A + B * atan(C-z)
    Aparm = np.linspace(0.05, 0.5, num=nstep).astype(np.float32)
    Bparm = np.linspace(0.05, 0.5, num=nstep).astype(np.float32)
    Cparm = np.linspace(1., 6., num=nstep).astype(np.float32)

    # Generate grids (float32)
    Agrid, Bgrid, Cgrid = np.meshgrid(Aparm, Bparm, Cparm, copy=False)
    C2grid = Cgrid**2

    # Sightline grid
    if verbose:
        print("Sightline calculation...")
    slgrid = calc_slgrid_atan(surveys, Agrid, Bgrid, Cgrid, C2grid)

    if bootstrap:
        if verbose:
            print("Bootstrapping!")
        sv_fits = []
        rN = np.random.poisson(ndla, size=nboot)
        # Boot me
        z_list = []
        for kk,irN in enumerate(rN):
            # Draw nPoisson
            rval = (np.random.uniform(size=irN)*ndla).astype(int)
            # Draw from all_z
            draw_z = all_z[rval]
            z_list.append(draw_z)
        # Run
        if nproc == 1:
            for draw_z in z_list:
                if verbose:
                    print("Working on iteration: {:d} of {:d}".format(kk, nboot))
                dfits, _, _ = Ln_lz_atan(Agrid, Bgrid, Cgrid, slgrid, draw_z, write=False)
                # Save
                sv_fits.append(dfits.copy())
        else:
            import multiprocessing
            pool = multiprocessing.Pool(nproc) # initialize thread pool N threads
            inp_list = []
            for ii in range(nboot):
                inp_list.append(
                    dict(A=Agrid, B=Bgrid, C=Cgrid, sl=slgrid, z=z_list[ii]))
            if verbose:
                print("Mapping...")
            sv_fits = pool.map(map_Ln_atan, inp_list)
        # Write
        boot_tbl = Table()
        for key in ['A', 'B', 'C']:
            boot_tbl[key] = [ifits['lz']['atan'][key] for ifits in sv_fits]
        boot_tbl.write(boot_out, overwrite=True)
        if verbose:
            print("Wrote {:s}".format(boot_out))
    else:
        boot_tbl = None
    # Best
    dfits, _, _ = Ln_lz_atan(Agrid, Bgrid, Cgrid, slgrid, all_z, write=True)

    # Finish
    return dfits, boot_tbl


def Ln_lz_atan(Agrid, Bgrid, Cgrid, slgrid, all_z, write=True, verbose=True):
    """ Likelihood function for arctan model

    Parameters
    ----------
    Agrid
    Bgrid
    Cgrid
    slgrid
    all_z
    write

    Returns
    -------
    dfits : dict
      Contains best fit model
    dlagrid : ndarray
      for debugging
    lngrid : ndarray

    """
    # z0 estimate from 21cm surveys
    lz_z0 = dict(value=np.mean([0.026, 0.045]), sig=0.01)
    # Init
    dlagrid = np.zeros_like(Agrid)
    # Generate Likelihood for DLAs
    np.seterr(invalid='ignore')
    for z in all_z:
        dlagrid += np.log(Agrid + Bgrid * np.arctan(z-Cgrid))
    bad = np.isnan(dlagrid)
    dlagrid[bad] = -1e9

    # Likelihood
    lngrid = dlagrid - slgrid

    # z=0
    model_z0 = Agrid + Bgrid * np.arctan(0.-Cgrid)
    lnP = -1 * (model_z0-lz_z0['value'])**2 / 2 / (lz_z0['sig']**2)
    lngrid += lnP
    # Best
    indices = np.where(lngrid == np.max(lngrid))
    best = Agrid[indices][0], Bgrid[indices][0], Cgrid[indices][0]
    if verbose:
        print('Best fit: A={}, B={}, C={}'.format(best[0], best[1], best[2]))
    # Load
    dfits = {}
    # Write
    dfits['lz'] = {}
    dfits['lz']['atan'] = dict(A=Agrid[indices][0], B=Bgrid[indices][0], C=Cgrid[indices][0],
                            form='A + B*atan(z-C)')
    # Return
    return dfits, dlagrid, lngrid


def map_Ln_atan(map_dict):
    """ For multiprocessing the bootstrap

    Parameters
    ----------
    map_dict

    Returns
    -------

    """
    dfits, _, _ = Ln_lz_atan(map_dict['A'], map_dict['B'], map_dict['C'],
                             map_dict['sl'], map_dict['z'], write=False,
                             verbose=False)
    return dfits



def fit_fN_dblpow(NHI, a3_mnx, a4_mnx, Nd_mnx, nstep=100,
                  Nmin=10**(20.3), Nmax=1e99, verbose=True):
    """  Fit a double power-law to an input NHI distribution
    Only does the shape

    Done in float32 to preserve memory

    Code from Prochaska & Neeleman (2017)  [and also PHW05]

    Parameters
    ----------
    NHI : ndarray
      log10 NHI values
    a3_mnx : tuple
      min/max of lower NHI power-law
    a4_mnx : tuple
      min/max of upper NHI power-law
    Nd_mnx : tuple
      min/max of break column in log10
    nstep : int, optional
    Nmin : float, optional
      Minimum NHI in the analysis [usually DLA criterion]
    Nmax : float, optional
      Maximum NHI in the analysis

    Returns
    -------
    dfits : dict
      Contains the fit
    best : tuple
      Best fit values in grid for Nd, a3, a4
    Ndgrid
    a3grid
    a4grid
    lik

    """
    # Generate 1D arrays
    a3stp = np.linspace(a3_mnx[0], a3_mnx[1], nstep).astype(np.float32)
    a4stp = np.linspace(a4_mnx[0], a4_mnx[1], nstep).astype(np.float32)
    Ndstp = np.linspace(Nd_mnx[0], Nd_mnx[1], nstep).astype(np.float32)

    # Generate grids (float32)
    a3grid, a4grid, Ndgrid = np.meshgrid(a3stp, a4stp, Ndstp, copy=False)

    # Linear
    Ns10 = 10.**Ndgrid

    # Calculate denominator
    denom = Ns10 * ((1. - (Nmin / Ns10)**(a3grid + 1.)) / (1. + a3grid) + (
        (Nmax / Ns10)**(a4grid + 1) - 1.) / (a4grid + 1.))

    num = np.zeros_like(Ns10)
    # Numerator
    # Loop on DLAs
    for iNHI10 in 10.**NHI:
        # Upper end
        high = iNHI10 > Ns10
        if np.sum(high) > 0:
            num[high] += a4grid[high] * np.log(iNHI10 / Ns10[high])
        # Low end
        if np.sum(~high) > 0:
            num[~high] += a3grid[~high] * np.log(iNHI10 / Ns10[~high])

    # Liklihood (Beware of Signs!)
    lik = num - NHI.size * np.log(denom)

    mxL = np.max(lik)
    indices = np.where(lik == mxL)
    best = Ndgrid[indices][0], a3grid[indices][0], a4grid[indices][0]
    if verbose:
        print('Best fit: Nd={}, a3={}, a4={}'.format(best[0], best[1], best[2]))

    # Load
    dfits = {}
    # Write
    dfits['fN'] = {}
    dfits['fN']['dpow'] = dict(Nd=Ndgrid[indices][0], a3=a3grid[indices][0], a4=a4grid[indices][0],
                            form='(N/Nd)**aa with aa=a3 if N<Nd else aa=a4')

    # KS Test
    ks_test = False
    if ks_test:
        ns10 = 10**best[0]
        dblpow_k = 1. / (ns10 * (1. - (Nmin / Ns10)**(best[1] + 1)) / (1. + best[1]) + (
            (Nmax / Ns10)**(best[2] + 1) - 1.) / (best[2] + 1))
        dblpow_b1 = best[1]
        dblpow_b2 = best[2]
        dblpow_nd = ns10
        dblpow_nmin = Nmin
        noise = 0.02
        dNHI = 10**(NHI + noise * np.random.uniform(size=NHI.size))
        #ksone, darr, 'x_maxdblpow_kscumf', d, ksprob

    return dfits, best, Ndgrid, a3grid, a4grid, lik
