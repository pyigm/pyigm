""" Class for f(N) modeling
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
from scipy import interpolate as scii
import pdb

from astropy.io import fits
from astropy.utils.misc import isiterable
from astropy import units as u
from astropy import constants as const
from astropy import cosmology

from pyigm import utils as pyigmu

# Path for pyigm
pyigm_path = imp.find_module('pyigm')[1]


class FNModel(object):
    """A Class for f(N,X) models

    Attributes
    ----------
    mtype : str
      Model type for the fN
      'Hspline' -- Hermite monotonic spline
      'Gamma' -- Gamma function [Following Inoue+14]
    zmnx : tuple
      Redshift range where this model applies (zmin,zmax)
    pivots : array
          log NHI values for the pivots
    param : dict
      Parameter dict
        * Spline points for the Hspline
        * Amplitudes, power-laws, etc. for Gamma
    zpivot : float, optional
          Pivot for redshift evolution (2.4)
    gamma : float, optional
          Power law for dN/dX, not dN/dz (1.5)
    """
    @classmethod
    def default_model(cls, use_mcmc=False, cosmo=None):
        """ Pass back a default fN_model from Prochaska+14

        Tested against XIDL code by JXP on 09 Nov 2014

        Parameters
        ----------
        recalc : boolean, optional (False)
          Recalculate the default model
        use_mcmc : boolean, optional (False)
          Use the MCMC chain to generate the model
        write : boolean, optional (False)
          Write out the model
        cosmo : astropy.cosmology, optional
          Cosmology for some calculations
        """
        if use_mcmc:
            from pyigm.fN import mcmc
            #raise ValueError("DEPRECATED")
            # MCMC Analysis (might put these on a website)
            chain_file = (os.environ.get('DROPBOX_DIR')+
                        'IGM/fN/MCMC/mcmc_spline_k13r13o13n12_8.fits.gz')
            outp = mcmc.chain_stats(chain_file)

            # Build a model
            NHI_pivots = [12., 15., 17.0, 18.0, 20.0, 21., 21.5, 22.]
            fN_model = cls('Hspline', zmnx=(0.5,3.0),
                            pivots=NHI_pivots,
                           param=dict(sply=np.array(outp['best_p'].flatten()))) # fN_data['FN']).flatten()))
                           #param=outp['best_p'])
        else:
            # Input the f(N) at z=2.4 from Prochaska+14
            fN_file = (pyigm_path+'/data/fN/fN_spline_z24.fits.gz')
            print('Using P14 spline values to generate a default model')
            print('Loading: {:s}'.format(fN_file))
            hdu = fits.open(fN_file)
            fN_data = hdu[1].data
            # Instantiate
            fN_model = cls('Hspline', zmnx=(0.5, 3.0), cosmo=cosmo,
                            pivots=np.array(fN_data['LGN']).flatten(),
                            param=dict(sply=np.array(fN_data['FN']).flatten()))
        # Return
        return fN_model

    def __init__(self, mtype, zmnx=(0., 0.), pivots=None, param=None,
                 zpivot=2.4, gamma=1.5, cosmo=None):
        """
        Parameters
        ----------
        mtype : str
          Type of model.  Gamma,  Hspline, PowerLaw
        zmnx : tuple, optional
          Redshift range for evaluation
        pivots : ndarray, optional
          NHI values for the Hspline model
        param : dict, optional
          Parameters of the model.  Format depends on model
        zpivot : float, optional
          Pivot redshift for (1+z)^gamma redshift evolution
        gamma : float, optional
          Power-law exponent for (1+z)^gamma redshift evolution
        cosmo : astropy.cosmology, optional
          Cosmology for some calculations

        Returns
        -------
        FNModel

        """
        if mtype not in ['Gamma', 'Hspline', 'PowerLaw']:
            raise IOError("Invalid mtype")
        self.mtype = mtype  # Should probably check the choice

        # Pivots
        if pivots is None:
            self.pivots = np.zeros(2)
        else:
            self.pivots = pivots
        self.npivot = len(self.pivots)

        # Param
        if param is None:
            self.param = {}
        else:
            self.param = param

        # Model
        if mtype == 'Gamma':  # I+14
            zmnx = (0., 10.)
            self.param = dict(common=dict(Nl=12., Nu=23., Nc=21., bval=28.),  # Common
                              Bi=[1.75828e8, 9.62288e-4],             # Bi values
                              LAF=dict(Aval=500, beta=1.7, zcuts=[1.2, 4.7], gamma=[0.2, 2.7, 4.5]), # LAF
                              DLA=dict(Aval=1.1, beta=0.9, zcuts=[2.0], gamma=[1.0, 2.0]), # DLA
                              )
        elif mtype == 'Hspline':
            if len(self.param) == 0:
                raise IOError("Need to specify sply parameters!")
            self.model = scii.PchipInterpolator(self.pivots, self.param['sply'],
                                                extrapolate=True)  # scipy 0.16
        elif mtype == 'PowerLaw':
            if len(self.param) == 0:
                raise IOError("Need to specify PowerLaw parameters!")

        #
        self.fN_mtype = mtype
        self.zmnx = zmnx

        # Redshift evolution (needs updating)
        self.zpivot = zpivot
        self.gamma = gamma
        self.cosmo = cosmo

    def update_parameters(self, parm):
        """ Update parameters (mainly used in the MCMC)

        Updates other things as needed

        Parameters
        ----------
        parm : ndarray
          Parameters for the f(N) model to update to
        """
        if self.mtype == 'Hspline':
            self.param['sply'] = np.array(parm)
            # Need to update the model too
            self.model = scii.PchipInterpolator(self.pivots, self.param['sply'])
        elif self.mtype == 'Gamma':
            if len(parm) == 4: # A,beta for LAF and A,beta for DLA
                self.param['LAF']['Aval'] = parm[0]
                self.param['LAF']['beta'] = parm[1]
                self.param['DLA']['Aval'] = parm[2]
                self.param['DLA']['beta'] = parm[3]
            else:
                raise ValueError('fN/model: Not ready for {:d} parameters'.format(
                    len(parm)))

    def calculate_lox(self, z, NHI_min, NHI_max=None, neval=10000,
                      cosmo=None, cumul=False):
        """ Calculate l(X) over an N_HI interval

        Parameters
        ----------
        z : float
          Redshift for evaluation
        NHI_min : float
          minimum log NHI value
        NHI_max : float, optional
          maximum log NHI value for evaluation (Infinity)
        neval : int, optional
          Discretization parameter (10000)
        cumul : bool, optional
          Return a cumulative array? (False)
        cosmo : astropy.cosmology, optional
          Cosmological model to adopt (as needed)

        Returns
        -------
        lX : float
          l(X) value
        """
        # Initial
        if NHI_max == None:
            NHI_max = 23.
            infinity=True
        else:
            infinity=False

        try:
            nz = len(z)
        except:
            nz=1
            z = np.array([z])

        # Brute force (should be good to ~0.5%)
        lgNHI = np.linspace(NHI_min,NHI_max,neval)  #NHI_min + (NHI_max-NHI_min)*np.arange(neval)/(neval-1.)
        dlgN = lgNHI[1]-lgNHI[0]

        # Evaluate f(N,X)
        lgfNX = self.evaluate(lgNHI, z, cosmo=cosmo)

        # Sum
        lX = np.zeros(nz)
        for ii in range(nz): 
            lX[ii] = np.sum(10.**(lgfNX[:, ii]+lgNHI)) * dlgN * np.log(10.)
        if cumul == True:
            if nz > 1:  #; Have not modified this yet
                raise ValueError('fN.model: Not ready for this model type %s' % self.mtype)
            cum_sum = np.cumsum(10.**(lgfNX[:, ii]+lgNHI)) * dlgN * np.log(10.)

        # Infinity?
        if infinity is True:
            # This is risky...
            # Best to cut it off
            neval2 = 1000L
            lgNHI2 = NHI_max + (99.-NHI_max)*np.arange(neval2)/(neval2-1.)
            dlgN = lgNHI2[1] - lgNHI2[0]
            lgfNX = np.zeros((neval2, nz))
            lX2 = np.zeros(nz)
            for ii in range(nz):
                lgfNX[:, ii] = self.evaluate(lgNHI2, z[ii], cosmo=cosmo).flatten()
                lX2[ii] = np.sum(10.**(lgfNX[:, ii]+lgNHI2)) * dlgN * np.log(10.)
            lX = lX + lX2

        # Return
        if nz == 1:
            lX = lX[0]
        if cumul == True:
            return lX, cum_sum, lgNHI
        else:
            return lX

    def calculate_rhoHI(self, z, NHI_mnx, neval=10000, cumul=False, cosmo=None):
        """ Calculate rho_HI over an N_HI interval

        Parameters
        ----------
        z : float
          Redshift for evaluation
        NHI_mnx : tuple
          of floats minimum/maximum log NHI values
        neval : int, optional
          Discretization parameter (10000)
        cumul : bool, optional
          Return lgNHI, cumul_rhoHI too
        cosmo : astropy.cosmology, optional
          Needed for H0 only

        Returns
        -------
        rho_HI: float
          rho_HI in units of Msun per comoving Mpc**3
        if cumul is True
          lgNHI : array of log NHI values
          cumul : array of cumulative rho_HI
        """
        # Cosmology
        if cosmo is None:
            cosmo = cosmology.core.FlatLambdaCDM(70., 0.3)

        # Initial
        NHI_min=NHI_mnx[0]
        NHI_max=NHI_mnx[1]
        try:
            nz = len(z)
        except:
            nz=1
            z = np.array([z])

        # Brute force (should be good to ~0.5%)
        lgNHI = NHI_min + (NHI_max-NHI_min)*np.arange(neval)/(neval-1.)
        dlgN = lgNHI[1]-lgNHI[0]

        # Evaluate f(N,X)
        lgfNX = self.evaluate(lgNHI, z, cosmo=cosmo)

        # Sum
        rho_HI = np.zeros(nz)
        for ii in range(nz): 
            rho_HI[ii] = np.sum(10.**(lgfNX[:, ii]+2*lgNHI)) * dlgN * np.log(10.)
        if cumul is True:
            if nz > 1:  #; Have not modified this yet
                raise ValueError('fN.model: Not ready for this model type %s' % self.mtype)
            cum_sum = np.cumsum(10.**(lgfNX[:, ii]+2*lgNHI)) * dlgN * np.log(10.)

        # Constants
        rho_HI = rho_HI * (const.m_p.cgs * cosmo.H0 /
            const.c.cgs / (u.cm**2)).to(u.Msun/u.Mpc**3)

        # Return
        if nz == 1:
            rho_HI = rho_HI[0]
        if cumul is True:
            return rho_HI, cum_sum, lgNHI
        else:
            return rho_HI

    # Evaluate
    def evaluate(self, NHI, z, vel_array=None, cosmo=None):
        """ Evaluate the f(N,X) model at a set of NHI values

        Parameters
        ----------
        NHI : array
          log NHI values
        z : float or array
          Redshift for evaluation
        vel_array : ndarray, optional
          Velocities relative to z
        cosmo : astropy.cosmology, optional


        Returns
        -------
        log_fNX : float, array, or 2D array
          f(NHI,X)[z] values
          Float if given one NHI,z value each. Otherwise 2D array
          If 2D, it is [NHI,z] on the axes

        """
        # Tuple?
        if isinstance(NHI, tuple):  # All values packed into NHI parameter
            z = NHI[1]
            NHI = NHI[0]
            flg_1D = 1
        else:  # NHI and z separate
            flg_1D = 0

        # NHI
        if isiterable(NHI):
            NHI = np.array(NHI)  # Insist on array
        else:
            NHI = np.array([NHI])
        lenNHI = len(NHI)

        # Redshift 
        if vel_array is not None:
            z_val = z + (1+z) * vel_array/(const.c.to('km/s').value)
        else:
            z_val = z
        if isiterable(z_val):
            z_val = np.array(z_val)
        else:
            z_val = np.array([z_val])
        lenz = len(z_val)

        # Check on zmnx
        bad = np.where((z_val < self.zmnx[0]) | (z_val > self.zmnx[1]))[0]
        if len(bad) > 0:
            raise ValueError(
                'fN.model.eval: z={:g} not within self.zmnx={:g},{:g}'.format(
                    z_val[bad[0]], *(self.zmnx)))

        if self.mtype == 'Hspline':
            # Evaluate without z dependence
            log_fNX = self.model.__call__(NHI)

            # Evaluate
            if (not isiterable(z_val)) | (flg_1D == 1):  # scalar or 1D array wanted
                log_fNX += self.gamma * np.log10((1+z_val)/(1+self.zpivot))
            else:
                # Matrix algebra to speed things up
                lgNHI_grid = np.outer(log_fNX, np.ones(len(z_val)))
                lenfX = len(log_fNX)
                #
                z_grid1 = 10**(np.outer(np.ones(lenfX)*self.gamma,
                                        np.log10(1+z_val)))  #; (1+z)^gamma
                z_grid2 = np.outer( np.ones(lenfX)*((1./(1+self.zpivot))**self.gamma), 
                            np.ones(len(z_val)))
                log_fNX = lgNHI_grid + np.log10(z_grid1*z_grid2) 

        # Gamma function (e.g. Inoue+14)
        elif self.mtype == 'Gamma':
            # Setup the parameters
            Nl, Nu, Nc, bval = [self.param['common'][key]
                                for key in ['Nl', 'Nu', 'Nc', 'bval']]
            # gNHI
            Bi = self.param['Bi']
            ncomp = len(Bi)
            log_gN = np.zeros((lenNHI, ncomp))
            beta = [self.param[itype]['beta'] for itype in ['LAF', 'DLA']]
            for kk in range(ncomp):
                log_gN[:, kk] += (np.log10(Bi[kk]) + NHI*(-1 * beta[kk])
                                + (-1. * 10.**(NHI-Nc) / np.log(10)))  # log10 [ exp(-NHI/Nc) ]
            # f(z)
            fz = np.zeros((lenz, 2))
            # Loop on NHI
            itypes = ['LAF', 'DLA']
            for kk in range(ncomp):
                if kk == 0:  # LyaF
                    zcuts = self.param['LAF']['zcuts']
                    gamma = self.param['LAF']['gamma']
                else:        # DLA
                    zcuts = self.param['DLA']['zcuts']
                    gamma = self.param['DLA']['gamma']
                zcuts = [0] + zcuts + [999.]
                Aval = self.param[itypes[kk]]['Aval']
                # Cut on z
                for ii in range(1,len(zcuts)):
                    izcut = np.where( (z_val < zcuts[ii]) &
                                      (z_val > zcuts[ii-1]) )[0]
                    liz = len(izcut)
                    # Evaluate (at last!)
                    if (ii <=2) & (liz > 0):
                        fz[izcut, kk] = Aval * ( (1+z_val[izcut]) /
                                                 (1+zcuts[1]) )**gamma[ii-1]
                    elif (ii == 3) & (liz > 0):
                        fz[izcut, kk] = Aval * ( ( (1+zcuts[2]) /
                                                   (1+zcuts[1]) )**gamma[ii-2] *
                                                    ((1+z_val[izcut]) / (1+zcuts[2]) )**gamma[ii-1] )
            # dX/dz
            if cosmo is None:
                cosmo = self.cosmo
            dXdz = pyigmu.cosm_xz(z_val, cosmo=cosmo, flg_return=1)

            # Final steps
            if flg_1D == 1:
                fnX = np.sum(fz * 10.**log_gN, 1) / dXdz
                log_fNX = np.log10(fnX)
            else: 
                # Generate the matrix
                fnz = np.zeros((lenNHI, lenz))
                for kk in range(ncomp):
                    fnz += np.outer(10.**log_gN[:, kk], fz[:, kk])
                # Finish up
                log_fNX = np.log10(fnz) - np.log10( np.outer(np.ones(lenNHI), dXdz))
        elif self.mtype == 'PowerLaw':
            log_fNX = self.param['B'] + self.param['beta'] * NHI
            #
            if (not isiterable(z_val)) | (flg_1D == 1):  # scalar or 1D array wanted
                log_fNX += self.gamma * np.log10((1+z_val)/(1+self.zpivot))
            else:
                lgNHI_grid = np.outer(log_fNX, np.ones(len(z_val)))
                lenfX = len(log_fNX)
                #
                z_grid1 = 10**(np.outer(np.ones(lenfX)*self.gamma,
                                        np.log10(1+z_val)))  #; (1+z)^gamma
                z_grid2 = np.outer( np.ones(lenfX)*((1./(1+self.zpivot))**self.gamma),
                                    np.ones(len(z_val)))
                log_fNX = lgNHI_grid + np.log10(z_grid1*z_grid2)
        else:
            raise ValueError('fN.model: Not ready for this model type {:%s}'.format(self.mtype))

        # Return
        if (lenNHI + lenz) == 2:
            return log_fNX.flatten()[0]  # scalar
        else:
            return log_fNX

    def mfp(self, zem, neval=5000, cosmo=None, zmin=0.6):
        """ Calculate mean free path

        Parameters
        ----------
        zem : float
          Redshift of source
        cosmo : astropy.cosmology, optional
          Cosmological model to adopt (as needed)
        neval : int, optional
          Discretization parameter (5000)
        zmin: float, optional
          Minimum redshift in the calculation (0.5)

        Returns
        -------
        mfp : Quantity
          Mean free path from zem (physical Mpc)
        """
        # Imports
        from pyigm.fN import tau_eff as pyteff

        # Cosmology
        if cosmo is None:
            cosmo = cosmology.core.FlatLambdaCDM(70., 0.3)

        # Calculate teff
        zval, teff_LL = pyteff.lyman_limit(self, zmin, zem, N_eval=neval, cosmo=cosmo)

        # Find tau=1
        zt_interp = scii.interp1d(teff_LL, zval)
        ztau1 = zt_interp(1.)  # Will probably break if 1 is not covered

        # MFP
        mfp = np.fabs(cosmo.lookback_distance(ztau1) -
                        cosmo.lookback_distance(zem))  # Mpc
        # Return
        return mfp


    ##
    # Output
    def __repr__(self):
        return ('<{:s}: {:s} zmnx=({:g},{:g})>'.format(
                self.__class__.__name__, self.mtype, self.zmnx[0], self.zmnx[1]))


