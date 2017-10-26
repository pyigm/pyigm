""" Class for working with CUBA
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
import warnings as warn
import pdb

from scipy.integrate import simps
from scipy.interpolate import interp1d

from astropy import units as u
from astropy import constants as const

# Path
pyigm_path = imp.find_module('pyigm')[1]
Ryd = const.Ryd.to('eV', equivalencies=u.spectral())

class CUBA(object):
    """
    Class for CUBA analysis

    JXP on 13 Oct 2015

    Attributes
    ----------
    fits_path : str, optional
      Path to the FITS data files for COS-Halos
    z : ndarray
      Array of z values from CUBA file
    energy : Quantity array
      Array of energy values, sorted (eV); converted from wave
    wave : Quantity array
      Array of wavelength values (reverse order) from CUBA file
    Jnu : Quantity 2D array [energy,z]
      Array of Jnu values from CUBA file
    """
    # Initialize with a .dat file
    def __init__(self, cuba_file=None):

        # CUBA file
        if cuba_file is None:
            cuba_file = pyigm_path+'/data/euvb/cuba_uvbapr2011_q1g01.hiz.out'
        self.cuba_file = cuba_file

        # Read
        self.read_cuba()

    def read_cuba(self):
        """ Read in a CUBA file 
        """
        # File
        # Read
        print('read_cuba: Using CUBA file -- {:s}'.format(self.cuba_file))
        with open(self.cuba_file,'r') as f:
            lines = f.readlines()
        # Parse
        flg_z = 0
        idx = 0
        nlin = len(lines)
        wave = []
        for qq, line in enumerate(lines):
            if line.strip()[0] == '#':
                continue
            # First good line has the redshifts
            if flg_z == 0:
                flg_z = 1
                self.z = np.array([float(val) for val in line.strip().split()])
                jnu = np.zeros( (nlin-qq-1, len(self.z)))
            else:
                parse = [float(val) for val in line.strip().split()]
                wave.append(parse[0])
                jnu[idx, :] = parse[1:]
                idx += 1

        # Unique values
        uni_wv, idx_wv = np.unique(np.array(wave), return_index=True)
        Jnu = np.zeros( (len(uni_wv), len(self.z)))
        # Sort
        srt = np.argsort(1./uni_wv)
        for ii in range(len(self.z)):
            Jnu[:, ii] = jnu[idx_wv[srt], ii] * u.erg/u.s/u.cm**2

        # Finish with units
        self.wave = np.array(uni_wv[srt])*u.AA
        self.energy = self.wave.to('eV', equivalencies=u.spectral())
        self.Jnu = Jnu * u.erg/u.s/u.cm**2

    def phi(self, zval, min_energy=None):
        """Calculate photon flux from a given minimum energy

        Parameters:
        ----------
        zval : float
          Redshift for evaluation
        min_energy : Quantity or Quantity array, optional
          Minimum energy for the calculation
          Default -- 1Ryd
        """
        # Init
        E_MAX = 1e10*Ryd
        if min_energy is None:
            min_energy = Ryd
            print('cuba.phi: Assuming minimum energy = {:g}'.format(min_energy))
        # Grab Jnu at the input redshift
        jnu = self.zinterp_jnu(zval)
        # Setup for log integral
        log_energy = np.log10(self.energy.to('eV')/Ryd)
        # Cut out high/low energies
        blue_energy = (self.energy >= min_energy) & (self.energy <= E_MAX)
        integrand = 4*np.pi*jnu.cgs/const.h.cgs # Note the factor of 4 pi
        integrand[~blue_energy] = 0.0
        # Integrate
        phi = np.log(10.)*simps(integrand.value, log_energy)
        # Return with Units
        return phi / u.s / u.cm**2

    def logU(self, zval, nH=1/u.cm**3, min_energy=1*Ryd):
        """ Estimate the ionization parameter at the given redshift
        for a default density of 1cc.

        Parameters
        ----------
        zval : float
        nH : Quantity, optional
          Gas density
        min_energy : Quantity, optional

        Returns
        -------
        logU : float
         log10 of the Ionization parameter, defined as U = Phi/c nH
        """
        Phi = self.phi(zval, min_energy=min_energy)
        U = (Phi / const.c.cgs / nH).decompose()
        if U.unit != u.dimensionless_unscaled:
            raise IOError("Bad units in your logU calculation..")
        else:
            return np.log10(U.value)

    def zinterp_jnu(self, zval, use_nearest=False):
        """Interpolate the Jnu grid at a given redshift

        Parameters
        ----------
        zval : float
          Redshift
        use_nearest : bool, optional
          Use nearest redshift instead??
        """
        # Do not interpolate beyond limits
        minz = np.min(self.z)
        maxz = np.max(self.z)
        if zval < minz:
            warn.warning('Input z was lower than z grid')
            print('Using z={:g}'.format(minz))
            return self.Jnu[:, 0].flatten()
        if zval > maxz:
            warn.warning('Input z was larger than z grid')
            print('Using z={:g}'.format(maxz))
            return self.Jnu[:, -1].flatten()

        # Find nearest? 
        if use_nearest:
            idx = np.argmin(np.abs(self.z-zval))
            return self.Jnu[:, idx].flatten()

        # Interpolate
        nval = self.energy.shape[0]
        jnu = np.zeros(nval)
        for ii in range(nval):
            jnu[ii] = interp1d(self.z, self.Jnu[ii, ])(zval)
        return jnu * self.Jnu.unit
        #

    def plot(self, zval, xlim=None):
        """Show the CUBA spectrum (Ryd vs. log Jnu)

        Parameters
        ----------
        zval : float
          Redshift
        xlim : tuple, optional
          xmin, xmax (Ryd)
        """
        import matplotlib as mpl
        mpl.rcParams['font.family'] = 'stixgeneral'
        from matplotlib import pyplot as plt

        plt.clf()
        xval = self.energy / const.Ryd.to('eV', equivalencies=u.spectral())
        yval = np.log10(self.zinterp_jnu(zval).value)
        #
        plt.plot(xval, yval, 'k-')
        plt.xlabel('Energy (Ryd)')
        plt.ylabel(r'$\log J_\nu$')
        plt.ylim(-25., -19.)
        if xlim is not None:
            plt.xlim(xlim)
        #
        plt.show()

    def __repr__(self):
        return ('[{:s}: cuba_file={:s}]'.format(
                self.__class__.__name__, self.cuba_file))

