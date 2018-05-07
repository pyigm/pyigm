""" Module containing the ClusteringField class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from scipy.ndimage import gaussian_filter as gf

from astropy.table import Table

from linetools import utils as ltu

from pyigm.field.igmfield import IgmGalaxyField
from pyigm.clustering.xcorr_utils import spline_sensitivity, random_gal, auto_pairs_rt, cross_pairs_rt
from pyigm.clustering.xcorr_utils import random_abs_zmnx


class ClusteringField(IgmGalaxyField):
    """The Field class is meant to contain information from a given
    galaxy field (with one or more QSO sightlines containing IGM
    information).  It sub-classes on IGMGalaxyField which is expected
    to hold *all* of the galaxies for the field of interest.

    Input parameters:
    ---
    radec : coordinate
    absreal:   catalog of real absorption lines (numpy rec array). It has to have
               dtype.names RA,DEC,ZABS,LOGN,B
    galreal:   catalog of real galaxies (numpy rec array). It has to have
               dtype.names RA, DEC, ZGAL, MAG
    wa:        numpy array with QSO spectral coverage (wavelength)
    fl:        numpy aray with the QSO spectrum flux it can be normalized or not)
    er:        numpy array with the QSO spectrum error (it can be normalized
               or not, but has to be consistent with fl)
    R:         resolution of the QSO spectrum spectrograph
    Ngal_rand: number of random galaxies per real one that will be created.
    Nabs_rand: number of random absorbers per real one that will be created.
    proper:    calculates everything in physical Mpc rather than co-moving (Boolean).
    wrapped : bool, optional
      Whether to wrap the pair counts along the line of sight to increase counts
        (if the signal is isotropic this is justified)


    Description of the Class:

    It first creates absrand and galrand using random_abs() and
    random_gal() functions. CRA and CDEC are the center position of the
    field, and so will define the coordinate system. Galaxies and
    absorbers (RA,DEC,z) are transformed to (X,Y,Z) co-moving
    coordinates (assuming proper=False). It will then calculate
    cross-pairs DaDg,DaRg,RaDg,RaRg and auto-pairs DgDg, DgRg, RgRg and
    DaDa, DaRa, RaRa where D means 'data' R means 'random', a means
    'absorber' and g means 'galaxy'. The class has also some
    implemented plots that are useful.

    """
    def __init__(self, radec, Ngal_rand=10, Nabs_rand=100, proper=False, field_name='',
                 wrapped=True, **kwargs):

        IgmGalaxyField.__init__(self, radec, **kwargs)

        self.Ngal_rand = Ngal_rand  # Ngal_rand x len(galreal) = NRANDOM (Gal)
        self.Nabs_rand = Nabs_rand  # Nabs_rand x len(absreal) = NRANDOM (Abs)

        self.absreal = None # absreal  # np rec array with absorber properties (single ion)
        #self.absrand = random_abs(self.absreal, self.Nabs_rand, wa, fl, er, R=R)

        # Central coordiantes
        self.CRA = self.coord.ra.value #np.mean(self.absreal.RA)
        self.CDEC = self.coord.dec.value # np.mean(self.absreal.DEC)

        if len(field_name) == 0:
            self.name = ltu.name_from_coord(self.coord)
        else:
            self.name = field_name
        self.proper = proper
        self.wrapped = wrapped

        self.galreal = None # Filled with addGal
        self.absreal = None # Filled with addAbs
        self.ion = None

    def addGal(self, gal_idx, mag_clm='MAG', z_clm='ZGAL', sens_galaxies=None,
               magbins=None, SPL=None):
        """ Adds a set of galaxies for clustering analysis from the main astropy Table

        A random sample is also generated

        Parameters
        ----------
        gal_idx : ndarray (int or bool)
          Indices of self.galaxy to add to the clustering analysis
        mag_clm : str, optional
          Name of the column for the galaxy magnitudes
        z_clm : str, optional
          Name of the column for the galaxy redshifts
        sens_galaxies : np.recarray
          Use for generating sensitivity function instead of the added galaxies
          Used when the subset is too small for an accurate sensitivity function
        magbins : ndarray (optional)
        SPL : dict (optional)
          Contains the sensitivity function (CubicSpline's)

        Returns
        -------
        Updates self.galreal internally

        """
        # Grab the galaxies
        sub_gal = self.galaxies[gal_idx]

        # Rename any columns here
        sub_gal.rename_column(mag_clm,'MAG')
        sub_gal.rename_column(z_clm,'ZGAL')

        # Strip down to important columns only here, as desired

        # Convert to numpy rec array
        galnew = sub_gal.as_array().view(np.recarray)

        # Calculate randoms
        if sens_galaxies is None:
            sens_galaxies = galnew
        # Sensitivity function
        if (SPL is None) or (magbins is None):
            magbins, _, SPL = spline_sensitivity(sens_galaxies)
        # Randoms
        rgal = random_gal(galnew, self.Ngal_rand, magbins, SPL)

        # Load me up
        if self.galreal is None:
            self.galreal = galnew  # np rec array with galaxy properties
            self.galrand = rgal
        else:
            galnew = galnew.astype(self.galreal.dtype)
            self.galreal = np.append(self.galreal, galnew)
            self.galreal = np.rec.array(self.galreal)
            self.galrand = np.append(self.galrand, rgal)
            self.galrand = np.rec.array(self.galrand)

    def addAbs(self, abs_input, zmnx, wrest, z_clm='ZABS'):
        # Convert to rec array, if need be
        if isinstance(abs_input, Table):
            # Rename any columns here
            abs_input.rename_column(z_clm,'ZABS')
            # Rename any columns here
            absnew = abs_input.as_array().view(np.recarray)
        else:  # Assuming rec array (what else would it be?)
            absnew = abs_input

        # Randoms
        absrand = random_abs_zmnx(absnew, self.Nabs_rand, zmnx, wrest)


        if self.absreal is None:
            self.absreal = absnew  # np rec array with galaxy properties

    def compute_pairs(self, tbinedges, rbinedges):
        """Computes relevant pairs dependent on what has been loaded
        into the object (i.e. galaxies, absorbers, galaxy-absorbers)

        Parameters
        ----------
        tbinedges : ndarray
          Defines the transverse bins in comoving Mpc
        rbinedges
          Defines the radial bins in comoving Mpc

        Returns
        -------
        Internal updates only, e.g. DgDg, RgRg

        """

        self.tbinedges = tbinedges
        self.rbinedges = rbinedges


        # Gal-gal only
        if self.galreal is not None:
            self.XYZ_gal()
            if self.proper:
                self.XYZ_proper()
            self.DgDg = auto_pairs_rt(self.xg, self.yg, self.zg, rbinedges, tbinedges, wrap=self.wrapped)
            self.RgRg = auto_pairs_rt(self.xgr, self.ygr, self.zgr, rbinedges, tbinedges, wrap=self.wrapped)
            self.DgRg = cross_pairs_rt(self.xg, self.yg, self.zg, self.xgr, self.ygr, self.zgr, rbinedges, tbinedges,
                                       wrapped=self.wrapped)

        # Abs-abs only
        if self.absreal is not None:
            self.XYZ_abs()
            if self.proper:
                self.XYZ_proper()
            self.DaDa = auto_pairs_rt(self.xa, self.ya, self.za, rbinedges, tbinedges, wrap=self.wrapped)
            self.RaRa = auto_pairs_rt(self.xar, self.yar, self.zar, rbinedges, tbinedges, wrap=self.wrapped)
            self.DaRa = cross_pairs_rt(self.xa, self.ya, self.za, self.xar, self.yar, self.zar, rbinedges, tbinedges,
                                       wrapped=self.wrapped)


            '''
            self.DaDg = cross_pairs_rt(self.xa, self.ya, self.za, self.xg, self.yg, self.zg, rbinedges, tbinedges,
                                       wrapped=wrapped)

            self.DaRg = cross_pairs_rt(self.xa, self.ya, self.za, self.xgr, self.ygr, self.zgr, rbinedges, tbinedges,
                                       wrapped=wrapped)
            self.RaRg = cross_pairs_rt(self.xar, self.yar, self.zar, self.xgr, self.ygr, self.zgr, rbinedges, tbinedges,
                                       wrapped=wrapped)
            self.RaDg = cross_pairs_rt(self.xar, self.yar, self.zar, self.xg, self.yg, self.zg, rbinedges, tbinedges,
                                       wrapped=wrapped)
            '''

    def XYZ_gal(self):
        """Calculates X,Y,Z coordinates for the galaxies, real and random
        from (RA,DEC,Z).

        Note: This assumes a small angle approximation which will not apply for
          very large separations.

        Internals filled only:
        xg,yg,zg
        xgr,ygr,zgr
        """
        from pyigm.clustering.coord_utils import radec_to_xyz

        # randomize order
        #self.galreal = shuffle_rec_array(self.galreal)
        #self.galrand = shuffle_rec_array(self.galrand)

        # Real galaxies
        #Rlos = np.array([self.cosmo.comoving_distance(self.galreal.ZGAL[i]) / PC / 1e6 for i in xrange(self.galreal.size)])
        Rlos = self.cosmo.comoving_distance(self.galreal.ZGAL).to('Mpc').value
        xyz = radec_to_xyz((self.galreal.RA - self.CRA) * np.cos(self.CDEC * np.pi / 180.),
                           self.galreal.DEC - self.CDEC).T * Rlos
        self.xg = xyz[0]
        self.yg = xyz[1]
        self.zg = xyz[2]

        # Random galaxies
        #Rlos = np.array([self.cosmo.Dc(self.galrand.ZGAL[i]) / PC / 1e6 for i in xrange(self.galrand.size)])
        Rlos = self.cosmo.comoving_distance(self.galrand.ZGAL).to('Mpc').value
        xyz = radec_to_xyz((self.galrand.RA - self.CRA) * np.cos(self.CDEC * np.pi / 180.),
                           self.galrand.DEC - self.CDEC).T * Rlos
        self.xgr = xyz[0]
        self.ygr = xyz[1]
        self.zgr = xyz[2]

    def XYZ_abs(self):
        """Calculates X,Y,Z coordinates for the galaxies, real and random
        from (RA,DEC,Z).

        Internals filled only:
        xg,yg,zg
        xgr,ygr,zgr
        """
        from pyigm.clustering.coord_utils import radec_to_xyz

        # randomize order
        #self.absreal = shuffle_rec_array(self.absreal)
        #self.absrand = shuffle_rec_array(self.absrand)

        #Rlos = np.array([cosmo.Dc(self.absreal.ZABS[i]) / PC / 1e6 for i in xrange(self.absreal.size)])
        Rlos = self.cosmo.comoving_distance(self.absreal.ZABS).to('Mpc').value
        xyz = radec_to_xyz((self.absreal.RA - self.CRA) * np.cos(self.CDEC * np.pi / 180.),
                           self.absreal.DEC - self.CDEC).T * Rlos
        self.xa = xyz[0]
        self.ya = xyz[1]
        self.za = xyz[2]

        #Rlos = np.array([cosmo.Dc(self.absrand.ZABS[i]) / PC / 1e6 for i in xrange(self.absrand.size)])
        Rlos = self.cosmo.comoving_distance(self.absrand.ZABS).to('Mpc').value
        xyz = radec_to_xyz((self.absrand.RA - self.CRA) * np.cos(self.CDEC * np.pi / 180.),
                           self.absrand.DEC - self.CDEC).T * Rlos
        self.xar = xyz[0]
        self.yar = xyz[1]
        self.zar = xyz[2]


    def XYZ_proper(self):
        """From co-moving coordinates goes to physical coordinates
        (dividing by 1+z). It is only meaningful at small distances
        (<~100 kpc)"""
        self.yg = self.yg / (1. + self.galreal.ZGAL)
        self.zg = self.zg / (1. + self.galreal.ZGAL)
        self.ygr = self.ygr / (1. + self.galrand.ZGAL)
        self.zgr = self.zgr / (1. + self.galrand.ZGAL)
        if self.absreal is not None:
            self.ya = self.ya / (1. + self.absreal.ZABS)
            self.za = self.za / (1. + self.absreal.ZABS)
            self.yar = self.yar / (1. + self.absrand.ZABS)
            self.zar = self.zar / (1. + self.absrand.ZABS)

    # Output
    def __repr__(self):
        rstr = '<{:s}: field={:s} '.format(
            self.__class__.__name__,
            self.name)

        if self.galreal is not None:
            rstr += 'ngreal={:d} '.format(len(self.galreal))
        rstr += '>'
        return rstr
