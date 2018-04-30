""" Module containing the ClusteringField class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from scipy.ndimage import gaussian_filter as gf

from pyigm.field.igmfield import IgmGalaxyField
from pyigm.clustering.xcorr_utils import random_gal, auto_pairs_rt, cross_pairs_rt


class ClusteringField(IgmGalaxyField):
    """The Field class is meant to contain information from a given
    galaxy field (with one or more QSO sightlines containing IGM
    information).

    Input parameters:
    ---
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
    def __init__(self, radec, R=20000, Ngal_rand=10, Nabs_rand=100, proper=False, field_name='',
                 wrapped=True, **kwargs):

        IgmGalaxyField.__init__(self, radec, **kwargs)

        self.Ngal_rand = Ngal_rand  # Ngal_rand x len(galreal) = NRANDOM (Gal)
        self.Nabs_rand = Nabs_rand  # Nabs_rand x len(absreal) = NRANDOM (Abs)

        self.absreal = None # absreal  # np rec array with absorber properties (single ion)
        #self.absrand = random_abs(self.absreal, self.Nabs_rand, wa, fl, er, R=R)

        # Central coordiantes
        self.CRA = self.coord.ra.value #np.mean(self.absreal.RA)
        self.CDEC = self.coord.dec.value # np.mean(self.absreal.DEC)

        self.name = field_name
        self.proper = proper
        self.wrapped = wrapped

        self.galreal = None # Filled with addGal

    def addGal(self, gal_idx):
        """ Adds a set of galaxies for clustering analysis from the main astropy Table

        A random sample is also generated

        Parameters
        ----------
        gal_idx : ndarray (int or bool)
          Indices of self.galaxy to add to the clustering analysis

        Returns
        -------
        Updates self.galreal internally

        """
        # Grab the galaxies
        sub_gal = self.galaxies[gal_idx]

        # Rename any columns here
        sub_gal.rename_column('mag','MAG')
        sub_gal.rename_column('Z','ZGAL')

        # Strip down to important columns only here, as desired

        # Convert to numpy rec array
        galnew = sub_gal.as_array().view(np.recarray)

        # Load me up
        if self.galreal is None:
            self.galreal = galnew  # np rec array with galaxy properties
            self.galrand = random_gal(self.galreal, self.Ngal_rand)
        else:
            galnew = galnew.astype(self.galreal.dtype)
            self.galreal = np.append(self.galreal, galnew)
            self.galreal = np.rec.array(self.galreal)
            aux = random_gal(galnew, self.Ngal_rand)
            self.galrand = np.append(self.galrand, aux)
            self.galrand = np.rec.array(self.galrand)

    def compute_pairs(self, tbinedges, rbinedges):
        """Computes relevant pairs"""

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
        """Goes from (RA,DEC,Z) to (X,Y,Z) coordinates"""
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
        """Goes from (RA,DEC,Z) to (X,Y,Z) coordinates"""
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
