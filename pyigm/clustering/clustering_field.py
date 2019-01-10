""" Module containing the ClusteringField class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15

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
      Sets self.coord (in IgmGalaxyfield)
    absreal:   catalog of real absorption lines (numpy rec array). It has to have
               dtype.names RA,DEC,ZABS,LOGN,B
    galreal:   catalog of real galaxies (numpy rec array). It has to have
               dtype.names RA, DEC, ZGAL, MAG
    wa:        numpy array with QSO spectral coverage (wavelength)
    fl:        numpy aray with the QSO spectrum flux it can be normalized or not)
    er:        numpy array with the QSO spectrum error (it can be normalized
               or not, but has to be consistent with fl)
    R:         resolution of the QSO spectrum spectrograph
    igal_rand: number of random galaxies per real one that will be created.
    iabs_rand: number of random absorbers per real one that will be created.
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
    def __init__(self, radec, igal_rand=10, iabs_rand=100, proper=False, field_name='',
                 wrapped=True, cosmo=None, **kwargs):

        IgmGalaxyField.__init__(self, radec, **kwargs)

        self.igal_rand = igal_rand
        self.iabs_rand = iabs_rand

        self.absreal = None # absreal  # np rec array with absorber properties (single ion)
        self.absrand = None

        # Cosmology
        if cosmo is None:
            cosmo = Planck15
        self.cosmo = cosmo

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

    @property
    def Ngal(self):
        if self.galreal is None:
            return 0
        else:
            return len(self.galreal)

    @property
    def Ngal_rand(self):
        if self.galrand is None:
            return 0
        else:
            return len(self.galrand)

    @property
    def Nabs(self):
        if self.absreal is None:
            return 0
        else:
            return len(self.absreal)

    @property
    def Nabs_rand(self):
        if self.absrand is None:
            return 0
        else:
            return len(self.absrand)

    def addGal(self, gal_idx, mag_clm='MAG', z_clm='ZGAL', sens_galaxies=None,
               magbins=None, SPL=None, Rcom_max=None, debug=False, tbinedges=None):
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
          array containing the shape of the sensitivity function of the galaxies being added
          Used for generating sensitivity function instead of the added galaxies
          Important to use when the subset is too small for an accurate sensitivity function
        magbins : ndarray (optional)
        SPL : dict (optional)
          Contains the sensitivity function (CubicSpline's)
        Rcom_max : float, optional
          Maximum comoving separation of the survey in Mpc
        tbinedges : ndarray
          Only used for debug=True

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
        galnewrand = random_gal(galnew, self.igal_rand, magbins, SPL)

        # Cut on Rcom_max?
        if Rcom_max is not None:
            rcoord = SkyCoord(ra=galnewrand.RA, dec=galnewrand.DEC, unit='deg')
            angsep = self.coord.separation(rcoord)
            # R comoving
            Rcom = self.cosmo.kpc_comoving_per_arcmin(galnewrand.ZGAL) * angsep.to('arcmin')
            # Cut
            goodr = Rcom.to('Mpc').value < Rcom_max
            galnewrand = galnewrand[goodr]
        if debug:
            gcoord = SkyCoord(ra=galnew.RA, dec=galnew.DEC, unit='deg')
            gangsep = self.coord.separation(gcoord)
            gRcom = self.cosmo.kpc_comoving_per_arcmin(galnew.ZGAL) * gangsep.to('arcmin')
            faintR = galnewrand.MAG > 19.
            faintg = galnew.MAG > 19.
            #
            from matplotlib import pyplot as plt
            import matplotlib.gridspec as gridspec
            plt.clf()
            if True:
                gs = gridspec.GridSpec(2,1)
                ax = plt.subplot(gs[0])
                ax.hist(gRcom[~faintg], color='k', bins=tbinedges, normed=1, label='DD', fill=False)
                ax.hist(Rcom[goodr][~faintR], edgecolor='red', bins=tbinedges, normed=1, label='RR', fill=False)
                # Faint
                ax = plt.subplot(gs[1])
                ax.hist(gRcom[faintg], color='k', bins=tbinedges, normed=1, label='DD', fill=False)
                ax.hist(Rcom[goodr][faintR], edgecolor='red', bins=tbinedges, normed=1, label='RR', fill=False)
                ax.set_ylabel('Faint')
                ax.set_xlabel('Rcom (Mpc)')
            else:
                ax = plt.gca()
                zbins = np.arange(0., 0.8, 0.025)
                ax.hist(galnew.ZGAL[faintg], color='k', bins=zbins, normed=1, label='DD', fill=False)
                ax.hist(galnewrand.ZGAL[faintR], edgecolor='red', bins=zbins, normed=1, label='RR', fill=False)
                ax.set_xlabel('zGAL')
            plt.show()

        # Load me up
        if self.galreal is None:
            self.galreal = galnew  # np rec array with galaxy properties
            self.galrand = galnewrand
        else:
            galnew = galnew.astype(self.galreal.dtype)
            self.galreal = np.append(self.galreal, galnew)
            self.galreal = np.rec.array(self.galreal)
            self.galrand = np.append(self.galrand, galnewrand)
            self.galrand = np.rec.array(self.galrand)


    def addAbs(self, abs_input, zmnx, wrest, z_clm='ZABS'):
        # Convert to rec array, if need be
        if isinstance(abs_input, Table):
            # Rename any columns here
            abs_input.rename_column(z_clm,'ZABS')
            absnew = abs_input.as_array().view(np.recarray)
        else:  # Assuming rec array (what else would it be?)
            absnew = abs_input

        # Randoms
        absrand = random_abs_zmnx(absnew, self.iabs_rand, zmnx, wrest)


        if self.absreal is None:
            self.absreal = absnew  # np rec array with galaxy properties
            self.absrand = absrand
        else:
            absnew = absnew.astype(self.absreal.dtype)
            self.absreal = np.append(self.absreal, absnew)
            self.absreal = np.rec.array(self.absreal)
            self.absrand = np.append(self.absrand, absrand)
            self.absrand = np.rec.array(self.absrand)


    def compute_pairs(self, tbinedges, rbinedges, debug=False):
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
            # Auto
            self.DaDa = auto_pairs_rt(self.xa, self.ya, self.za, rbinedges, tbinedges, wrap=self.wrapped)
            self.RaRa = auto_pairs_rt(self.xar, self.yar, self.zar, rbinedges, tbinedges, wrap=self.wrapped)
            self.DaRa = cross_pairs_rt(self.xa, self.ya, self.za, self.xar, self.yar, self.zar, rbinedges, tbinedges,
                                       wrapped=self.wrapped)
            #
            if self.galreal is not None:
                self.DaDg = cross_pairs_rt(self.xa, self.ya, self.za, self.xg, self.yg, self.zg, rbinedges, tbinedges,
                                           wrapped=self.wrapped, debug=debug)

                self.DaRg = cross_pairs_rt(self.xa, self.ya, self.za, self.xgr, self.ygr, self.zgr, rbinedges, tbinedges,
                                           wrapped=self.wrapped)
                self.RaRg = cross_pairs_rt(self.xar, self.yar, self.zar, self.xgr, self.ygr, self.zgr, rbinedges, tbinedges,
                                           wrapped=self.wrapped)
                self.RaDg = cross_pairs_rt(self.xar, self.yar, self.zar, self.xg, self.yg, self.zg, rbinedges, tbinedges,
                                           wrapped=self.wrapped)

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
        """Calculates X,Y,Z coordinates for the absorbers, real and random
        from (RA,DEC,Z).

        Internals filled only:
        xa,ya,za
        xar,yar,zar
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
        (<~200 kpc)"""
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

        rstr += 'ngreal={:d} '.format(self.Ngal)

        if self.absreal is not None:
            rstr += 'nareal={:d} '.format(len(self.absreal))
        rstr += '>'
        return rstr
