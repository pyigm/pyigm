"""Module containing 2PCF Survey Class"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
from scipy.ndimage import gaussian_filter as gf

from astropy.coordinates import SkyCoord
from astropy import units as u

from pyigm.clustering.xcorr_utils import random_gal, auto_pairs_rt, cross_pairs_rt, W3, collapse_along_LOS

class Survey2D2PCF(object):
    """
    Class for 2D-2PCF for absorber-galaxy two-point
    cross and auto-correlations. The Class loads several independent
    IgmGalaxyField objects together. It internally calculates the total numbers
    of cross- and auto-pairs of real data and produce randoms. It calculates
    the 2D2PCF using the Landy & Szalay estimator (or others), and it estimates
    the uncertainty with a bootstrap, jacknife or Landy & Szalay approximation.

    Parameters
    ----------
    field : ClusteringField object
      The field which initializes this object sets the binning.  It also
      sets whether galaxy-galaxy, absorber-absorber, and galaxy-absorber analysis is on
      So, load a good one first!
    """

    def __init__(self, field):
        #f = copy.deepcopy(field)
        self.fields = [field]

        if field.galreal is not None:
            self.DgDg = field.DgDg
            self.DgRg = field.DgRg
            self.RgRg = field.RgRg
            self.Ngal = len(field.galreal)
            self.gal_anly = True
        else:
            self.gal_anly = False

        if field.absreal is not None:
            self.DaDa = field.DaDa
            self.DaRa = field.DaRa
            self.RaRa = field.RaRa
            self.Nabs = len(field.absreal)
            self.abs_anly = True
            # Both!
            if field.galreal is not None:
                self.DaDg = field.DaDg
                self.DaRg = field.DaRg
                self.RaDg = field.RaDg
                self.RaRg = field.RaRg
        else:
            self.abs_anly = False


        # Field edges
        self.rbinedges = field.rbinedges # Should check that any new field has the same edges
        self.tbinedges = field.tbinedges

        # Binning
        self.rcbins = 0.5*(self.rbinedges[:-1] + self.rbinedges[1:])
        self.tcbins = 0.5*(self.tbinedges[:-1] + self.tbinedges[1:])
        self.rdiff  = self.rbinedges[1:] - self.rbinedges[:-1]
        self.tdiff  = self.tbinedges[1:] - self.tbinedges[:-1]

    def addField(self, new_field, SEP_TOL=0.5*u.deg):
        """ Add a new field

        Parameters
        ----------
        new_field : ClusteringField
        SEP_TOL : Angle
          Separation tolerance.  If the new field is within SEP_TOL
          of any existing ones, the code will barf

        Returns
        -------
        All internal appends

        """
        # Check edges
        if not np.all(np.isclose(new_field.rbinedges, self.rbinedges)):
            raise IOError("rbinedges of new field does not match Survey!!")
        if not np.all(np.isclose(new_field.tbinedges, self.tbinedges)):
            raise IOError("tbinedges of new field does not match Survey!!")
        # Check coordiantes (this is not a perfect check..)
        fcoord = SkyCoord(ra=[field.CRA for field in self.fields],
                          dec=[field.CDEC for field in self.fields], unit='deg')
        if np.any(new_field.coord.separation(fcoord) < SEP_TOL):
            raise IOError('New field is within {} of an already input field!  We assume a mistake was made'.format(SEP_TOL))


        #f = copy.deepcopy(field)
        self.fields.append(new_field)

        if new_field.galreal is not None:
            self.DgDg += new_field.DgDg
            self.DgRg += new_field.DgRg
            self.RgRg += new_field.RgRg
            self.Ngal += len(new_field.galreal)

        if new_field.absreal is not None:
            self.DaDa += new_field.DaDa
            self.DaRa += new_field.DaRa
            self.RaRa += new_field.RaRa
            self.Nabs += len(new_field.absreal)
            # Both?
            if new_field.galreal is not None:
                self.DaDg += new_field.DaDg
                self.DaRg += new_field.DaRg
                self.RaDg += new_field.RaDg
                self.RaRg += new_field.RaRg

    def calc_xi_gg(self, **kwargs):
        """ Simple wrapper to calc_xi
        Parameters
        ----------
        kwargs

        Returns
        -------

        """
        W, err_W = self.calc_xi('galaxy-galaxy', **kwargs)
        # Set em
        self.xi_gg = W
        self.xi_gg_err_ls = err_W
        # Return them too
        return W, err_W

    def calc_xi_ag(self, **kwargs):
        """ Simple wrapper to calc_xi
        Parameters
        ----------
        kwargs

        Returns
        -------

        """
        W, err_W = self.calc_xi('absorber-galaxy', **kwargs)
        # Set em
        self.xi_ag = W
        self.xi_ag_err_ls = err_W
        # Return them too
        return W, err_W

    def calc_xi(self, xi_type, sigma=0, error=None):
        """ Calculate xi_gg with the Landy-Szalay estimator (W3)

        Parameters
        ----------
        xi_type : str
          galaxy-galaxy

        sigma : float, optional
          Gaussian standard deviation smoothing;
          defined in units of the spatial grid.
        error : str, optional
          Type of error analysis
          'jk' -- jackknife
          'bs' -- bootstrap

        Returns
        -------
        Wgg : ndarray
        err_W : ndarray
        xi_gg and xi_gg_err_ls are also set internally

        """
        if xi_type == 'galaxy-galaxy':
            xa, xb = 'g', 'g'
            auto = True
        elif xi_type == 'absorber-galaxy':
            xa, xb = 'a', 'g'
            auto = False
        else:
            raise IOError("Not ready for xi_type={:s}".format(xi_type))

        # Wrappers
        def DD(obj):
            return getattr(obj,'D{:s}D{:s}'.format(xa,xb))
        def RR(obj):
            return getattr(obj,'R{:s}R{:s}'.format(xa,xb))
        def DR(obj):
            return getattr(obj,'D{:s}R{:s}'.format(xa,xb))
        def RD(obj):
            if auto:
                return DR(obj)
            else:
                return getattr(obj,'R{:s}D{:s}'.format(xa,xb))
        # Number
        nDD = getattr(self,'nD{:s}D{:s}'.format(xa,xb))
        nRR = getattr(self,'nR{:s}R{:s}'.format(xa,xb))
        nDR = getattr(self,'nD{:s}R{:s}'.format(xa,xb))
        if auto:
            nRD = nDR
        else:
            nRD = getattr(self,'nR{:s}D{:s}'.format(xa,xb))

        # Calculate me
        s = sigma
        W, err_W = W3(gf(DD(self), s), gf(RR(self), s), gf(DR(self), s), gf(RD(self), s),
                        Ndd=nDD, Nrr=nRR, Ndr=nDR, Nrd=nRD)
        #

        # jacknife error
        if error == 'jk':
            err_Wjk = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
            psd_val = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
            N = len(self.fields)
            for field in self.fields:
                DD_aux = DD(self) - DD(field) #self.DgDg - field.DgDg
                RR_aux = RR(self) - RR(field) #self.RgRg - field.RgRg
                DR_aux = DR(self) - DR(field) #self.DgRg - field.DgRg
                if auto:
                    RD_aux = DR_aux
                else:
                    RD_aux = RD(self) - RD(field) #self.DgRg - field.DgRg
                W_aux, _ = W3(gf(DD_aux, s), gf(RR_aux, s), gf(DR_aux, s), gf(RD_aux, s),
                                Ndd=nDD, Nrr=nRR, Ndr=nDR, Nrd=nRD)
                                #Ndd=self.nDgDg, Nrr=self.nRgRg, Ndr=self.nDgRg, Nrd=self.nDgRg)
                psd_val += N * W - (N - 1) * W_aux
            mean_psd = psd_val / N
            for field in self.fields:
                DD_aux = DD(self) - DD(field) #self.DgDg - field.DgDg
                RR_aux = RR(self) - RR(field) #self.RgRg - field.RgRg
                DR_aux = DR(self) - DR(field) #self.DgRg - field.DgRg
                if auto:
                    RD_aux = DR_aux
                else:
                    RD_aux = RD(self) - RD(field) #self.DgRg - field.DgRg
                W_aux, _ = W3(gf(DD_aux, s), gf(RR_aux, s), gf(DR_aux, s), gf(RD_aux, s),
                                Ndd=nDD, Nrr=nRR, Ndr=nDR, Nrd=nRD)
                err_Wjk += (mean_psd - W_aux) ** 2

            err_Wjk = err_Wjk / N / (N - 1)
            err_Wjk = np.sqrt(err_Wjk)
            err_W = err_Wjk

        if error == 'bs':
            Nf = len(self.fields)
            err_Wbt = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
            W_sum = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
            W_sum2 = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)

            N = self.Nbt
            for i in range(N):
                inds = np.random.randint(0, Nf, Nf)
                DD_aux = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
                RR_aux = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
                DR_aux = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
                RD_aux = np.zeros((len(self.rbinedges) - 1, len(self.tbinedges) - 1), float)
                for ind in inds:
                    auxfield = self.fields[ind]
                    DD_aux += DD(auxfield) #auxfield.DgDg
                    RR_aux += RR(auxfield) #auxfield.RgRg
                    DR_aux += DR(auxfield) #auxfield.DgRg
                    if auto:
                        RD_aux = DR_aux
                    else:
                        RD_aux += RD(auxfield) #auxfield.DgRg
                W_aux, _ = W3(gf(DD_aux, s), gf(RR_aux, s), gf(DR_aux, s), gf(RD_aux, s),
                                Ndd=nDD, Nrr=nRR, Ndr=nDR, Nrd=nRD)
                                #Ndd=self.nDgDg, Nrr=self.nRgRg, Ndr=self.nDgRg, Nrd=self.nDgRg)
                W_sum += W_aux
                W_sum2 += W_aux ** 2
            W_sum = W_sum / N
            W_sum2 = W_sum2 / N
            err_Wbt = np.sqrt(W_sum2 - W_sum ** 2)
            err_W = err_Wbt

        # Set
        # Return
        return W, err_W

    def calc_xi_gg_transverse(self, s, nbins=None):
        """  Calculate xi_gg in the transverse dimension
        Estimated with Landay-Szalay  (W3)

        Parameters
        ----------
        s : float
          Smoothing parameter
        nbins : int, optional

        Returns
        -------
        xi_gg_T
        xi_gg_T_err
          Also save in the object

        """

        # Transverse pairs (may wish to make this a method)
        self.DgDg_T = collapse_along_LOS(self.DgDg,nbins,s=s)
        self.DgRg_T = collapse_along_LOS(self.DgRg,nbins,s=s)
        self.RgRg_T = collapse_along_LOS(self.RgRg,nbins,s=s)


        # Calculate
        self.xi_gg_T,self.xi_gg_T_err = W3(self.DgDg_T,self.RgRg_T,self.DgRg_T,self.DgRg_T,
                                           Ndd=self.nDgDg,Nrr=self.nRgRg,Ndr=self.nDgRg,Nrd=self.nDgRg)
        pdb.set_trace()

        return self.xi_gg_T, self.xi_gg_T_err

    def calc_xi_ag_transverse(self, s, nbins=None):
        """  Calculate xi_ag in the transverse dimension
        Estimated with Landay-Szalay  (W3)

        Parameters
        ----------
        s : float
          Smoothing parameter
        nbins : int, optional

        Returns
        -------
        xi_ag_T
        xi_ag_T_err
          Also save in the object

        """

        # Transverse pairs (may wish to make this a method)
        self.DaDg_T = collapse_along_LOS(self.DaDg,nbins,s=s)
        self.DaRg_T = collapse_along_LOS(self.DaRg,nbins,s=s)
        self.RaDg_T = collapse_along_LOS(self.RaDg,nbins,s=s)
        self.RaRg_T = collapse_along_LOS(self.RaRg,nbins,s=s)

        # Calculate
        self.xi_ag_T,self.xi_ag_T_err = W3(self.DaDg_T,self.RaRg_T,self.DaRg_T,self.RaDg_T,
                                           Ndd=self.nDaDg,Nrr=self.nRaRg,Ndr=self.nDaRg,Nrd=self.nRaDg)

        return self.xi_ag_T, self.xi_ag_T_err

    def set_normalization(self, norm, Ngal_rand=None, Nabs_rand=None):
        """ Set normalization for the pair counts

        Parameters
        ----------
        norm : bool
          False -- Internal normalization will be performed by the estimator
        Ngal_rand : int, optional
        Nabs_rand : int, optional

        Returns
        -------
        Attributes like nDgDg, nRgRg are set internally
        """

        # normalization factors
        if norm:
            if self.gal_anly:
                self.nDgDg = self.Ngal * (self.Ngal - 1) / 2.
                self.nDgRg = self.Ngal * self.Ngal * Ngal_rand
                self.nRgRg = self.Ngal * Ngal_rand * (self.Ngal * Ngal_rand - 1) / 2.

            if self.abs_anly:
                self.nDaDa = self.Nabs * (self.Nabs - 1) / 2.
                self.nDaRa = self.Nabs * self.Nabs * Nabs_rand
                self.nRaRa = self.Nabs * Nabs_rand * (self.Nabs * Nabs_rand - 1) / 2.

            if self.gal_anly & self.abs_anly:
                self.nRaDg = self.Nabs * self.Ngal * Nabs_rand
                self.nRaRg = self.Nabs * self.Ngal * Nabs_rand * Ngal_rand
                self.nDaDg = self.Nabs * self.Ngal
                self.nDaRg = self.Nabs * self.Ngal * Ngal_rand
        else:
            self.nDgDg = None
            self.nDaDa = None
            self.nDaDg = None
            self.nDaRg = None
            self.nRaDg = None
            self.nRaRg = None
            self.nRgRg = None
            self.nRaRa = None
            self.nDaRa = None
            self.nDgRg = None

    # Output
    def __repr__(self):
        rstr = '<{:s}: nfield={:d} gal_anly={:s} abs_analy={:s} \n'.format(
            self.__class__.__name__,
            len(self.fields),
            str(self.gal_anly),
            str(self.abs_anly))

        for field in self.fields:
            rstr += 'field at ({},{}) \n'.format(field.coord.ra.value, field.coord.dec.value)
        rstr += '>'
        return rstr

