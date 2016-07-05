""" Class for correlating galaxies and IGM within a given field on the sky
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import warnings
import pdb

from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy.cosmology import Planck15

import pyigm.field.utils as pfu

#from xastropy.xutils import xdebug as xdb


class IgmGalaxyField(object):
    """ Class for a field associating galaxies to the IGM/CGM

    Parameters
    ----------
    radec : tuple or SkyCoord
        (RA,DEC) in deg or astropy.coordinate
    name : str; optional
        Default is set from sexagesimal coordinates
    cosmo : Cosmology; optional
        Default is astropy.cosmology.Plack15

    """

    # Initialize
    def __init__(self, radec, name=None, cosmo=None, verbose=False):
        # coord
        if isinstance(radec, (tuple)):
            self.coord = SkyCoord(ra=radec[0], dec=radec[1])
        elif isinstance(radec, SkyCoord):
            self.coord = radec

        # Field
        if name is None:
            self.name = 'IGMFIELD_J{:s}{:s}'.format(
                    self.coord.ra.to_string(unit=u.hour, pad=True),
                    self.coord.dec.to_string(pad=True, alwayssign=True))
        self.name = name

        # Cosmology
        if cosmo is None:
            cosmo = Planck15
            if verbose is True:
                print('IgmGalxyField: Using Planck15 cosmology')
        self.cosmo = cosmo

        # Init
        self.igm = None
        self.targets = None
        self.galaxies = None
        self.observing = None
        self.selection = None

    def calc_rhoimpact(self, obj, comoving=True, los_coord=None):
        """Calculate impact parameter from LOS RA/DEC for a set of objects

        Parameters
        ----------
        obj : Table
          (can be anything that takes 'Z', and 'RA','DEC' in degrees)
          Sources for calculation
        comoving : bool, optional
           If True then comoving, else physical
        los_radec : SkyCoord, optional
          Defaults to field RA/DEC


        Returns
        -------
        rho : Quantity (array usually)
          Impact parameter(s) in kpc
        """
        if los_coord is None:
            los_coord = self.coord
        # Coord
        if isinstance(obj['RA'], Quantity):
            ora = obj['RA']
            odec = obj['DEC']
        else:
            ora = obj['RA']*u.deg
            odec = obj['DEC']*u.deg

        o_coord = SkyCoord(ra=ora, dec=odec)
        ang_sep = o_coord.separation(los_coord).to('arcmin')
        # Cosmology (kpc per arcmin)
        if comoving:
            kpc_amin = self.cosmo.kpc_comoving_per_arcmin(obj['Z'])
        else:
            kpc_amin = self.cosmo.kpc_proper_per_arcmin(obj['Z'])
        rho = ang_sep * kpc_amin
        # Return
        return rho

    def get_associated_galaxies(self, z, los_coord=None, R=300*u.kpc, dv_tol=500*u.km/u.s):
        """Return a Table of associated galaxies for a given redshift and separation

        Parameters
        ----------
        z : float
          Redshift for association (usually IGM/CGM absorber redshift)
        R : Quantity
          Radius of impact parameter for association [300kpc]
        dv_tol : Quantity
          Velocity window for association [500km/s]
        los_coord : SkyCoord, optional
          Line-of-sight coordinates

        Returns
        -------
        assoc_gal : Table
          Table of associated galaxies (if any)
        rho : Quantity
          Impact parameters [kpc]
        """
        # los_radec
        if los_coord is None:
            los_coord = self.coord
        # Cut on z first
        dv_gal = const.c.to('km/s') * (self.galaxies['Z']-z)/(1+z)  # Approximate
        gdz = np.where(np.abs(dv_gal)<dv_tol)[0]
        if len(gdz) == 0:
            return None
        #
        gdz_gal = self.galaxies[gdz]
        rho = self.calc_rhoimpact(gdz_gal, los_coord)  # Could add this to Table
        #
        gd_rho = np.where(rho < R)[0]
        if len(gd_rho) == 0:
            return None
        # Return
        return gdz_gal[gd_rho], rho[gd_rho]

    def get_observed(self, theta=None, subtab=None):
        """Generate a Table of observed targets

        Optionally to an angular distance from field center

        Parameters
        ----------
        theta : Quantity or Angle, optional
          Angular radius
        subtab : Table, optional
          User may input a table for processing
          theta is ignored

        Returns
        -------
        obs_targ : Table
          Sub-table of targets that have been observed within theta (if given)
          and/or within subtab (if given)
        obs_dates : List
          List of observing dates [eventually might add to Table]
        indices : array
          Indices from the target table
        """
        if (self.targets is None) or (self.observing is None):
            raise ValueError('IgmGalaxyField: Need to fill the target and/or observing table first!')
        if subtab is None:
            if theta is None:
                subtab = self.targets
            else:
                # Trim on angular cut first
                targ_coord = SkyCoord(ra=self.targets['TARG_RA']*u.deg,
                    dec=self.targets['TARG_DEC']*u.deg)
                sep = self.coord.separation(targ_coord)
                gdsep = np.where(sep < theta)[0]
                if len(gdsep) == 0:
                    return None
                subtab = self.targets[gdsep]
        else:
            gdsep = np.arange(len(subtab)) # For indexing below
        # Generate mask (set all to False; True is masked in numpy)
        tmsk = np.array([False]*len(subtab))
        # Grab those with a MASK_NAME
        have_mask = np.where(~subtab['MASK_NAME'].mask)[0]
        if len(have_mask) == 0:
            warnings.warn("No sources with a MASK_NAME")
            pdb.set_trace()
            return None
        # Get unique mask values
        all_masks = subtab['MASK_NAME'][have_mask]
        uni_masks = np.unique(all_masks)
        obs_dict = {}
        # Loop on these
        for mask in uni_masks:
            obs_dates = self.get_mask_obsdate(mask)
            if len(obs_dates) > 0:
                mt2 = np.where(subtab['MASK_NAME'][have_mask]==mask)
                tmsk[have_mask[mt2]] = True
                obs_dict[mask] = obs_dates
        # Finish
        return subtab[tmsk], obs_dict, gdsep[tmsk]

    def get_unobserved(self, theta=None):
        """Generate a Table of unobserved targets within an angular distance

        Parameters
        ----------
        theta : Quantity, optional
          Angular distance

        Returns
        -------
        unobs_targ : Table
          Sub-table of targets that have been not been observed within theta (if given)
        """
        if (self.targets is None) or (self.observing is None):
            raise ValueError('IgmGalaxyField: Need to fill the target and/or observing table first!')
        # Trim on angular cut first
        if theta is None:
            targ_coord = SkyCoord(ra=self.targets['TARG_RA']*u.deg,
                dec=self.targets['TARG_DEC']*u.deg)
            sep = self.coord.separation(targ_coord)
            gdsep = np.where(sep < theta)[0]
            if len(gdsep) == 0:
                return None
            # Set all to False to start
            subtab = self.targets[gdsep]
        else:
            subtab = self.targets
        tmsk = np.array([True]*len(subtab))
        # Grab observed (short cut!)
        obs_tab, odict, _ = self.get_observed(theta, subtab=subtab)
        # Remove those
        for kk, row in enumerate(subtab):
            if row['TARG_RA'] in obs_tab['TARG_RA']: # Could use DEC too
                tmsk[kk] = False
        # Return
        return subtab[tmsk]


    def get_mask_obsdate(self, mask_name):
        """Given a mask name, find the observing dates

        Parameters
        ----------
        mask_name : str
          Name of the mask

        Returns
        -------
        obs_dates : List
          List of the observing dates (can be empty)
        """
        if self.observing is None:
            raise ValueError('Need to fill observing info!')
        #
        mt = np.where(self.observing['MASK_NAME'] == mask_name)[0]
        if self.observing['DATE_OBS'].mask[mt[0]]:
            return []
        obs_dates = [self.observing['DATE_OBS'][imt] for imt in mt]
        # Return
        return obs_dates


    def clean_duplicates(self, table, tol=1*u.arcsec, method='first'):
        """ Clean duplicates in table based on (ra,dec) coordinates

        Parameters
        ----------
        table : Table
            Table to clean duplicates based on (ra, dec)
        tol : Angle, optional
            Angular tolerance for considering duplicates
        method : str, optional
            Method to use. Current options are:
            ``'first'``: if duplicates exist keep only the first one

        Returns
        -------
        cleaned_table : Table
            A version of `table` without duplicates
        """

        # TODO: add more methods for merging/cleaning duplicates

        if method not in ['first']:
            raise RuntimeError('Not ready for this method=`{}`'.format(method))

        isdup, idx, dcoord = pfu.check_dup_table(table, tol=tol)
        dup_inds = np.where(isdup == True)[0]
        keep = []

        if method == 'first':
            for ii in dup_inds:
                mtch = np.where(dcoord[ii].separation(dcoord) < tol)[0]
                keep.append(min(mtch))

            first_dup = np.unique(np.array(keep))

            no_dup = np.arange(len(idx))[~isdup]
            clean_inds = np.append(no_dup, first_dup)
            clean_inds = np.sort(clean_inds)

        #return
        return table[clean_inds]

    #    
    def __repr__(self):
        return ('<{:s}: {:s} {:s} {:s}>'.format(
                self.__class__.__name__,
                 self.name,
                 self.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                 self.coord.dec.to_string(sep=':', pad=True, alwayssign=True)))




