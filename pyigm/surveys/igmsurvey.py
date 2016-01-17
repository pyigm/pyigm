""" Class for IGMSurvey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import json
from abc import ABCMeta
import warnings
import pdb

from astropy import coordinates as coords
from astropy.io import ascii
from astropy import units as u
from astropy.table import QTable, Column, Table, vstack
from astropy.units.quantity import Quantity
from astropy.coordinates import SkyCoord

from linetools.spectra import io as lsio
from linetools.isgm import utils as ltiu

from pyigm.abssys.igmsys import IGMSystem
from pyigm.abssys.utils import class_by_type


class IGMSurvey(object):
    """
    Class for a survey of absorption line systems.

    Attributes
    ----------
    abs_type : str, unicode
      Type of Absorption system (DLA, LLS)
    ref : str, optional
      Reference(s) to the Survey
    _abs_sys : list
      List of AbsSystem objects
    mask : bool array, optional
      Defines a subset of the systems (e.g. statistical)
    sightlines : Table, optional
      Table of the sightlines in the survey
    """

    __metaclass__ = ABCMeta

    @classmethod
    def from_flist(cls, flist, tree=None, **kwargs):
        """ Read from list of .dat files (historical JXP format)

        Parameters
        ----------
        flist : str
          ASCII file including list of .dat files
        tree : str, optional
          Path to .dat files
        kwargs :
          Passed to __init__
        """
        if tree is None:
            tree = ''
        # Load up (if possible)
        data = ascii.read(tree+flist, data_start=0,
                          guess=False, format='no_header')
        slf = cls(**kwargs)
        slf.tree = tree
        slf.flist = flist

        # Load up
        slf.dat_files = list(data['col1'])
        # Generate IGMSys list
        for dat_file in slf.dat_files:
            slf._abs_sys.append(class_by_type(slf.abs_type).from_datfile(dat_file, tree=slf.tree))
        print('Read {:d} files from {:s} in the tree {:s}'.format(
            slf.nsys, slf.flist, slf.tree))

        return slf

    @classmethod
    def from_sfits(cls, summ_fits, **kwargs):
        """Generate the Survey from a summary FITS file

        Handles SPEC_FILES too.

        Parameters
        ----------
        summ_fits : str or Table or QTable
          Summary FITS file
        **kwargs : dict
          passed to __init__
        """
        # Init
        slf = cls(**kwargs)
        # Read
        if isinstance(summ_fits, Table):
            systems = summ_fits
        else:
            systems = QTable.read(summ_fits)
        nsys = len(systems)
        # Dict
        kdict = dict(NHI=['NHI', 'logNHI'],
                     sig_NHI=['sig(logNHI)', 'SIGNHI'],
                     name=['Name'], vlim=['vlim'],
                     zabs=['Z_LLS', 'ZABS'], zem=['Z_QSO', 'QSO_ZEM'],
                     RA=['RA'], Dec=['DEC', 'Dec'])
        # Parse the Table
        inputs = {}
        for key in kdict.keys():
            vals, tag = lsio.get_table_column(kdict[key], [systems],idx=0)
            if vals is not None:
                inputs[key] = vals
        # vlim
        if 'vlim' not in inputs.keys():
            default_vlim = [-1000, 1000.]* u.km / u.s
            inputs['vlim'] = [default_vlim]*nsys
        # Generate
        for kk in range(nsys):
            # Generate keywords
            kwargs = {}
            args = {}
            for key in inputs.keys():
                if key in ['vlim', 'zabs', 'RA', 'Dec']:
                    args[key] = inputs[key][kk]
                else:
                    kwargs[key] = inputs[key][kk]
            # Instantiate
            abssys = class_by_type(slf.abs_type)((args['RA'], args['Dec']), args['zabs'], args['vlim'], **kwargs)
            # spec_files
            try:
                abssys.spec_files += systems[kk]['SPEC_FILES'].tolist()
            except (KeyError, AttributeError):
                pass
            slf._abs_sys.append(abssys)
        # Mask
        slf.init_mask()
        # Return
        return slf

    def __init__(self, abs_type, ref=''):
        # Expecting a list of files describing the absorption systems
        """  Initiator

        Parameters
        ----------
        abs_type : str, unicode
          Type of IGMSystem in the Survey, e.g.  MgII, DLA, LLS
        ref : string, optional
          Reference(s) for the survey
        """
        self.abs_type = abs_type
        self.ref = ref
        self._abs_sys = []
        self.sightlines = None

        #

        # Mask
        self.mask = None
        self.init_mask()

        # Init
        self.flist = None

    @property
    def nsys(self):
        """ Number of systems
        Returns
        -------
        nsys : int
        """
        return len(self._abs_sys)

    def init_mask(self):
        """ Initialize the mask for abs_sys
        """
        if self.nsys > 0:
            self.mask = np.array([True]*self.nsys)

    def abs_sys(self):
        # Recast as an array
        return lst_to_array(self._abs_sys, mask=self.mask)

    def add_abs_sys(self, abs_sys):
        """ Add an IGMSys to the Survey

        Enables one to add checks

        Parameters
        ----------
        abs_sys : IGMSystem
        """
        assert self.chk_abs_sys(abs_sys)
        # Might check to see if a duplicate exists..

        # Append
        self._abs_sys.append(abs_sys)

    def calculate_gz(self, zstep=1e-4):
        """ Uses sightlines table to generate a g(z) array

        Parameters
        ----------
        zstep : float, optional
          Step size for g(z) array

        Returns
        -------
        zeval : ndarray
          Redshifts where g(z) is evaluate
        gz : ndarray
          g(z)
        """
        if self.sightlines is None:
            raise IOError("calculate_gz: Need to set sightlines table")
        # zeval
        zmin = np.min(self.sightlines['Z_START'])
        zmax = np.max(self.sightlines['Z_END'])
        zeval = np.arange(zmin, zmax, step=zstep)
        gz = np.zeros_like(zeval).astype(int)
        # Evaluate
        for row in self.sightlines:
            gd = (zeval >= row['Z_START']) & (zeval <= row['Z_END'])
            gz[gd] += 1
        # Return
        return zeval, gz


    def chk_abs_sys(self, abs_sys):
        """ Preform checks on input abs_sys

        Parameters
        ----------
        abs_sys : IGMSystem

        Returns
        -------
        bool

        """
        if not isinstance(abs_sys, IGMSystem):
            raise IOError("Must be an IGMSystem object")
        return True

    def fill_ions(self, use_Nfile=False, jfile=None, use_components=False):
        """ Loop on systems to fill in ions

        Parameters
        ----------
        jfile : str, optional
          JSON file containing the information
        use_Nfile : bool, optional
          Use (historic) .clm files?
        use_components : bool, optional
          Load up the Table with components (recommended)
        """
        if jfile is not None:
            # Load
            with open(jfile) as data_file:    
                ions_dict = json.load(data_file)
            # Loop on systems
            for abs_sys in self._abs_sys:
                abs_sys.get_ions(idict=ions_dict[abs_sys.name])
        elif use_Nfile:
            for abs_sys in self._abs_sys:
                abs_sys.get_ions(use_Nfile=True)
        elif use_components:
            for abs_sys in self._abs_sys:
                abs_sys._ionN = ltiu.iontable_from_components(abs_sys._components,
                                                              ztbl=abs_sys.zabs)
        else:
            raise ValueError("Not sure how to load the ions")

    # Get ions
    def ions(self, iZion, Ej=0., skip_null=False):
        """
        Generate a Table of columns and so on
        Restrict to those systems where flg_clm > 0

        Parameters
        ----------
        iZion : tuple
           Z, ion   e.g. (6,4) for CIV
        Ej : float [1/cm]
           Energy of the lower level (0. is resonance)
        skip_null : boolean (False)
           Skip systems without an entry, else pad with zeros 

        Returns
        -------
        Table of values for the Survey
        """
        if len(self.abs_sys()[0]._ionN) == 0:
            raise IOError("ionN table not set.  Use fill_ionN")
        #
        keys = [u'name', ] + self.abs_sys()[0]._ionN.keys()
        t = Table(self.abs_sys()[0]._ionN[0:1]).copy()   # Avoids mixin trouble
        t.add_column(Column(['dum'], name='name', dtype='<U32'))
        t = t[keys]
        if 'Ej' not in keys:
            warnings.warn("Ej not in your ionN table.  Ignoring. Be careful..")

        # Loop on systems (Masked)
        for abs_sys in self.abs_sys():
            # Grab
            if 'Ej' in keys:
                mt = ((abs_sys._ionN['Z'] == iZion[0])
                      & (abs_sys._ionN['ion'] == iZion[1])
                      & (abs_sys._ionN['Ej'] == Ej))
            else:
                mt = ((abs_sys._ionN['Z'] == iZion[0])
                      & (abs_sys._ionN['ion'] == iZion[1]))
            if np.sum(mt) == 1:
                irow = abs_sys._ionN[mt]
                # Cut on flg_clm
                if irow['flag_N'] > 0:
                    row = [abs_sys.name] + [irow[key] for key in keys[1:]]
                    t.add_row(row)   # This could be slow
                else:
                    if skip_null is False:
                        row = [abs_sys.name] + [0 for key in keys[1:]]
                        t.add_row(row)
            elif np.sum(mt) == 0:
                if skip_null is False:
                    row = [abs_sys.name] + [0 for key in keys[1:]]
                    t.add_row( row )
                continue
            else:
                raise ValueError("Multple entries")


        # Return
        return t[1:]

    # Mask
    def update_mask(self, mask, increment=False):
        """ Update the Mask for the abs_sys

        Parameters
        ----------
        mask : array (usually Boolean)
           Mask of systems 
        increment : bool, optional
           Increment the mask (i.e. keep False as False)
        """
        if len(mask) == len(self._abs_sys):  # Boolean mask
            if increment is False:
                self.mask = mask
            else:
                self.mask = self.mask & mask
        else:
            raise ValueError('abs_survey: Needs developing!')

    def __getattr__(self, k):
        """ Generate an array of attribute 'k' from the IGMSystems

        Mask is applied

        Parameters
        ----------
        k : str
          Attribute

        Returns
        -------
        numpy array
        """
        try:
            lst = [getattr(abs_sys, k) for abs_sys in self._abs_sys]
        except ValueError:
            raise ValueError("Attribute does not exist")
        # Special cases
        if k == 'coord':
            ra = [coord.ra for coord in lst]
            dec = [coord.dec for coord in lst]
            lst = SkyCoord(ra=ra, dec=dec)
            return lst[self.mask]
        # Recast as an array
        return lst_to_array(lst, mask=self.mask)

    def __add__(self, other, toler=2*u.arcsec):
        """ Combine one or more IGMSurvey objects

        Routine does a number of checks on the abstype,
        the uniqueness of the sightlines and systems, etc.

        Parameters
        ----------
        other : IGMSurvey
        toler : Angle or Quantity
          Tolerance for uniqueness

        Returns
        -------
        combined : IGMSurvey

        """
        # Check the Surveys are the same type
        if self.abs_type != other.abs_type:
            raise IOError("Combined surveys need to be same abs_type")

        # Init
        combined = IGMSurvey(self.abs_type)
        combined.ref = self.ref + ',' + other.ref

        # Check for unique systems
        other_coord =other.coord
        for abssys in self._abs_sys:
            if np.sum((abssys.coord.separation(other_coord) < toler) & (
                        np.abs(abssys.zabs-other.zabs) < (1000*(1+abssys.zabs)/3e5))) > 0:
                raise NotImplementedError("Need to deal with this")
        # Combine systems
        combined._abs_sys = self._abs_sys + other._abs_sys
        combined.mask = np.concatenate((self.mask, other.mask))

        # Sightlines?
        if self.sightlines is not None:
            slf_scoord = SkyCoord(ra=self.sightlines['RA']*u.deg,
                                  dec=self.sightlines['DEC']*u.deg)
            oth_scoord = SkyCoord(ra=other.sightlines['RA']*u.deg,
                                  dec=other.sightlines['DEC']*u.deg)
            idx, d2d, d3d = coords.match_coordinates_sky(slf_scoord,
                                                     oth_scoord, nthneighbor=1)
            mt = d2d < toler
            if np.sum(mt) > 0:
                raise NotImplementedError("Need to deal with this")
            else:
                # Combine systems
                combined.sightlines = vstack([self.sightlines,
                                              other.sightlines])
        # Return
        return combined

    def __repr__(self):
        if self.flist is not None:
            return '<IGMSurvey: {:s} {:s}, nsys={:d}, type={:s}, ref={:s}>'.format(
                self.tree, self.flist, self.nsys, self.abs_type, self.ref)
        else:
            repr = '<IGMSurvey: nsys={:d}, type={:s}, ref={:s}'.format(
                self.nsys, self.abs_type, self.ref)
            if self.sightlines is not None:
                repr = repr + ', nsightlines={:d}'.format(len(self.sightlines))
            repr = repr +'>'
            return repr


class GenericIGMSurvey(IGMSurvey):
    """A simple absorption line survey
    """
    def __init__(self, **kwargs):
        IGMSurvey.__init__(self, 'Generic', **kwargs)


def lst_to_array(lst, mask=None):
    """ Simple method to convert a list to an array

    Allows for a list of Quantity objects

    Parameters
    ----------
    lst : list
      Should be number or Quantities
    mask : boolean array, optional

    Returns
    -------
    array or Quantity array

    """
    if mask is None:
        mask = np.array([True]*len(lst))
    if isinstance(lst[0], Quantity):
        return Quantity(lst)[mask]
    else:
        return np.array(lst)[mask]

