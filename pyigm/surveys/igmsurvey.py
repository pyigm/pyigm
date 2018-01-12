""" Class for IGMSurvey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import json
from abc import ABCMeta
import warnings
import pdb

from collections import OrderedDict

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

try:
    unic = unicode
except:
    unic = str


class IGMSurvey(object):
    """ Class for a survey of absorption line systems.

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
    _data : Table, optional
      Table of 'key' data
    _dict : OrderedDict, optional
      Nested data
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
    def from_sfits(cls, summ_fits, coords=None, **kwargs):
        """Generate the Survey from a summary FITS file or Table

        Handles SPEC_FILES too.

        Parameters
        ----------
        summ_fits : str or Table or QTable
          Summary FITS file
        coords : SkyCoord array
          Contains all the coords for all the systems
        **kwargs : dict
          passed to __init__
        """
        # Init
        slf = cls(**kwargs)
        # Read
        if isinstance(summ_fits, Table):
            systems = summ_fits
        else:
            systems = Table.read(summ_fits)
        nsys = len(systems)
        # Special keys
        kdict = dict(NHI=['NHI', 'logNHI'],
                     sig_NHI=['sig(logNHI)', 'SIGNHI', 'NHI_ERR'],
                     name=['Name'], vlim=['vlim'],
                     zabs=['Z_LLS', 'ZABS', 'zabs'],
                     zem=['Z_QSO', 'QSO_ZEM', 'ZEM'],
                     RA=['RA'], DEC=['DEC', 'Dec'])
        # Parse the Table to make uniform the keys used
        for key in kdict.keys():
            for ikey in kdict[key]:
                if ikey in systems.keys():
                    if ikey == key:
                        pass
                    else:
                        systems.rename_column(ikey, key)
        # Set
        slf._data = systems
        # vlim
        if 'vlim' not in slf._data.keys():
            default_vlim = [-1000, 1000.]* u.km / u.s
            slf._data['vlim'] = [default_vlim]*nsys
        # Coords
        if coords is None:
            coords = SkyCoord(ra=slf._data['RA'], dec=slf._data['DEC'], unit='deg')
        slf.coords = coords
        # Mask
        slf.mask = None
        slf.init_mask()
        # Return
        return slf

    def __init__(self, abs_type, ref='', verbose=False):
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
        self._data = Table()
        self._dict = OrderedDict()
        self.sightlines = None
        self.coords = None  # Intended to be a SkyCoord obj with *all* of the system coordinates
        self.verbose=verbose

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
          Number of statistical if mask is set
        """
        if self.mask is not None:
            return np.sum(self.mask)
        elif len(self._data) > 0:
            return len(self._data)
        elif len(self._dict) > 0:
            return len(self._dict)
        else:
            return len(self._abs_sys)

    def sys_idx(self, abssys_name):
        """ System index"""
        # Index
        try:
            idx = list(self._dict.keys()).index(abssys_name)
        except ValueError:
            raise ValueError("System {:s} is not in the _dict".format(abssys_name))
        # Return
        return idx


    def abs_sys(self, inp, fill_coord=True):
        """ Return an abs_system by index
        Instantiate as needed
        Returns
        -------
        inp : int

        """
        # Mask?
        if self.mask is not None:
            idx = np.where(self.mask)[0][inp]
        else:
            idx = inp
        # Pull the system
        isys = self._abs_sys[idx]
        # Add coord
        if fill_coord:
            isys.coord = self.coords[idx]
        return isys

    def init_abs_sys(self, clobber=False):
        """ Initialize the abs_sys list
        """
        if (len(self._abs_sys) == 0) or clobber:
            self._abs_sys = [None]*self.nsys
        else:
            warnings.warn("abs_sys list is already initialized.  Use clobber=True to reset")

    def init_mask(self):
        """ Initialize the mask for abs_sys
        """
        if self.nsys > 0:
            self.mask = np.array([True]*self.nsys)

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

    def build_all_abs_sys(self, linelist=None, **kwargs):
        """
        Build all of the AbsSystem objects from the _dict (maybe _data soon)
        """
        # This speeds things up a bunch
        if linelist is None:
            linelist = default_linelist(self.verbose)
        # Loop me
        print("Starting the AbsSystem build.  Be patient..")
        for key in self._dict.keys():
            _ = self.build_abs_sys(key, linelist=linelist, **kwargs)
        print("Done!")

    def build_abs_sys(self, abssys_name, **kwargs):
        """ Build an AbsSystem from the _dict (and maybe _data
        soon enough)
        The item in self._abs_sys is filled and
        the systems is also returned

        Parameters
        ----------
        abssys_name : str
          Needs to match a key in the dict
        **kwargs
          Passed to compoenents_from_dict()

        Returns
        -------
        abs_sys : AbsSystem

        """
        # Index
        idx = self.sys_idx(abssys_name)
        # Components first
        comps = self.components_from_dict(abssys_name, coord=self.coords[idx], **kwargs)
        # Add HI if needed
        NHI = None
        if 'flag_NHI' in self._dict[abssys_name].keys():
            Zions = [comp.Zion for comp in comps]
            if (1,1) not in Zions:
                NHI = self._dict[abssys_name]['NHI']
        # Now the System
        vlim=np.array(self._dict[abssys_name]['vlim'])*u.km/u.s
        abssys = class_by_type(self.abs_type).from_components(comps, vlim=vlim, NHI=NHI)
        # Fill
        if len(self._abs_sys) == 0:
            self.init_abs_sys()
        self._abs_sys[idx] = abssys
        # Return too
        return abssys

    def calculate_gz(self, zstep=1e-4, zmin=None, zmax=None):
        """ Uses sightlines table to generate a g(z) array

        Parameters
        ----------
        zstep : float, optional
          Step size for g(z) array
        zmin : float, optional
          Minimum redshift of evaluated array.  Default is minimum in the sightlines
        zmax : float, optional
          Maximum redshift of evaluated array.  Default is maximum in the sightlines

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
        if zmin is None:
            zmin = np.min(self.sightlines['Z_START'])
        if zmax is None:
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

    def components_from_dict(self, abssys_name, coord=None, linelist=None):
        """ Build and return a list of AbsComponent objects
        from the dict for a given system

        Parameters
        ----------
        abssys_name : str
        coord : SkyCoord, optional
          coordinates to use for the components
        linelist : LineList, optional

        Returns
        -------
        compllist : list of AbsComponent objects

        """
        # Do it
        if linelist is None:
            linelist = default_linelist(self.verbose)
        # Components
        comps = ltiu.build_components_from_dict(self._dict[abssys_name],
                                                coord=coord, linelist=linelist)
        # Return
        return comps

    def data_from_dict(self):
        """ Generate the data Table from the internal dict
        """
        from astropy.table import Column
        # Table columns
        key0 = list(self._dict.keys())[0]
        tab_clms = list(self._dict[key0].keys())
        # Remove unwanted ones
        rmv_keys = ['abs_type', 'components', 'kin', 'Refs']
        for rkey in rmv_keys:
            if rkey in tab_clms:
                tab_clms.remove(rkey)
        # Build it
        for tclm in tab_clms:
            values = []
            for key in self._dict.keys():
                values.append(self._dict[key][tclm])
            # Add column
            clm = Column(values, name=tclm)
            self._data.add_column(clm)

    def fill_ions(self, use_Nfile=False, jfile=None, use_components=False,
                  verbose=True):
        """ Loop on systems to fill in _ionN Table

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
                abs_sys.get_ions(use_Nfile=True, verbose=verbose)
        elif use_components:
            for abs_sys in self._abs_sys:
                abs_sys._ionN = ltiu.table_from_complist(abs_sys._components)
        else:
            raise ValueError("Not sure how to load the ions")

    # Get ions
    def ions(self, Zion, Ej=0., skip_null=False):
        """ Generate a Table of columns and so on
        Restrict to those systems where flg_clm > 0

        Parameters
        ----------
        Zion : tuple
           Z, ion   e.g. (6,4) for CIV
        Ej : float [1/cm]
           Energy of the lower level (0. is resonance)
        skip_null : boolean (False)
           Skip systems without an entry, else pad with zeros 

        Returns
        -------
        tbl : MaskedTable of values for the Survey
          Systems without the ion have rows masked
        """
        if self._abs_sys[0]._ionN is None:
            raise IOError("ionN tables are not set.  Use fill_ionN")

        # Loop me!
        tbls = []
        names = []
        for kk,abs_sys in enumerate(self._abs_sys):
            if len(abs_sys._ionN) == 0:
                names.append('MASK_ME')
                tbls.append(None)
            #
            mt = (abs_sys._ionN['Z'] == Zion[0]) & (abs_sys._ionN['ion'] == Zion[1]) & (
                abs_sys._ionN['Ej'] == Ej)
            if np.any(mt):
                tbls.append(abs_sys._ionN[mt])
                names.append(abs_sys.name)
            else:
                tbls.append(None)
                names.append('MASK_ME')
        # Fill in the bad ones
        names = np.array(names)
        idx = np.where(names != 'MASK_ME')[0]
        if len(idx) == 0:
            return None
        bad = np.where(names == 'MASK_ME')[0]
        for ibad in bad:
            tbls[ibad] = tbls[idx[0]]
        # Stack me
        tbl = vstack(tbls)
        tbl['abssys_name'] = names
        # Mask
        tbl = Table(tbl, masked=True)
        mask = names == 'MASK_ME'
        for key in tbl.keys():
            if key == 'flag_N':
                tbl[key][mask] = 0
            else:
                tbl[key].mask = mask
        '''
        #
        keys = [u'abssys_name', ] + list(self._abs_sys[kk]._ionN.keys())
        t = Table(self._abs_sys[kk]._ionN[0:1]).copy()   # Avoids mixin trouble
        t.add_column(Column(['dum']*len(t), name='name', dtype='<U32'))
        t = t[keys]
        if 'Ej' not in keys:
            warnings.warn("Ej not in your ionN table.  Ignoring. Be careful..")

        # Loop on systems (Masked)
        for abs_sys in self._abs_sys:
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
                    row = [abs_sys.name] + [irow[key][0] for key in keys[1:]]
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
                pdb.set_trace()
                raise ValueError("Multple entries")
        '''
        # Reorder
        all_keys = list(tbl.keys())
        all_keys.remove('abssys_name')
        all_keys = ['abssys_name']+all_keys
        # Return
        return tbl[all_keys]

    def trans(self, inp):
        """ Generate a Table of Data on a given transition, e.g. SiIII 1206

        Parameters
        ----------
        inp : str or Quantity
          str -- Name of the transition, e.g. 'CII 1334'
          Quantity -- Rest wavelength of the transition, e.g. 1334.53*u.AA
            to 0.01 precision

        Returns
        -------
        tbl : astropy.Table
        """
        attrib = ['sys', 'z', 'flag_EW', 'EW', 'sig_EW', 'flag_N', 'logN', 'sig_logN']
        nattrib = len(attrib)
        clms = []
        for ii in range(nattrib):
            clms.append([])
        for abs_sys in self._abs_sys:
            # Name
            clms[0].append(abs_sys.name)
            #
            aline = abs_sys.get_absline(inp)
            if aline is None:
                for jj in range(1,nattrib):
                    clms[jj].append(0)
            else:
                for jj in range(1,nattrib):
                    try:  # Deal with Quantity
                        clms[jj].append(aline.attrib[attrib[jj]].value)
                    except AttributeError:
                        clms[jj].append(aline.attrib[attrib[jj]])
                    except KeyError:
                        clms[jj].append(0)
        # Generate the Table
        tbl = Table(clms, names=attrib)
        # Return
        return tbl

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

    def write_survey(self, outfile='tmp.tar', tmpdir = 'IGM_JSON', chk_dict=True):
        """ Generates a gzipped tarball of JSON files, one per system

        Parameters
        ----------
        outfile : str, optional
          Output filename
        tmpdir : str, optional
        chk_dict : bool, optional
          Check that the _dict matches what is in the _abs_sys list

        Returns
        -------

        """
        import os, io
        import subprocess
        try:
            os.mkdir(tmpdir)
        except OSError:
            pass
        jfiles = []

        # Loop on systems
        for igm_abs in self._abs_sys:
            # Dict
            idict = igm_abs.to_dict()
            # Temporary JSON file
            json_fil = tmpdir+'/'+igm_abs.name+'.json'
            jfiles.append(json_fil)
            with io.open(json_fil, 'w', encoding='utf-8') as f:
                f.write(unic(json.dumps(idict, sort_keys=True, indent=4,
                                           separators=(',', ': '))))
        # Tar
        subprocess.call(['tar', '-czf', outfile, tmpdir])
        print('Wrote: {:s}'.format(outfile))

        # Clean up
        for jfile in jfiles:
            try:
                os.remove(jfile)
            except OSError:  # Likely a duplicate.  This can happen
                pass
        os.rmdir(tmpdir)

    def __getattr__(self, k):
        """ Generate an array of attribute 'k' from the IGMSystems
        NOTE: We only get here if the Class doesn't have this attribute set already

        The Mask will be applied

        Order of search is:
          _data
          _dict
          _abs_sys

        Parameters
        ----------
        k : str
          Attribute

        Returns
        -------
        numpy array
        """
        # Special case(s)
        if k == 'coord':
            lst = self.coords
            if self.mask is not None:
                return lst[self.mask]
            else:
                return lst
        elif k in self._data.keys():  # _data
            lst = self._data[k]
        else:
            lst = None
        # Now try _dict
        if lst is None:
            if len(self._dict) > 0:
                if k in next(iter(self._dict.items()))[1].keys():
                    lst = [self._dict[key][k] for key in self._dict.keys()]
        # AbsSystem last!
        if lst is None:
            if len(self._abs_sys) == 0:
                raise ValueError("Attribute does not exist anywhere!")
            try:
                lst = [getattr(abs_sys, k) for abs_sys in self._abs_sys]
            except ValueError:
                raise ValueError("Attribute does not exist")
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
        if self.ref is not None:
            combined.ref = self.ref + ',' + other.ref
        else:
            combined.red = None

        # Check for unique systems
        other_coord =other.coord
        for abssys in self._abs_sys:
            if np.sum((abssys.coord.separation(other_coord) < toler) & (
                        np.abs(abssys.zabs-other.zabs) < (1000*(1+abssys.zabs)/3e5))) > 0:
                raise NotImplementedError("Need to deal with this")
        # Combine systems
        combined._abs_sys = self._abs_sys + other._abs_sys
        if self.mask is not None:
            combined.mask = np.concatenate((self.mask, other.mask)).flatten()
        else:
            combined.mask = None

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
        # Generate the Table
        tbl = Table(clms, names=attrib)
        # Return
        return tbl


def default_linelist(verbose=True):
    from linetools.lists.linelist import LineList
    if verbose:
        print("No LineList input.  Assuming you want the ISM list")
    linelist = LineList('ISM')
    return linelist
