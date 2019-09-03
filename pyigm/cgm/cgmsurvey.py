""" Classes for CGM Surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import warnings
from IPython import embed
import json, io

from pkg_resources import resource_filename

from collections import OrderedDict

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import  units as u

from linetools.lists.linelist import LineList
from linetools import utils as ltu

from pyigm.utils import lst_to_array
from pyigm.surveys.igmsurvey import GenericIGMSurvey
from pyigm.cgm.cgm import CGMAbsSys

try:
    basestring
except NameError:  # For Python 3
    basestring = str

class CGMAbsSurvey(object):
    """A CGM Survey class in absorption

    Attributes
    ----------
    survey : str, optional
      Survey name
    ref : str, optional
      Reference(s)
    """

    @classmethod
    def from_json(cls, jfile, **kwargs):
        """ Instantiate from a single JSON file
        Parameters
        ----------
        jfile : str
        kwargs

        Returns
        -------
        CGMAbsSurvey

        """
        llist = LineList('ISM')
        #
        slf = cls(**kwargs)
        slf.load_json(jfile, llist=llist, **kwargs)
        return slf

    @classmethod
    def from_tarball(cls, tfile, **kwargs):
        """ Load the COS-Halos survey from a tarball of JSON files
        Parameters
        ----------
        tfile : str

        Returns
        -------

        """
        warnings.warn("This is being deprecated..  One will eventually only load/write JSON files.", DeprecationWarning)

        slf = cls(**kwargs)

        # Dict
        llist = LineList('ISM')
        slf.load_tarball(tfile, llist=llist, **kwargs)

        # Return
        return slf

    @classmethod
    def from_cgmabssys(cls, cgmlist, **kwargs):
        """ Instantiate new survey from list of CgmAbsSys objects

        Parameters
        ----------
        cgmlist : list CgmAbsSys


        """
        if not isinstance(cgmlist,list):
            raise IOError("Input must be list of CGMAbsSys")
        elif not isinstance(cgmlist[0],CGMAbsSys):
            raise IOError("Input must be list of CGMAbsSys")

        slf = cls(**kwargs)
        slf.cgm_abs.extend(cgmlist)
        return slf

    @classmethod
    def load_B16(cls,select_method='rvir'):
        """ Load the Burchett+16 z<0.015 samples

        Parameters
        ----------
        select_method : str, optional
            Selection scheme to associate galaxies with absorbers.  By default,
            the sample of galaxies with smallest impact parameters relative
            to their Rvir is loaded.  Otherwise, the closest in proper distance.

        """
        if select_method=='rvir':
            b16_tarfile = resource_filename('pyigm', '/data/CGM/z0/B16_vir_sys.tar')
        else:
            b16_tarfile = resource_filename('pyigm', '/data/CGM/z0/B16_kpc_sys.tar')
        print('Loading Burchett+16 using {:s} selection method'.format(select_method))
        # Load
        b16 = CGMAbsSurvey.from_tarball(b16_tarfile, chk_lowz=True, chk_z=False, build_sys=True)
        return b16

    @classmethod
    def load_J15(cls):
        """ Load the Johnson+15 sample
        """
        j15_tarfile = resource_filename('pyigm', '/data/CGM/z0/J15_sys.tar')
        j15 = CGMAbsSurvey.from_tarball(j15_tarfile, chk_lowz=False, chk_z=False, build_sys=True)
        return j15

    def __init__(self, survey='', ref='', **kwargs):
        """
        Parameters
        ----------
        survey : str, optional
        ref : str, optional

        Returns
        -------

        """
        # Name of survey
        self.survey = survey
        self.ref = ref

        # Several ways to hold the data
        self._dict = OrderedDict()
        self._data = Table()
        self.cgm_abs = []
        self.mask = None

    @property
    def nsys(self):
        """ Number of systems
        Returns
        -------
        nsys : int
        """
        if self.mask is not None:
            return np.sum(self.mask)
        elif len(self._data) > 0:
            return len(self._data)
        elif len(self._dict) > 0:
            return len(self._dict)
        else:
            return len(self.cgm_abs)

    def add_ion_to_data(self, inp_ion):
        """ Add columns for ion info (column density)
        to the Table

        Parameters
        ----------
        inp_ion

        Returns
        -------

        """
        from linetools.analysis import absline as ltaa
        ion = inp_ion.replace(' ','')
        # Check for previous
        if 'flag_N_{:s}'.format(ion) in self._data.keys():
            print("Ion data is already in the _data table")
            return
        # Loop on the systems
        flagNs, Ns, sigNs = [], [], []
        for key in self._dict.keys():
            igm_comp = self._dict[key]['igm_sys']['components']
            comps = [] # Allow for more than one
            for comp in igm_comp.keys():
                sion = comp.split('_')[0]
                if sion == ion:
                    if 'attrib' in igm_comp[comp].keys():
                        attrib = igm_comp[comp]['attrib']
                        attrib['sig_logN'] = np.array(attrib['sig_logN'])
                        comps.append(attrib.copy())
                    else: # Deprecated
                        comps.append(dict(logN=igm_comp[comp]['logN'],
                                          flag_N=igm_comp[comp]['flag_N'],
                                          sig_logN=np.array([igm_comp[comp]['sig_logN']]*2)))
            # Now sum em up
            if len(comps) == 0:
                flagNs.append(0)
                Ns.append(0.)
                sigNs.append(np.array([0.]*2))
                continue
            obj = dict(flag_N=comps[0]['flag_N'], logN=comps[0]['logN'], sig_logN=comps[0]['sig_logN'])
            for comp in comps[1:]:
                if comp['flag_N'] != 0:
                    obj['flag_N'], obj['logN'], obj['sig_logN'] = ltaa.sum_logN(obj, comp)
            # Save
            flagNs.append(obj['flag_N'])
            Ns.append(obj['logN'])
            sigNs.append(obj['sig_logN'])
        # Add to Table
        self._data.add_column(Column(flagNs, name='flag_N_{:s}'.format(ion)))
        self._data.add_column(Column(Ns, name='logN_{:s}'.format(ion)))
        self._data.add_column(Column(sigNs, name='sig_logN_{:s}'.format(ion)))

    def build_sys_from_dict(self, sys_name, llist=None, **kwargs):
        """ Build CGMAbsSys from the internal dict

        Parameters
        ----------
        sys_name : str
          Name of the system
        llist : LineList, optional
        kwargs : passed to CGMAbsSys.from_dict

        Returns
        -------
        CGMAbsSys

        """
        tdict = self._dict[sys_name]
        cgmsys = CGMAbsSys.from_dict(tdict, chk_vel=False, chk_sep=False, chk_data=False,
                                 use_coord=True, use_angrho=True,
                                     linelist=llist, **kwargs)
        # Return
        return cgmsys

    def build_systems_from_dict(self, **kwargs):
        """ Build all of the CGMAbsSys objects from the interal _dict

        Parameters
        ----------
        kwargs : passed to build_sys_from_dict

        Returns
        -------

        """
        # Check for LineList
        if 'llist' not in kwargs.keys():
            warnings.warn("It will likely speed things up greatly to pass a LineList in to build the systems")
        #
        for key in self._dict.keys():
            cgmsys = self.build_sys_from_dict(key, **kwargs)
            self.cgm_abs.append(cgmsys)
        return

    def build_data_from_dict(self):
        """ Build the _data Table from the internal _dict
        Columns are:  RA_IGM, DEC_IGM, RA, DEC, zabs, and likely a few others (e.g. rho)

        Other methods add on Ion info, etc.
        """
        # Table columns
        key0 = list(self._dict.keys())[0]
        tab_clms = list(self._dict[key0].keys())
        # Handle RA/DEC
        IGM_RA = []
        IGM_DEC = []
        zabs = []
        for key in self._dict.keys():
            IGM_RA.append(self._dict[key]['igm_sys']['RA'])
            IGM_DEC.append(self._dict[key]['igm_sys']['DEC'])
            zabs.append(self._dict[key]['igm_sys']['zabs'])
        self._data.add_column(Column(IGM_RA, name='RA_IGM'))
        self._data.add_column(Column(IGM_DEC, name='DEC_IGM'))
        self._data.add_column(Column(zabs, name='zabs'))
        # Remove unwanted ones
        rmv_keys = ['CreationDate', 'cosmo', 'ebv', 'user', 'Refs', 'igm_sys', 'galaxy']
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

    def to_json(self, outfile, overwrite=True):
        """ Generates a JSON file of the survey

        Parameters
        ----------
        outfil : str

        """
        survey_dict = OrderedDict()
        # Loop on systems
        for cgm_abs in self.cgm_abs:
            # Dict from copy
            cdict = cgm_abs.to_dict()
            # Use galaxy name for key;  Should be unique
            survey_dict[cgm_abs.galaxy.name+'_'+cgm_abs.igm_sys.name] = cdict.copy()

        # JSON
        clean_dict = ltu.jsonify(survey_dict)
        ltu.savejson(outfile, clean_dict, overwrite=overwrite)
        print("Wrote: {:s}".format(outfile))
        print("You may now wish to compress it..")

    def to_json_tarball(self, outfil):
        """ Generates a gzipped tarball of JSON files, one per system

        Parameters
        ----------
        outfil : str

        """
        import subprocess
        warnings.warn("This method is likely to be Deprecated", DeprecationWarning)
        print("Continue at your own peril...")
        embed(header='330')

        tmpdir = 'CGM_JSON'
        try:
            os.mkdir(tmpdir)
        except OSError:
            pass
        jfiles = []

        # Loop on systems
        for cgm_abs in self.cgm_abs:
            # Dict from copy
            cabscopy = cgm_abs.copy()
            cdict = cabscopy.to_dict()

            # Temporary JSON file
            json_fil = tmpdir+'/'+cabscopy.name+'.json'
            jfiles.append(json_fil)
            with io.open(json_fil, 'w', encoding='utf-8') as f:
                #try:
                f.write(json.dumps(cdict, sort_keys=True, indent=4,
                                           separators=(',', ': ')))
        # Tar
        subprocess.call(['tar', '-czf', outfil, tmpdir])
        print('Wrote: {:s}'.format(outfil))

        # Clean up
        for jfile in jfiles:
            os.remove(jfile)
        os.rmdir(tmpdir)

    def load_json(self, jfile, build_data=True, build_sys=False, verbose=True, **kwargs):
        """
        Parameters
        ----------
        jfile : str
        build_data : bool, optional
          Generate the internal _data Table  [very fast and recommended]
        build_sys : bool, optional
          Generate the list of cgm_abs objects from the internal _dict [May be slow]
        kwargs

        Returns
        -------

        """
        # Load
        self._dict = ltu.loadjson(jfile)
        # Generate
        if build_sys:
            self.build_systems_from_dict(**kwargs)

        # Galaxy coords
        ras = [self._dict[key]['RA'] for key in self._dict.keys()]
        decs = [self._dict[key]['DEC'] for key in self._dict.keys()]
        self.coords = SkyCoord(ra=ras, dec=decs, unit='deg')

        # Sightline coords
        ras = [self._dict[key]['igm_sys']['RA'] for key in self._dict.keys()]
        decs = [self._dict[key]['igm_sys']['DEC'] for key in self._dict.keys()]
        self.scoords = SkyCoord(ra=ras, dec=decs, unit='deg')

        # Data table
        if build_data:
            self.build_data_from_dict()
        # Return
        return

    def load_tarball(self, tfile, build_data=True, build_sys=False, llist=None, verbose=True, **kwargs):
        """
        Parameters
        ----------
        tfile
        build_data
        build_sys
        llist
        kwargs

        Returns
        -------

        """
        import tarfile
        # Load
        tar = tarfile.open(tfile)
        for kk, member in enumerate(tar.getmembers()):
            if '.json' not in member.name:
                print('Skipping a likely folder: {:s}'.format(member.name))
                continue
            # Extract
            f = tar.extractfile(member)
            try:
                tdict = json.load(f)
            except:
                if verbose:
                    print('Unable to load {}'.format(member))
                continue
            # Build dict
            self._dict[tdict['Name']] = tdict
            # Generate
            if build_sys:
                cgmsys = CGMAbsSys.from_dict(tdict, chk_vel=False, chk_sep=False, chk_data=False,
                                         use_coord=True, use_angrho=True,
                                         linelist=llist, **kwargs)
                self.cgm_abs.append(cgmsys)
        tar.close()

        # Galaxy coords
        ras = [self._dict[key]['RA'] for key in self._dict.keys()]
        decs = [self._dict[key]['DEC'] for key in self._dict.keys()]
        self.coords = SkyCoord(ra=ras, dec=decs, unit='deg')

        # Sightline coords
        ras = [self._dict[key]['igm_sys']['RA'] for key in self._dict.keys()]
        decs = [self._dict[key]['igm_sys']['DEC'] for key in self._dict.keys()]
        self.scoords = SkyCoord(ra=ras, dec=decs, unit='deg')

        # Data table
        if build_data:
            self.build_data_from_dict()
        # Return
        return


    def ion_tbl(self, Zion, fill_ion=True, vrange=None, extra_galaxy_attr=[], **kwargs):
        """ Generate a Table of Ionic column densities for an input ion

        Parameters
        ----------
        Zion : tuple or str
        fill_ion : bool, optional
          Fill each ionN table in the survey (a bit slow)
        vrange : Quantity, optional
          Velocity range of components to sum column densities

        Returns
        -------
        tbl : astropy.Table
           Returns None if there are no matches to input Zion
        """
        from linetools.abund.ions import name_to_ion
        if isinstance(Zion, basestring):
            Zion = name_to_ion(Zion)

        # Generate dummy IGMSurvey
        dumb = GenericIGMSurvey()
        names = []
        rhos = []
        gal_dict = dict(gal_ra=[], gal_dec=[])
        for gattr in extra_galaxy_attr:
            gal_dict[gattr] = []

        for cgmabs in self.cgm_abs:
            if fill_ion:
                cgmabs.igm_sys.fill_ionN(vrange=vrange,summed_ion=True)
            if cgmabs.igm_sys._ionN is not None:
                dumb._abs_sys.append(cgmabs.igm_sys)
                # Names
                names.append(cgmabs.name)
                # Impact parameters
                rhos.append(cgmabs.rho.to(u.kpc).value)
                # Galaxy info
                gal_dict['gal_ra'].append(cgmabs.galaxy.coord.ra.value)
                gal_dict['gal_dec'].append(cgmabs.galaxy.coord.dec.value)
                for gattr in extra_galaxy_attr:
                    gal_dict[gattr].append(getattr(cgmabs.galaxy, gattr))
        # Run ions
        tbl = dumb.ions(Zion,skip_null=False,**kwargs)
        if tbl is None:
            return None
        tbl.add_column(Column(names, name='cgm_name'))
        # Add impact parameter
        tbl.add_column(Column(rhos*u.kpc, name='rho_impact'))
        # Galaxy info
        for key in gal_dict:
            tbl[key] = gal_dict[key]
        # Return
        return tbl

    def component_tbl(self, Zion):
        """ Generate a Table of line measurements for an input ion broken
        down by component.

        Parameters
        ----------
        Zion : tuple or str
            E.g., (8,6) or 'OVI'

        Returns
        -------
        tbl : astropy.Table
        """
        from linetools.abund.ions import name_to_ion
        from linetools.isgm import utils as ltiu
        from astropy.table import Table, vstack, Column

        if isinstance(Zion, basestring):
            Zion = name_to_ion(Zion)
        newcomptab = Table()
        rhos = []
        names = []
        for i,isys in enumerate(self.cgm_abs):
            isys.igm_sys.update_component_vel()
            comptab = ltiu.table_from_complist(isys.igm_sys._components)
            comptab = comptab[(comptab['Z']==Zion[0]) & (comptab['ion']==Zion[1])]
            if len(comptab)==0:
                continue
            newcomptab = vstack([newcomptab,comptab])
            rhos.extend([isys.rho.value] * len(comptab))
            names.extend([isys.name] * len(comptab))
        newcomptab.add_column(Column(names, name='cgm_name'))
        newcomptab.add_column(Column(rhos*u.kpc, name='rho_impact'))
        return newcomptab

    def trans_tbl(self, inp, fill_ion=True):
        """ Generate a Table of Data on a given transition, e.g. SiIII 1206

        Parameters
        ----------
        inp : str or Quantity
          str -- Name of the transition, e.g. 'CII 1334'
          Quantity -- Rest wavelength of the transition, e.g. 1334.53*u.AA to 0.01 precision

        Returns
        -------
        tbl : astropy.Table
        """
        # Generate dummy IGMSurvey
        dumb = GenericIGMSurvey()
        names = []
        for cgmabs in self.cgm_abs:
            dumb._abs_sys.append(cgmabs.igm_sys)
            # Names
            names.append(cgmabs.name)
        # Run ions
        tbl = dumb.trans(inp)
        # Add CGM name
        tbl.add_column(Column(names, name='cgm_name'))
        # Return
        return tbl

    def abs_kin(self, lbl):
        """  Create a Table of the Kinematic info

        Parameters
        ----------
        lbl : string
          Label for the Kinematics dict
        TODO:
          Add wrest!!
        """
        from astropy.table import Table

        # The folowing line will fail if the first system is not populated..
        keys = list(self.cgm_abs[0].igm_sys.kin[lbl].keys())
        clms = []
        for ii in range(len(keys)):
            clms.append([])

        cgm_names = []
        for isys, cgm_abs in enumerate(self.cgm_abs):
            # Populated?
            try:
                kdict = cgm_abs.igm_sys.kin[lbl]
            except KeyError:
                continue
            if len(kdict) == 0:
                continue
            #
            cgm_names.append(cgm_abs.name)
            for jj,key in enumerate(keys):
                try:  # Deal with Quantity
                    clms[jj].append(kdict[key].value)
                except AttributeError:
                    clms[jj].append(kdict[key])
                except KeyError:
                    embed(header='606')

        # Table me
        t = Table(clms, names=keys)
        t['name'] = cgm_names
        return t

    def get_cgmsys(self, cgmname, return_index=False):
        """Convenience method to return CGMAbsSys by name

        Parameters
        ----------
        cgmname : str
            Name of CGMAbsSys to return
        return_index : bool
            If True, return index into self.cgm_abs of match

        Returns
        -------
        cgmsys : CGMAbsSys
            CGMAbsSys with name matching 'cgmname'
        index : int, optional
            Index into self.cgm_abs of match
        """
        names = np.array([str(system.name) for system in self.cgm_abs])
        index = np.where(names==cgmname)[0]
        if len(index) == 0:
            raise IOError("No CGMAbsSys with a matching name!")
        else:
            index=index[0]
        if return_index is True:
            return self.cgm_abs[index],index
        else:
            return self.cgm_abs[index]

    def __getattr__(self, k):
        # Try Self first
        try:
            lst = [getattr(cgm_abs, k) for cgm_abs in self.cgm_abs]
        except AttributeError:
            # Galaxy?
            try:
                lst = [getattr(cgm_abs.galaxy, k) for cgm_abs in self.cgm_abs]
            except AttributeError:
                # Try AbsLine_Sys last
                try:
                    lst = [getattr(cgm_abs.igm_sys, k) for cgm_abs in self.cgm_abs]
                except AttributeError:
                    print('cgm.core: Attribute not found!')
                    embed(header='655')
        # Special cases
        if k == 'coord':
            ra = [coord.icrs.ra.value for coord in lst]
            dec = [coord.icrs.dec.value for coord in lst]
            lst = SkyCoord(ra=ra, dec=dec, unit='deg')
            if self.mask is not None:
                return lst[self.mask]
            else:
                return lst
        elif k == 'scoord':  # Sightline coordinates
            lst = [getattr(cgm_abs.igm_sys, 'coord') for cgm_abs in self.cgm_abs]
            ra = [coord.icrs.ra.value for coord in lst]
            dec = [coord.icrs.dec.value for coord in lst]
            lst = SkyCoord(ra=ra, dec=dec, unit='deg')
            if self.mask is not None:
                return lst[self.mask]
            else:
                return lst
        # Return array
        return lst_to_array(lst, mask=self.mask)

    def __repr__(self):
        str1 = '<CGM_Survey: {:s} nsys={:d}, ref={:s}>\n'.format(self.survey, self.nsys, self.ref)
        return str1

    def __add__(self, cgmsurvey2, only_overlapping=True):

        if len(cgmsurvey2._dict) == 0:
            warnings.warn("You're second survey doesn't have a '._dict'. You need this to combine surveys.")
            warnings.warn("exiting")
            assert len(cgmsurvey2._dict) != 0

        if only_overlapping:
            print('Careful. The order of the objects matters here. Acts like a "LEFT JOIN"')
            # We only want to keep the cgm objects in the second survey that are in the fields of the first survey
            c = cgmsurvey2.coords
            catalog = self.coords

            # match on sky position
            idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, 15 * u.arcmin)

            # get the names of the matched systems
            names = np.unique(cgmsurvey2._data[idxc]['Name'].data)

            # make a new dict where we only keep the matches
            newhalos_dict = OrderedDict()
            for kk in cgmsurvey2._dict.keys():
                if cgmsurvey2._dict[kk]['Name'] in names:
                    newhalos_dict[kk] = cgmsurvey2._dict[kk]

            matched_survey = CGMAbsSurvey()
            matched_survey._dict = newhalos_dict

            # inherit
            matched_survey.survey = cgmsurvey2.survey
            matched_survey.ref = cgmsurvey2.ref

            # rename so the rest follows the same conventions
            cgmsurvey2 = matched_survey

        # add a label so you know which dict the data is coming from
        for ii, key in enumerate(cgmsurvey2._dict.keys()):
            cgmsurvey2._dict[key]['survey'] = cgmsurvey2.survey

        # add a label so you know which dict the data is coming from
        for ii, key in enumerate(self._dict.keys()):
            if self.survey is not None:
                self._dict[key]['survey'] = self.survey
            else:
                warnings.warn('You really should name your survey')

        # Merge the dicts
        combined_dict = {**self._dict, **cgmsurvey2._dict}

        # Instantiate a new CGMSurvey and give it a name
        combined_survey = CGMAbsSurvey()
        combined_survey.survey = "{0} + {1}".format(self.survey, cgmsurvey2.survey)

        # fill the new survey with a dict
        combined_survey._dict = combined_dict

        # Generate the data table
        combined_survey.build_data_from_dict()

        # return the new object
        return combined_survey