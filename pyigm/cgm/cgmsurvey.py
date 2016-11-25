""" Classes for CGM Surveys
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import warnings
import pdb
import json, io

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord

from pyigm.utils import lst_to_array
from pyigm.surveys.igmsurvey import GenericIGMSurvey
from pyigm.cgm.cgm import CGMAbsSys

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
    def from_tarball(cls, tfile, debug=False, **kwargs):
        """ Load the COS-Halos survey from a tarball of JSON files
        Parameters
        ----------
        tfile : str

        Returns
        -------

        """
        import tarfile
        import json
        from linetools.lists.linelist import LineList
        llist = LineList('ISM')

        slf = cls(**kwargs)
        # Load
        tar = tarfile.open(tfile)
        for kk, member in enumerate(tar.getmembers()):
            if '.' not in member.name:
                print('Skipping a likely folder: {:s}'.format(member.name))
                continue
            # Debug
            if debug and (kk == 5):
                break
            # Extract
            f = tar.extractfile(member)
            tdict = json.load(f)
            # Generate
            cgmsys = CGMAbsSys.from_dict(tdict, chk_vel=False, chk_sep=False, chk_data=False,
                                         use_coord=True, use_angrho=True,
                                         linelist=llist, **kwargs)
            slf.cgm_abs.append(cgmsys)
        tar.close()
        # Return
        return slf

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

        self.cgm_abs = []
        self.mask = None

    @property
    def nsys(self):
        """ Number of systems
        Returns
        -------
        nsys : int
        """
        return len(self.cgm_abs)

    def to_json_tarball(self, outfil):
        """ Generates a gzipped tarball of JSON files, one per system

        Parameters
        ----------
        outfil : str

        """
        import subprocess
        tmpdir = 'CGM_JSON'
        try:
            os.mkdir(tmpdir)
        except OSError:
            pass
        jfiles = []

        # Loop on systems
        for cgm_abs in self.cgm_abs:
            # Dict
            cdict = cgm_abs.to_dict()
            # Temporary JSON file
            json_fil = tmpdir+'/'+cgm_abs.name+'.json'
            jfiles.append(json_fil)
            with io.open(json_fil, 'w', encoding='utf-8') as f:
                #try:
                f.write(unicode(json.dumps(cdict, sort_keys=True, indent=4,
                                           separators=(',', ': '))))
        # Tar
        warnings.warn("Modify to write directly to tar file")
        subprocess.call(['tar', '-czf', outfil, tmpdir])
        print('Wrote: {:s}'.format(outfil))

        # Clean up
        for jfile in jfiles:
            os.remove(jfile)
        os.rmdir(tmpdir)

    def ion_tbl(self, Zion, fill_ion=True):
        """ Generate a Table of Ionic column densities for an input ion

        Parameters
        ----------
        Zion : tuple
        fill_ion : bool, optional
          Fill each ionN table in the survey (a bit slow)

        Returns
        -------
        tbl : astropy.Table
        """
        # Generate dummy IGMSurvey
        dumb = GenericIGMSurvey()
        names = []
        for cgmabs in self.cgm_abs:
            if fill_ion:
                cgmabs.igm_sys.fill_ionN()
            dumb._abs_sys.append(cgmabs.igm_sys)
            # Names
            names.append(cgmabs.name)
        # Run ions
        tbl = dumb.ions(Zion)
        # Add CGM name
        tbl.add_column(Column(names, name='cgm_name'))
        # Return
        return tbl

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

        keys = self.cgm_abs[0].igm_sys.kin[lbl].keys
        t = Table(names=keys,
                  dtype=self.cgm_abs[0].igm_sys.kin[lbl].key_dtype)

        for cgm_abs in self.cgm_abs:
            try:
                kdict = cgm_abs.igm_sys.kin[lbl]
            except KeyError:
                # No dict.  Filling in zeros
                row =  [0 for key in keys]
                t.add_row( row )
                continue
            # Filling
            row = [kdict[key] for key in keys]
            t.add_row( row )
        return t

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
                    pdb.set_trace()
        # Special cases
        if k == 'coord':
            ra = [coord.ra.value for coord in lst]
            dec = [coord.dec.value for coord in lst]
            lst = SkyCoord(ra=ra, dec=dec, unit='deg')
            return lst[self.mask]
        # Return array
        return lst_to_array(lst, mask=self.mask)

    def __repr__(self):
        str1 = '<CGM_Survey: {:s} nsys={:d}, ref={:s}>\n'.format(self.survey, self.nsys, self.ref)
        for ii in range(self.nsys):
            str1 = str1+self.cgm_abs[ii].igm_sys.__repr__()+'\n'
        return str1
