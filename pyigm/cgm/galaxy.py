"""  Module for CGM of the Milky Way
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob
import pdb
import warnings
import h5py
import json, yaml

from astropy.io import fits, ascii
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord


from linetools.spectra import io as lsio
from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.analysis import abskin as laak
from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm.abssystem import GenericAbsSystem

from pyigm.metallicity.pdf import MetallicityPDF, DensityPDF, GenericPDF
from pyigm.cgm.cgmsurvey import CGMAbsSurvey
from pyigm.field.galaxy import Galaxy
from .cgm import CGM, CGMAbsSys
from pyigm.abssys.igmsys import IGMSystem
import pyigm

c_kms = const.c.to('km/s').value

class GalaxyCGM(CGM):
    """Inherits CGM Abs Survey

    Parameters:
    -----------
    fits_path : str, optional
      Path to the FITS data files for COS-Halos
    cdir : str, optional
      Path to the COS-Halos Dropbox
    kin_init_file : str, optional
      Path to kinematics file
    """
    def __init__(self, load=True, **kwargs):
        # Generate the Milky Way
        milkyway = Galaxy((0.,0.), z=0.)
        CGM.__init__(self, milkyway)
        self.refs = ''
        # Hot gas
        if load:
            self.load_hotgas()

    def load_hotgas(self):
        """ Load data on hot gas (e.g. OVII, OVIII)
        Fang+15
        """
        from linetools.lists.linelist import LineList
        llist = LineList('EUV',use_ISM_table=False)
        ovii = AbsLine('OVII 21', linelist=llist)

        # Fang+15  Table 1  [OVII]
        self.fang15 = Table.read(pyigm.__path__[0]+'/data/CGM/Galaxy/fang15_table1.dat', format='cds')
        # Reference
        if len(self.refs) > 0:
            self.refs += ','
        self.refs += 'Fang+15'
        # Generate the systems
        # # (should check to see if low-ion ones exist already)
        for row in self.fang15:
            # Coordinates
            gc = SkyCoord(l=row['GLON']*u.degree, b=row['GLAT']*u.degree, frame='galactic')
            # Limits
            # OVII line
            aline = ovii.copy()
            aline.attrib['coord'] = gc
            if isinstance(row['b'], float):
                aline.attrib['b'] = row['b'] * u.km / u.s
                aline.attrib['flag_EW'] = 1
                aline.attrib['EW'] = row['EW1'] / 1e3 * u.AA
                aline.attrib['sig_EW'] = row['e_EW1'] / 1e3 * u.AA
                vlim = np.array([-1,1]) * (2 * row['b'] + 2 * row['E_b']) * u.km/u.s
                # N_OVII
                aline.attrib['flag_N'] = 1
                aline.attrib['logN'] = row['logNO']
                aline.attrib['sig_logN'] = [row['e_logNO'], row['E_logNO']]
                z = row['Vel']/c_kms
            else:
                vlim = np.array([-300,300]) * u.km/u.s
                aline.attrib['flag_EW'] = 3
                aline.attrib['flag_N'] = 0  # Might be able to set an upper limit
                aline.attrib['EW'] = row['EW1'] / 1e3 * u.AA
                aline.attrib['sig_EW'] = 99. * u.AA
                z=0.
            # OVII
            aline.setz(z)
            aline.limits.set(vlim)
            # Generate component and add
            comp = AbsComponent.from_abslines([aline])
            # Instantiate
            abssys = GenericAbsSystem(gc, z, vlim)
            abssys.add_component(comp, chk_sep=False)
            # Add to cgm_abs
            self.cgm_abs.append(abssys)


    def write(self, outfil='COS-Halos_sys.tar.gz'):
        """ Write the survey to a tarball of JSON files

        Parameters
        ----------
        outfil : str, optional
        """
        self.to_json_tarball(outfil)

    def __getitem__(self, inp):
        """Grab CgmAbs Class from the list

        Parameters:
        -----------
        ion: tuple
          tuple:  (field,gal_id)
          str: field_gal_id

        Returns:
        ----------
        cgm_abs
        """
        if isinstance(inp,int):
            return self.cgm_abs[inp]
        elif isinstance(inp,tuple):
            if not isinstance(inp[0], basestring):
                raise IOError("Bad input")
            if not isinstance(inp[1], basestring):
                raise IOError("Bad input")
            field = inp[0]
            galid = inp[1]
        elif isinstance(inp, basestring):
            # Split at the first _
            under = inp.find('_')
            field = inp[:under]
            galid = inp[under+1:]
        else:
            raise IOError("Bad input")
        # Generate lists
        fields = np.array([cgm_abs.galaxy.field for cgm_abs in self.cgm_abs])
        galids = np.array([cgm_abs.galaxy.gal_id for cgm_abs in self.cgm_abs])
        #
        mt = np.where( (fields == field) & (galids == galid))[0]
        if len(mt) == 0:
            warnings.warn('CosHalos: CGM not found')
            return None
        elif len(mt) == 1:
            return self.cgm_abs[mt[0]]
        else:
            raise ValueError("Multiple hits.  Should not happen")


