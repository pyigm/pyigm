"""  Module for the QPQ survey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob
import pdb
import warnings
import h5py
import json, yaml

from pkg_resources import resource_filename

#from astropy.io import fits, ascii
#from astropy import units as u
from astropy.table import Table, Column

#from linetools.spectra import io as lsio
#from linetools.spectralline import AbsLine
#from linetools.analysis import absline as ltaa
#from linetools.analysis import abskin as laak
#from linetools.isgm.abscomponent import AbsComponent

#from pyigm.metallicity.pdf import MetallicityPDF, DensityPDF, GenericPDF
#from pyigm.cgm.cgmsurvey import CGMAbsSurvey
#from pyigm.field.galaxy import Galaxy
#from .cgm import CGMAbsSys
#from pyigm.abssys.igmsys import IGMSystem
#import pyigm

try:
    basestring
except NameError:  # For Python 3
    basestring = str


def load_qpq(v):
    """ Reads table with data
    Parameters
    ----------
    v : int
      version (6,7,8)

    Returns
    -------
    qpqdata : Table with data

    """

    if v == 8:
        q8file = resource_filename('pyigm', 'data/CGM/QPQ/qpq8_all_measured.dat')
        pdb.set_trace()
        qpqdata = Table.read(q8file, format="ascii") #, guess=False, header_start=None)

    if v == 7:
        q7file = resource_filename('pyigm', 'data/CGM/QPQ/qpq7_pairs.fits.gz')
        qpqdata = Table.read(q7file)

    if v == 6:
        q6file = resource_filename('pyigm', 'data/CGM/QPQ/qpq6_final_cut_2015nov16.fits')
        qpqdata = Table.read(q6file)

    return qpqdata


