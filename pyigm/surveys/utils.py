""" Utilities for IGMSurvey code
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import glob
import json
import pdb

from astropy.coordinates import SkyCoord

from linetools import utils as ltu
from linetools.lists.linelist import LineList

from pyigm.abssys import utils as pyasu

from .igmsurvey import IGMSurvey

# Load here to speed up line making
llist = LineList('ISM')


def load_sys_files(inp, type, ref=None, sys_path=False, build_abs_sys=False, **kwargs):
    """ Load up a set of SYS files from the hard-drive (JSON files)

    Parameters
    ----------
    inp : str
      Name of JSON tarball or if sys_path=True then the path to a folder of JSON files
    type : str
      type of IGMSystem, e.g. LLS
    ref : str, optional
      Reference label
    sys_path : str, optional
      indicates that inp is a path to a set of JSON SYS files
      otherwise, inp should be the filename of a tarball of JSON files
    build_abs_sys : bool, optional
      Build a list of AbsSystem's?  Can always be instantiated later
    **kwargs :
      Passed to system

    Returns
    -------
    survey : IGMSurvey
    """
    import tarfile
    #
    survey = class_by_type(type)(ref=ref)
    system = pyasu.class_by_type(type)
    if sys_path:
        pdb.set_trace()  # THIS NEEDS TO BE UPDATED AS WAS DONE FOR THE TARBALL
        # Individual files
        files = glob.glob(inp+'*.json')
        files.sort()
        for ifile in files:
            tdict = ltu.loadjson(ifile)
            abssys = system.from_dict(tdict, linelist=llist)
            survey._abs_sys.append(abssys)
    else:  # tarball
        print('Loading systems from {:s}'.format(inp))
        tar = tarfile.open(inp)
        for member in tar.getmembers():
            if '.' not in member.name:
                print('Skipping a likely folder: {:s}'.format(member.name))
                continue
            # Extract
            f = tar.extractfile(member)
            tdict = json.load(f)
            # Add keys (for backwards compatability)
            if ('NHI' in tdict.keys()) and ('flag_NHI' not in tdict.keys()):
                tdict['flag_NHI'] = 1
            # Add to list of dicts
            survey._dict[tdict['Name']] = tdict
        tar.close()

    # Set coordinates
    ras = [survey._dict[key]['RA'] for key in survey._dict.keys()]
    decs = [survey._dict[key]['DEC'] for key in survey._dict.keys()]
    survey.coords = SkyCoord(ra=ras, dec=decs, unit='deg')

    # Dummy abs_sys
    if build_abs_sys:
        survey.build_all_abs_sys()

    # Generate the data table
    print("Building the data Table from the internal dict")
    survey.data_from_dict()

    # Return
    return survey


def class_by_type(type):
    """ Enable IGMSurvey init by type

    Parameters
    ----------
    type : str

    Returns
    -------

    """
    from .llssurvey import LLSSurvey
    from .dlasurvey import DLASurvey

    if type == 'LLS':
        survey = LLSSurvey
    elif type == 'DLA':
        survey = DLASurvey
    else:
        raise IOError("Bad survey type!")
    # Return
    return survey

