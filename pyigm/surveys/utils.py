""" Utilities for IGMSurvey code
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import glob
import json
import pdb

from linetools import utils as ltu

from pyigm.abssys import utils as pyasu

from .igmsurvey import IGMSurvey


def load_sys_files(inp, type, ref=None, sys_path=False, **kwargs):
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

    Returns
    -------
    survey : IGMSurvey
    """
    import tarfile
    #
    survey = class_by_type(type)(ref=ref)
    system = pyasu.class_by_type(type)
    if sys_path:
        # Individual files
        files = glob.glob(inp+'*.json')
        files.sort()
        for ifile in files:
            tdict = ltu.loadjson(ifile)
            abssys = system.from_dict(tdict)
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
            # Generate
            abssys = system.from_dict(tdict, chk_sep=False, **kwargs)   # Consider use_coord=True as default
            survey._abs_sys.append(abssys)
        tar.close()
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

