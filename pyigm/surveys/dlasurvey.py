""" Class for DLA Surveys
"""
import numpy as np
import imp, glob
import pdb
import urllib2
import h5py


#from astropy.table import QTable, Column, Table
#from astropy import units as u
#from astropy.coordinates import SkyCoord


from pyigm.surveys.igmsurvey import IGMSurvey

pyigm_path = imp.find_module('pyigm')[1]


# Class for DLA Survey
class DLASurvey(IGMSurvey):
    """An DLA Survey class

    Attributes:

    """
    @classmethod
    def default_sample(cls):
        """
        Returns
        -------
        dlasurvey : IGMSurvey
        """
        # Default sample of DLA:  Neeleman
        if os.getenv('DLA') is None:
            print('Need to grab the DLA tree from JXP')
            return None
        dlasurvey = cls.from_flist('Lists/Neeleman13.lst', tree=os.environ.get('DLA'))
        dlasurvey.ref = 'Neeleman+13'

        # Return
        return dlasurvey

    def __init__(self, **kwargs):
        IGMSurvey.__init__(self, 'DLA', **kwargs)

