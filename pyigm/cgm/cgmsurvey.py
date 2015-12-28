""" Classes for CGM Surveys
"""
import pdb

from pyigm.utils import lst_to_array

class CGMAbsSurvey(object):
    """A CGM Survey class in absorption

    Attributes
    ----------
    survey : str, optional
      Survey name
    ref : str, optional
      Reference(s)
    """
    # Initialize with a .dat file
    def __init__(self, survey='', ref=''):

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

    # Extend attributes
    def __getattr__(self, k):
        # Try Self first
        try:
            lst = [getattr(cgm_abs, k) for cgm_abs in self.cgm_abs]
        except AttributeError:
            # Try AbsLine_Sys next
            try:
                lst = [getattr(cgm_abs.igm_sys, k) for cgm_abs in self.cgm_abs]
            except AttributeError:
                # Galaxy?
                try:
                    lst = [getattr(cgm_abs.galaxy, k) for cgm_abs in self.cgm_abs]
                except AttributeError:
                    print('cgm.core: Attribute not found!')
                    pdb.set_trace()
        # Return array
        return lst_to_array(lst, mask=self.mask)

    # Kinematics (for convenience)
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

    # Printing
    def __repr__(self):
        str1 = '<CGM_Survey: {:s} nsys={:d}, ref={:s}>\n'.format(self.survey, self.nsys, self.ref)
        for ii in range(self.nsys):
            str1 = str1+self.cgm_abs[ii].igm_sys.__repr__()+'\n'
        return str1
