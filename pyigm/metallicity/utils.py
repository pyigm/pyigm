""" Module for simple calculations with MCMC outputs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import warnings

import numpy as np
import os
import pdb


def calc_logNH(hdf_file, modl=None, sys_error=0., NH_interpol=None, init_interpol=False,
                test=False, ret_flags=['NH'], debug=False, min_NHI=None):
    """ Use the outputs from a MCMC chain output file to generate an  array of log10 N_H values
    Can also return NO assuming 8.69

    Parameters
    ----------
    hdf_file : str
    sys_error : float, optional
      Add a systematic but random (Gaussian) error to the values
    init_interpol : bool, optional
      Simply return the interpolator?
    NH_interpol : RegularGridInterpolator, optional
      Input to speed up
    test : bool, optional
    min_NHI : float, optional
      Minimum NHI of the grid.  If input, correct any values to this level
      and then correct back
    ret_flags : list, optional
      Return list
      'NH' = NH_values
      'NHI' = NHI_values
      'NO' = NO_values  -- Assumes 8.69 but that might not be quite right
      'NSi' = NSi_values  -- 7.51 from Asplund 2009 as used in most Cloudy grids

    Returns
    -------
    NH_values : ndarray
     log10(NH) values for the input chain
    NHI_values : ndarray
     log10(NHI) values for the input chain
    NO_values : ndarray
     log10(NO) values for the input chain

    """
    import h5py, pickle
    from scipy import interpolate

    # Read the model grid
    if (modl is None) and (NH_interpol is None):
        model_file = os.getenv('DROPBOX_DIR')+'/cosout/grid_minextended.pkl'
        fil=open(model_file)
        modl=pickle.load(fil)
        fil.close()

    # Setup
    if NH_interpol is None:
        mod_axistag=modl[0]
        mod_axisval=[]
        xhi = modl[2]['HI']
        nhigrid=modl[3]['N_HI']
        NHgrid = nhigrid - xhi

        # Interpolator
        #define the dimension of the problem
        ndim=len(mod_axistag)
        nmodels=1
        for tt in mod_axistag:
            nmodels=nmodels*(modl[1][tt]).size
            #append axis value in a list
            mod_axisval.append(modl[1][tt])

        # Generate an Interpolation Grid
        NH_interpol = interpolate.RegularGridInterpolator(mod_axisval,NHgrid,
                                                           method='linear',bounds_error=False,fill_value=-np.inf)
        if init_interpol:
            return NH_interpol

    # PDFs
    fh5 = h5py.File(hdf_file, 'r')
    pdfs = fh5['outputs']['pdfs'].value
    tags = fh5['outputs/tags'].value
    fh5.close()

    # Deal with NHI off the grid?
    off = None
    if min_NHI is not None:
        # Checks
        if min_NHI > 16:
            warnings.warn("Approaching the optically thick limit.  Be warned..")
        assert tags[0] == 'col'
        #
        minNHI_pdf = np.min(pdfs[:,0])  # ASSUME NHI is in 0, aka col
        if minNHI_pdf < min_NHI:
            warnings.warn("Some NHI values are off the grid.  Correcting..")
            off = min_NHI - minNHI_pdf
            pdfs[:,0] += off

    # Here we go!
    if debug:
        pdb.set_trace()
    NH_values = NH_interpol(pdfs)

    # Offset back?
    if off is not None:
        NH_values -= off

    # NO? -- Before systematic error
    if 'NO' in ret_flags:
        assert tags[2] == 'met'
        NO_values = NH_values + pdfs[:,2] - 12 + 8.69  # Might not be right
    if 'NSi' in ret_flags:
        assert tags[2] == 'met'
        NSi_values = NH_values + pdfs[:,2] - 12 + 7.51

    # Systematic error (recommended)
    if sys_error > 0.:
        rand = np.random.normal(size=NH_values.size) * sys_error
        NH_values += rand

    if test:
        pdb.set_trace()
    # Finish
    retval = []
    for option in ret_flags:
        if option == 'NH':
            retval.append(NH_values)
        if option == 'NHI':
            retval.append(pdfs[:,0].flatten())
        if option == 'NO':
            retval.append(NO_values)
        if option == 'NSi':
            retval.append(NSi_values)
    return retval



