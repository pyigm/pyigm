""" Module for simple calculations with MCMC outputs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import warnings

import numpy as np
import os
import pdb


def calc_logNH(hdf_file, modl=None, sys_error=0., NH_interpol=None, init_interpol=False,
                test=False, ret_flags=['NH'], debug=False):
    """ Use the outputs from a MCMC chain output file to generate an  array of log10 N_H values
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
    ret_flag : int, optional
      Return bitwise flag
      1 = NH_values
      2 = NHI_values

    Returns
    -------
    NH_values : ndarray
     log10(NH) values for the input chain

    """
    import h5py, pickle
    from scipy import interpolate

    # Read the model grid
    if modl is None:
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
    fh5.close()

    # Here we go!
    if debug:
        pdb.set_trace()
    NH_values = NH_interpol(pdfs)

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
    return retval



