""" Methods specific to DLAs, e.g. l(X)
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import imp
import pdb

from scipy.signal import medfilt
from scipy import interpolate

from astropy import cosmology
from astropy.table import Table

# Path for pyigm
pyigm_path = imp.find_module('pyigm')[1]


def lX(z):
    """ Returns l(X) from a preferred data source
    Currently Sanchez-Ramirez+16 with smoothing and interpolation

    Parameters
    ----------
    z : float or ndarray
      Redshift

    Returns
    -------
    lX : float or ndarray
      Float if z is a float
      Returns 0 if beyond the interpolated range
    """
    # Input
    if isinstance(z, float):
        z = np.array([z])
        flg_float = True
    else:
        flg_float = False

    # Sanchez-Ramirez et al. 2016 (Vanilla cosmology:  0.3, 0.7, 70)
    tab7_fil = pyigm_path+'/data/DLA/XQ-100/sramirez16_table7.dat'
    sr16_tab7 = Table.read(tab7_fil, format='ascii')

    # Smooth
    lxmed = medfilt(sr16_tab7['lx'], 21)

    # Interpolate
    flX = interpolate.interp1d(sr16_tab7['z'], lxmed, bounds_error=False, fill_value=0.)

    # Return
    eval = flX(z)
    if flg_float:
        return eval[0]
    else:
        return eval




