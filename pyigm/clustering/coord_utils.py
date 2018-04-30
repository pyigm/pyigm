"""Utilities for coordinate calculations
Taken from pyntejos on GitHub (via Crighton)
"""
import numpy as np


# constants
DEG_PER_HR = 360. / 24.             # degrees per hour
DEG_PER_MIN = DEG_PER_HR / 60.      # degrees per min
DEG_PER_S = DEG_PER_MIN / 60.       # degrees per sec
DEG_PER_AMIN = 1./60.               # degrees per arcmin
DEG_PER_ASEC = DEG_PER_AMIN / 60.   # degrees per arcsec
RAD_PER_DEG = np.pi / 180.             # radians per degree

def radec_to_xyz(ra_deg, dec_deg):
    """ Convert RA and Dec to xyz positions on a unit sphere.
    Parameters
    ----------
    ra_deg, dec_deg : float or arrays of floats, shape (N,)
         RA and Dec in degrees.
    Returns an array of floats with shape (N, 3).
    """
    ra  = np.asarray(ra_deg) * RAD_PER_DEG
    dec = np.asarray(dec_deg) * RAD_PER_DEG
    cosd = np.cos(dec)
    xyz = np.array([cosd * np.cos(ra),
                    cosd * np.sin(ra),
                    np.sin(dec)]).T

    return np.atleast_2d(xyz)
