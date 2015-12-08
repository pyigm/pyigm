""" Utilities for IGM calculations
"""

def cosm_xz(z, cosmo=None, flg_return=0):
    """ Calculates X(z) -- absorption path length or dXdz

    Parameters
    ----------
    z : float or ndarray
       Redshift to evaluate at.  May be an array
    cosmo : astropy.cosmology
       Cosmological model to adopt
    flg_return : int, optional (0)
       Flag controlling the output
         0 = X(z)
         1 = dX/dz at z

    Returns
    -------
    Xz or dXdz :
      Absorption path or its derivative
      Set by flag
    """

    # Cosmology
    if cosmo is None:
        from astropy.cosmology import core as acc
        cosmo = acc.FlatLambdaCDM(70., 0.3)

    # Flat?
    if cosmo.Ok(0.) == 0:
        if flg_return == 0:  # X(z)
            rslt = cosmo.absorption_distance(z)
        elif flg_return == 1:  # dX/dz
            rslt = cosmo.abs_distance_integrand(z)
        else:
            raise ValueError('pyigm.utils.cosm_xz: Bad flg {:d}'.format(flg_return))
    else:
        raise ValueError('pyigm.utils.cosm_xz: Not prepared for non-flat cosmology')

    #
    return rslt