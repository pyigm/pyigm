""" Utilities for IGM calculations
"""
import numpy as np
import pdb

from astropy import constants as const
from astropy import units as u
from astropy.units.quantity import Quantity

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
            try:
                rslt = cosmo.abs_distance_integrand(z)
            except AttributeError:
                rslt = cosmo._xfunc(z)
        else:
            raise ValueError('pyigm.utils.cosm_xz: Bad flg {:d}'.format(flg_return))
    else:
        raise ValueError('pyigm.utils.cosm_xz: Not prepared for non-flat cosmology')

    #
    return rslt

def mk_ew_lyman_spline(bval, ew_fil=None, chk=False):
    """ Generate a pickle file of a Spline of EW vs NHI for the Lyman series

    Parameters
    ----------
    bval : float or Quantity
      Doppler parameter (km/s)
    ew_fil : string ('EW_SPLINE_b##.p')
      Name of output pickle file
    chk : bool, optional
        Check the calculation

    """
    import pickle
    import yaml
    from linetools.lists.linelist import LineList
    from linetools.analysis import voigt as lav
    from scipy import interpolate

    # Outfil
    if ew_fil == None:
        ew_fil = 'EW_SPLINE_b'+str(int(bval))+'.p'

    # Units
    if not isinstance(bval,u.quantity.Quantity):
        bval = bval * u.km/u.s  # km/s

    # NHI
    nspl = 100
    log_NHI = 11.0 + 11*np.arange(nspl)/(nspl-1.)

    # Lines
    HI = LineList('HI')
    wrest = HI._data['wrest']

    # Output
    owrest = [row['wrest'].value for row in HI._data]
    outp = {'wrest': owrest, 'tck': []}

    # Setup
    nvel = 60001
    velo = (-30000. + np.arange(nvel,dtype='float64'))*u.km/u.s # km/s
    dvel = 1. * u.km/u.s # km/s
    uval = velo / bval

    # Loop
    for cnt, line in enumerate(HI._data):

        # Wave array
        dwv = dvel.to(u.cm/u.s) * line['wrest'] / const.c.cgs  # Ang

        # Voigt
        vd = (bval/line['wrest']).to(u.Hz)  # Frequency
        a = line['gamma'].value / (12.56637 * vd.value)
        vgt = lav.voigt_wofz(uval.value, a)

        # tau
        tau = 0.014971475*line['f']*vgt/vd  # Normalized to N_HI = 1 cm^-2

        # Flux
        tau_array = np.outer(tau, 10.**log_NHI)
        fx = np.exp(-1.*tau_array)

        # EW
        EW = np.sum(1.-fx, 0) * dwv

        # Spline
        tck = interpolate.splrep(log_NHI, EW)

        # Check?
        if chk:
            from matplotlib import pyplot as plt
            plt.clf()
            plt.plot(log_NHI, EW, 'o')
            # Spline
            xnew = np.linspace(np.amin(log_NHI),np.amax(log_NHI), nspl*10)
            ynew = interpolate.splev(xnew, tck, der=0)
            plt.plot(xnew, ynew, '-')
            plt.show()

        # Output
        print('line = {:g}'.format(line['wrest']))
        outp['tck'].append(tck)

    # Write
    print('Writing {:s}'.format(ew_fil))
    with open(ew_fil, 'w') as yamlf:
        yamlf.write(yaml.dump(outp))

def lst_to_array(lst, mask=None):
    """ Simple method to convert a list to an array

    Allows for a list of Quantity objects

    Parameters
    ----------
    lst : list
      Should be number or Quantities
    mask : boolean array, optional

    Returns
    -------
    array or Quantity array

    """
    if mask is None:
        mask = np.array([True]*len(lst))
    if isinstance(lst[0], Quantity):
        return Quantity(lst)[mask]
    else:
        return np.array(lst)[mask]
