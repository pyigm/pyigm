""" Module for core continuum codes
"""
from __future__ import print_function, absolute_import, division, unicode_literals


def init_conti_dict(Norm=0., tilt=0., tilt2=0.,
                    piv_wv=0., piv_wv2=None, igm='None',
                    fN_gamma=-1., LL_flatten='True'):
    """Initialize a continuum conti_dict

    Parameters
    ----------
    Norm : float, optional
      Normaliztion
    tilt : float, optional
      Power-law tilt to continuum
    piv_wv : float, optional
      Pivot wave for tilt.  Best kept *without* units
    piv_wv2 : float, optional
      Pivot wave for a second tilt. Better be at wavelength < piv_wv
    igm : str, optional
      Adopt average IGM model? ['None']
    LL_flatten : bool, optional
      Set Telfer to a constant below the LL?

    Returns
    -------
    conti_dict : dict
      Useful for simple modeling.  Keep as a dict for JSON writing
    """
    conti_dict = dict(Norm=Norm, tilt=tilt, piv_wv=piv_wv, piv_wv2=piv_wv2,
                      tilt2=tilt2, igm=igm, fN_gamma=fN_gamma,
                      LL_flatten=LL_flatten)
    # Checks
    if piv_wv2 is None:
        conti_dict.pop('piv_wv2')
    else:
        if piv_wv is None:
            raise IOError("piv_wv required if piv_wv2 set")
        else:
            if piv_wv2 > piv_wv:
                raise ValueError("piv_wv2 < piv_wv required!")
    #
    return conti_dict



