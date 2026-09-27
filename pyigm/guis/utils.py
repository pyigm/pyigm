""" Helpers shared by the pyigm GUIs that do not themselves need Qt

Keeping these here means headless code paths (e.g.
:meth:`pyigm.igm.igmsightline.IGMSightline.write_to_igmguesses`) can use them
without importing a GUI toolkit.
"""

from astropy import units as u


def comp_init_attrib(comp):
    """ Initialize the attrib dict of an AbsComponent to IGMGuesses defaults

    Parameters
    ----------
    comp : AbsComponent
      Modified in place
    """
    # Attributes
    comp.attrib = {'N': 0./u.cm**2, 'Nsig': 0./u.cm**2, 'flagN': 0,  # Column
               'logN': 0., 'sig_logN': 0.,
               'b': 0.*u.km/u.s, 'bsig': 0.*u.km/u.s,  # Doppler
               'z': comp.zcomp, 'zsig': 0.,
               'reliability': 'none'}
