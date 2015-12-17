""" Simple Class for a Galaxy
  Likely to be replaced with an external Class
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from astropy.coordinates import SkyCoord
from astropy import units as u

class Galaxy(object):
    """A Galaxy Class

    Parameters
    ----------
    radec : tuple or SkyCoord
      (RA,DEC) in deg or astropy.coordinate
    z : float, optional
      Redshift

    Attributes
    ----------
    name : str
        Name(s)
    z : float, optional
       Adopted redshift
    coord : SkyCoord
    """
    # Initialize with a .dat file
    def __init__(self, radec, z=None):
        if isinstance(radec, (tuple)):
            self.coord = SkyCoord(ra=radec[0], dec=radec[1])
        elif isinstance(radec, SkyCoord):
            self.coord = radec
        # Redshift
        self.z = z
        
        # Name
        self.name = ('J'+self.coord.ra.to_string(unit=u.hour, sep='', pad=True)+
                    self.coord.dec.to_string(sep='', pad=True, alwayssign=True))

    # #############
    def __repr__(self):
        return ('<Galaxy: {:s} {:s} {:s}, z={:g}>'.format(
                 self.name, 
                 self.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                 self.coord.dec.to_string(sep=':', pad=True, alwayssign=True),
                 self.z))



