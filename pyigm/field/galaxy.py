""" Simple Class for a Galaxy
  Likely to be replaced with an external Class
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import pdb

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.utils import radec_to_coord
from linetools import utils as ltu

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
    @classmethod
    def from_dict(cls, idict):
        """ Generate a Galaxy object from a dict

        Parameters
        ----------
        idict : dict

        """
        slf = cls(SkyCoord(ra=idict['RA'], dec=idict['DEC'], unit=u.deg),
                  name=idict['Name'])
        # Attributes
        for key in idict.keys():
            if key in ['RA', 'DEC', 'CreationDate', 'user']:
                continue
            try:
                setattr(slf,key,idict[key])
            except KeyError:
                pass
        # Return
        return slf

    def __init__(self, radec, z=None, name=None):
        self.coord = radec_to_coord(radec)
        # Redshift
        self.z = z
        
        # Name
        if name is None:
            self.name = ('J'+self.coord.ra.to_string(unit=u.hour, sep='', pad=True)+
                        self.coord.dec.to_string(sep='', pad=True, alwayssign=True))
        else:
            self.name = name

    def to_dict(self):
        """ Convert the galaxy to a JSON-ready dict for output

        Returns
        -------
        gdict : dict

        """
        import datetime
        import getpass
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        user = getpass.getuser()
        # Generate the dict
        gdict = dict(Name=self.name,
                       RA=self.coord.ra.value,
                       DEC=self.coord.dec.value,
                       CreationDate=date,
                       user=user
                       )
        # Attributes (e.g. SFR)
        for key in self.__dict__.keys():
            if key in ['coord', 'name']:
                continue
            gdict[key] = getattr(self, key)
        # Polish
        gdict = ltu.jsonify(gdict)
        # Return
        return gdict

    # #############
    def __repr__(self):
        return ('<Galaxy: {:s} {:s} {:s}, z={:g}>'.format(
                 self.name, 
                 self.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                 self.coord.dec.to_string(sep=':', pad=True, alwayssign=True),
                 self.z))



