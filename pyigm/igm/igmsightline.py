""" Class for IGM Absorption sightline
 Uses AbsSightline from linetools
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssightline import AbsSightline
from linetools.isgm.abssystem import add_comps_from_dict


class IGMSightline(AbsSightline):
    """Class for IGM Absorption Sightline
    """
    @classmethod
    def from_dict(cls, idict, coord=None, **kwargs):
        """ Instantiate from a dict

        Parameters
        ----------
        idict : dict

        Returns
        -------

        """
        slf = cls(SkyCoord(ra=idict['RA'], dec=idict['DEC'], unit='deg'),
                  idict['zem'], name=idict['name'], **kwargs)
        # Other
        for key in idict.keys():
            if key in ['RA', 'DEC', 'zem', 'name', 'components']:
                continue
            else:
                setattr(slf, key, idict[key])
        # Components
        add_comps_from_dict(slf, idict, **kwargs)

        # Return
        return slf

    def __init__(self, radec, zem, **kwargs):
        AbsSightline.__init__(self, radec, sl_type='IGM', **kwargs)
        self.zem = zem

    def __repr__(self):
        txt = '<{:s}: {:s} {:s}, zem={:f}'.format(
                self.__class__.__name__, self.coord.ra.to_string(unit=u.hour,sep=':', pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True), self.zem)

        # Type?
        if self.em_type is not None:
            txt = txt + ', em_type={:s}'.format(self.em_type)

        # Finish
        txt += '>'
        return (txt)
