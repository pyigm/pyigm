""" Class for IGM Absorption sightline
Uses AbsSightline from linetools
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from astropy import units as u
from linetools.isgm.abssightline import AbsSightline


class IGMSightline(AbsSightline):
    """Class for IGM Absorption Sightline
    """
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
