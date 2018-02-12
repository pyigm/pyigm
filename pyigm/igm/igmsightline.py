""" Class for IGM Absorption sightline
 Uses AbsSightline from linetools
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import warnings

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssightline import AbsSightline
from linetools.isgm.abssystem import add_comps_from_dict
from linetools.isgm.utils import build_systems_from_components
from linetools import utils as ltu


class IGMSightline(AbsSightline):
    """Class for IGM Absorption Sightline
    """
    @classmethod
    def from_json(cls, jsonfile, **kwargs):
        """ Instantiate from a JSON file

        Parameters
        ----------
        jsonfile : str
          Filename
          See from_dict for required keys
        kwargs : passed to from_dict

        Returns
        -------

        """
        jdict = ltu.loadjson(jsonfile)
        slf = cls.from_dict(jdict, **kwargs)
        # Return
        return slf

    @classmethod
    def from_igmguesses_old(cls, radec, zem, igmgfile, name=None, **kwargs):
        """ Instantiate from a JSON file from IGMGuesses
        The input coordinates are used for all the components

        Parameters
        ----------
        radec : RA/DEC input
          See ltu.radec_to_coord for options
        zem : float
          Emission redshift of sightline
        igmgfile : str
          Filename

        Returns
        -------

        """
        # Read
        jdict = ltu.loadjson(igmgfile)   # cmps, specfile, meta
        # Add in additional keys
        coord = ltu.radec_to_coord(radec)
        jdict['RA'] = coord.fk5.ra.deg
        jdict['DEC'] = coord.fk5.dec.deg
        jdict['zem'] = zem
        # Name
        if name is None:
            name = 'J{:s}{:s}_z{:0.3f}'.format(
                coord.fk5.ra.to_string(unit=u.hour,sep='',pad=True)[0:4],
                coord.fk5.dec.to_string(sep='',pad=True,alwayssign=True)[0:5],
                zem)
        jdict['name'] = name
        jdict['components'] = jdict.pop('cmps')
        kwargs['use_coord'] = True
        slf = cls.from_dict(jdict, **kwargs)
        # Return
        return slf

    @classmethod
    def from_igmguesses(cls, igmgfile, name=None, **kwargs):
        """ Instantiate from a JSON file from IGMGuesses

        Parameters
        ----------
        igmgfile : str
          Filename

        Returns
        -------

        """
        # Read
        jdict = ltu.loadjson(igmgfile)  # cmps, specfile, meta
        # Add in additional keys

        jdict['RA'] = jdict['meta']['RA'] * u.deg
        jdict['DEC'] = jdict['meta']['DEC'] * u.deg
        coord = SkyCoord(jdict['RA'], jdict['DEC'], unit='deg')
        jdict['zem'] = jdict['meta']['zem']
        if jdict['zem'] == 0.:  # conforming IGMGuesses current rule zem = 0. means zem not set
            jdict['zem'] = None
        # Name
        if name is None:
            if jdict['zem'] is None:
                zem = 'z-unknown'
            else:
                zem = 'z{:0.3f}'.format(jdict['zem'])
            name = 'J{:s}{:s}_{:s}'.format(
                coord.fk5.ra.to_string(unit=u.hour, sep='', pad=True)[0:4],
                coord.fk5.dec.to_string(sep='', pad=True, alwayssign=True)[0:5], zem)
        jdict['name'] = name
        jdict['components'] = jdict.pop('cmps')
        kwargs['use_coord'] = True
        slf = cls.from_dict(jdict, **kwargs)
        # Return
        return slf

    @classmethod
    def from_dict(cls, idict, **kwargs):
        """ Instantiate from a dict

        Parameters
        ----------
        idict : dict
          Required keys are:
           'RA' -- float (deg)
           'DEC' -- float(deg)
           'zem' -- float
           'name' -- str
           'components' -- list
         Other keys are added as attributes to the IgmSightline object

        Returns
        -------

        """
        slf = cls(SkyCoord(ra=idict['RA'], dec=idict['DEC'], unit='deg'),
                  zem=idict['zem'], name=idict['name'], **kwargs)
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

    def __init__(self, radec, zem=None, **kwargs):
        AbsSightline.__init__(self, radec, sl_type='IGM', **kwargs)
        if zem is None:
            self.zem = None
            warnings.warn("You really should set zem to create IGMSightline")
        else:
            self.zem = zem

    def make_igmsystems(self, igmsystem=None, **kwargs):
        """ Use the component list to generate a list of IGMSystems

        Returns
        -------
        igm_systems : list
          list of IGMSystem objects

        """
        if igmsystem is None:
            from pyigm.abssys.igmsys import IGMSystem
            igmsystem = IGMSystem
        # Main call
        igm_sys = build_systems_from_components(self._components, systype=igmsystem, **kwargs)
        # Return
        return igm_sys

    def write_to_igmguesses(self, outfile, fwhm=3, specfilename=None, creator='unknown', instrument='unknown',
                            altname='unknown'):
        import json
        """
        Writes an IGMGuesses formatted JSON file

        Parameters
        ----------
        outfile : str
            Name of the IGMGuesses JSON file to write to.
        fwhm : int
            FWHM for IGMguesses
        specfilename : str
            Name of the spectrum file these guesses are associated to.
        altname : str
            Alternative name for the sightline. e.g. 3C273
        creator : str
            Name of the person who is creating the file. Important for tracking.
        instrument : str
            String indicating the instrument and its configuration associated to the specfilename.
            e.g. HST/COS/G130M+G160M/LP2


        Returns
        -------

        """
        import datetime
        # components
        comp_list = self._components

        # coordinates and meta
        coord_ref = comp_list[0].coord
        RA = coord_ref.ra.to('deg').value
        DEC = coord_ref.dec.to('deg').value
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        jname = ltu.name_from_coord(coord_ref, precision=(2, 1))
        if self.zem is None:
            zem = 0.  # IGMGuesses rules
        else:
            zem = self.zem
        # Create dict of the components
        out_dict = dict(cmps={},
                        spec_file=specfilename,
                        fwhm=fwhm, bad_pixels=[],
                        meta={'RA': RA, 'DEC': DEC, 'ALTNAME': altname,
                              'zem': zem, 'Creator': creator,
                              'Instrument': instrument, 'Date': date, 'JNAME': jname})

        for comp in comp_list:
            key = comp.name
            out_dict['cmps'][key] = comp.to_dict()
            # import pdb; pdb.set_trace()
            # check coordinate
            if comp.coord != coord_ref:
                raise ValueError("All AbsComponent objects must have the same coordinates!")
            out_dict['cmps'][key]['zcomp'] = comp.zcomp
            out_dict['cmps'][key]['zfit'] = comp.zcomp
            out_dict['cmps'][key]['Nfit'] = comp.logN
            out_dict['cmps'][key]['bfit'] = comp.attrib['b']
            out_dict['cmps'][key]['wrest'] = comp._abslines[0].wrest.value
            out_dict['cmps'][key]['vlim'] = list(comp.vlim.value)
            out_dict['cmps'][key]['reliability'] = str(comp.reliability)
            out_dict['cmps'][key]['comment'] = str(comp.comment)
            # out_dict['cmps'][key]['mask_abslines'] = comp.mask_abslines

        # JSONify
        gd_dict = ltu.jsonify(out_dict)

        # Write file
        # with io.open(outfile, 'w', encoding='utf-8') as f:
        f = open(outfile, 'w')
        f.write(unicode(json.dumps(gd_dict, sort_keys=True, indent=4, separators=(',', ': '))))
        print('Wrote: {:s}'.format(outfile))


    def __repr__(self):
        if self.zem is None:
            zem = 'zem = unknown'
        else:
            zem = 'zem={:f}'.format(self.zem)
        txt = '<{:s}: {:s} {:s}, {:s}'.format(
                self.__class__.__name__, self.coord.fk5.ra.to_string(unit=u.hour,sep=':', pad=True),
                self.coord.fk5.dec.to_string(sep=':',pad=True,alwayssign=True), zem)

        # Type?
        if self.em_type is not None:
            txt = txt + ', em_type={:s}'.format(self.em_type)

        # Finish
        txt += '>'
        return (txt)
