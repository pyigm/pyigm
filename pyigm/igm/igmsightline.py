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
    def from_igmguesses(cls, igmgfile, name=None, radec=None, zem=None, **kwargs):
        """ Instantiate from a JSON file from IGMGuesses

        Parameters
        ----------
        igmgfile : str
          Filename
        name : str, optional
          Name of the IGMSightline
        radec : RA/DEC input, optional
          See ltu.radec_to_coord for options on format
          If given, it will overwrite the RA/DEC in IGMGuesses file (if any)
        zem : float, optional
          Emission redshift of sightline
          If given, it will overwrite the zem in IGMGuesses file (if any)

        Returns
        -------

        """
        from linetools.isgm import abscomponent
        # Read
        jdict = ltu.loadjson(igmgfile)  # cmps, specfile, meta
        # Add in additional keys
        # Coords
        if radec is not None:
            coord = ltu.radec_to_coord(radec)
        else:
            coord = SkyCoord(jdict['meta']['RA'], jdict['meta']['DEC'], unit='deg')
        jdict['RA'] = coord.icrs.ra.deg
        jdict['DEC'] = coord.icrs.dec.deg
        # zem
        if zem is not None:
            jdict['zem'] = zem
        else:
            jdict['zem'] = jdict['meta']['zem']
            if jdict['zem'] == 0.:  # conforming IGMGuesses current rule zem = 0. means zem not set
                jdict['zem'] = None
        # Name
        if name is None:
            if jdict['zem'] is None:
                zem_name = 'z-unknown'
            else:
                zem_name = 'z{:0.3f}'.format(jdict['zem'])
            name = 'J{:s}{:s}_{:s}'.format(
                coord.icrs.ra.to_string(unit=u.hour, sep='', pad=True)[0:4],
                coord.icrs.dec.to_string(sep='', pad=True, alwayssign=True)[0:5], zem_name)
        jdict['name'] = name
        # Components
        jdict['components'] = jdict.pop('cmps')

        kwargs['use_coord'] = True
        slf = cls.from_dict(jdict, **kwargs)

        # Separate IGMGuesses attributes from AbsComponent
        acomp_keys = list(abscomponent.init_attrib.keys())
        for comp in slf._components:
            comp.igmg_attrib = {}
            akeys = list(comp.attrib.keys())
            for key in akeys:
                if key in acomp_keys:
                    pass
                else:
                    comp.igmg_attrib[key] = comp.attrib.pop(key)
            # Add other IGMGuesses specific attributes
            comp.igmg_attrib['top_level'] = {}
            for key in ['Nfit', 'bfit', 'zfit', 'mask_abslines', 'wrest', 'vlim']:  # I added vlim because these aren't always consistent. Must be something bad in igmguesses
                if key in jdict['components'][comp.name].keys():
                    comp.igmg_attrib['top_level'][key] = jdict['components'][comp.name][key]
        # Slurp a few other IGMGuesses things
        slf.igmg_dict = {}
        for key in ['spec_file', 'meta', 'fwhm']:
            slf.igmg_dict[key] = jdict[key]
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
        if 'linelist' not in kwargs.keys():
            from linetools.lists.linelist import LineList
            ism = LineList('ISM')
            kwargs['linelist'] = ism
        from pyigm.abssys.utils import class_by_type
        ism = LineList('ISM')
        kwargs['linelist'] = ism
        # Load ISM to speed things up
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

        # Systems
        if 'systems' in idict.keys():
            for key in idict['systems'].keys():
                asys = class_by_type(idict['systems'][key]['abs_type']).from_dict(idict['systems'][key], **kwargs)
                slf._abssystems.append(asys)
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

    def write_to_json(self, outfile):
        # Generate the dict
        igms_dict = self.to_dict()
        # Jsonify
        clean_dict = ltu.jsonify(igms_dict)
        # Write
        ltu.savejson(outfile, clean_dict, overwrite=True)

    def write_to_igmguesses(self, outfile, fwhm=3., specfilename=None, creator=None,
                            instrument='unknown',
                            altname='unknown', date=None, overwrite=False):
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
        overwrite : bool
            Wthether to overwrite output


        Returns
        -------

        """
        import datetime
        # Slurp IGMGuesses component attributes
        from pyigm.guis.igmguesses import comp_init_attrib
        from linetools.isgm.abscomponent import AbsComponent
        tmp = AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s)
        comp_init_attrib(tmp)
        igmg_akeys = list(tmp.attrib.keys())
        # components
        comp_list = self._components
        coord_ref = comp_list[0].coord
        if date is None:
            date = str(datetime.date.today().strftime('%Y-%b-%d'))

        # spec_file, meta
        if hasattr(self, 'igmg_dict'):
            spec_file = self.igmg_dict['spec_file']
            meta = self.igmg_dict['meta']
            # Updates
            meta['Date'] = date
            if creator is not None:
                meta['Creator'] = creator
            #
            fwhm = self.igmg_dict['fwhm']
        else:
            spec_file = specfilename
            # coordinates and meta
            RA = coord_ref.ra.to('deg').value
            DEC = coord_ref.dec.to('deg').value
            jname = ltu.name_from_coord(coord_ref, precision=(2, 1))
            if self.zem is None:
                zem = 0.  # IGMGuesses rules
            else:
                zem = self.zem
            if creator is None:
                creator = 'unknown'
            meta = {'RA': RA, 'DEC': DEC, 'ALTNAME': altname,
                    'zem': zem, 'Creator': creator,
                    'Instrument': instrument, 'Date': date, 'JNAME': jname}

        # Create dict of the components
        out_dict = dict(cmps={},
                        spec_file=spec_file,
                        fwhm=fwhm, bad_pixels=[],
                        meta=meta)

        for comp in comp_list:
            key = comp.name
            out_dict['cmps'][key] = comp.to_dict()
            # import pdb; pdb.set_trace()
            # check coordinate
            if comp.coord != coord_ref:
                raise ValueError("All AbsComponent objects must have the same coordinates!")
            out_dict['cmps'][key]['zcomp'] = comp.zcomp
            # IGMGuesses specific component attr
            for igm_key in ['zfit','Nfit','bfit', 'wrest', 'mask_abslines', 'vlim']:
                out_dict['cmps'][key][igm_key] = comp.igmg_attrib['top_level'][igm_key]
            # IGMGuesses attribute dict
            out_dict['cmps'][key]['attrib'] = {}
            for igm_key in igmg_akeys:
                try:
                    out_dict['cmps'][key]['attrib'][igm_key] = comp.igmg_attrib[igm_key]
                except KeyError:
                    out_dict['cmps'][key]['attrib'][igm_key] = comp.attrib[igm_key]
            #out_dict['cmps'][key]['vlim'] = list(comp.vlim.value)
            out_dict['cmps'][key]['reliability'] = str(comp.reliability)
            out_dict['cmps'][key]['comment'] = str(comp.comment)
            # Compatability on sig_logN
            out_dict['cmps'][key]['attrib']['sig_logN'] = comp.attrib['sig_logN'][0]

        # JSONify
        gd_dict = ltu.jsonify(out_dict)

        # Write file
        ltu.savejson(outfile, gd_dict, overwrite=overwrite, sort_keys=True, indent=4, separators=(',', ': '))
        print('Wrote: {:s}'.format(outfile))


    def __repr__(self):
        if self.zem is None:
            zem = 'zem = unknown'
        else:
            zem = 'zem={:f}'.format(self.zem)
        txt = '<{:s}: {:s} {:s}, {:s}'.format(
                self.__class__.__name__, self.coord.icrs.ra.to_string(unit=u.hour,sep=':', pad=True),
                self.coord.icrs.dec.to_string(sep=':',pad=True,alwayssign=True), zem)

        # Type?
        if self.em_type is not None:
            txt = txt + ', em_type={:s}'.format(self.em_type)

        # Finish
        txt += '>'
        return (txt)
