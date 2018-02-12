""" Utilities for pyigm.igm
 Best to keep these separate from the Class modules
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str


def write_igmg_from_components(comp_list, specfile, fwhm, zem= 0., outfile='IGM_model.json'):
        """ Write to a JSON file of the IGMGuesses format from a list of AbsComponents.
        All coordinates must be equal!

        Parameters
        ----------
        comp_list : list
            list of AbsComponents
        specfile : str
            Name of spectrum file associated to these components
        fwhm : int
            FWHM of the spectrum
        outfile : str, optional
            Name of the output json file

        """
        import datetime
        # coordinates and meta
        coord_ref = comp_list[0].coord
        RA = coord_ref.ra.to('deg').value
        DEC = coord_ref.dec.to('deg').value
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        jname = ltu.name_from_coord(coord_ref, precision=(2, 1))
        # Create dict of the components
        out_dict = dict(cmps={},
                        spec_file=specfile,
                        fwhm=fwhm, bad_pixels=[],
                        meta={'RA': RA, 'DEC': DEC, 'ALTNAME':'unknown',
                              'zem': zem, 'Creator': 'pyigm.igm.utils.write_igmg_from_components()',
                              'Instrument': 'unknown', 'Date': date, 'JNAME': jname})

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
