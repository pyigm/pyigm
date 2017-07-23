#!/usr/bin/env python

"""
Script to generate a quick and dirty IGMSystem on disk (JSON file)
"""

import pdb
 
try:
    ustr = unicode
except NameError:
    ustr = str


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Show contents of a JSON file, assuming one of several formats. (v1.1)')
    parser.add_argument("itype", help="Type of IGMSystem: dla, lls, mgii")
    parser.add_argument("zabs", type=float, help="Absorption redshift")
    parser.add_argument("outfile", type=str, help="Name of JSON file to create")
    parser.add_argument("--NHI", type=float, help="log10 NHI value")
    parser.add_argument("--jcoord", type=str, help="Coordinates in JXXXXXXXX.X+XXXXXX.X format")
    parser.add_argument("--zem", type=float, help="Emission redshift")
    parser.add_argument("--sigNHI", type=float, help="Error in NHI")
    parser.add_argument("--vlim", type=str, help="Velocity limits in format ###,###")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from linetools import utils as ltu
    from pyigm.abssys.dla import DLASystem
    from pyigm.abssys.lls import LLSSystem
    from pyigm.abssys.igmsys import MgIISystem

    if args is None:
        pargs = parser()
    else:
        pargs = args

    # Coordinates
    if pargs.jcoord is not None:
        coord = ltu.radec_to_coord(pargs.jcoord)
    else:
        coord = SkyCoord(ra=0., dec=0., unit='deg')

    # vlim
    if pargs.vlim is not None:
        vlims = [float(vlim) for vlim in pargs.vlim.split(',')]*u.km/u.s
    else:
        vlims = None

    # go
    if pargs.itype.lower() == 'dla':
        isys = DLASystem(coord, pargs.zabs, vlims, pargs.NHI, zem=pargs.zem, sig_NHI=pargs.sigNHI)
    elif pargs.itype.lower() == 'lls':
        isys = LLSSystem(coord, pargs.zabs, vlims, NHI=pargs.NHI, zem=pargs.zem, sig_NHI=pargs.sigNHI)
    elif pargs.itype.lower() == 'mgii':
        isys = MgIISystem(coord, pargs.zabs, vlims, NHI=pargs.NHI, zem=pargs.zem, sig_NHI=pargs.sigNHI)
    else:
        raise IOError("Not prepared for this type of IGMSystem")

    # Write
    isys.write_json(pargs.outfile)

