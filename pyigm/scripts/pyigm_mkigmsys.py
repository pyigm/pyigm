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
        description='Show contents of a JSON file, assuming one of several formats.')
    parser.add_argument("itype", help="Type of IGMSystem: dla, lls")
    parser.add_argument("zabs", type=float, help="Absorption redshift")
    parser.add_argument("NHI", type=float, help="log10 NHI value")
    parser.add_argument("outfile", type=str, help="Name of JSON file to create")
    parser.add_argument("--jcoord", type=str, help="Coordiantes in JXXXXXXXX.X+XXXXXX.X format")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):
    from astropy.coordinates import SkyCoord
    from linetools import utils as ltu
    from pyigm.abssys.dla import DLASystem
    from pyigm.abssys.lls import LLSSystem

    pargs = parser()

    # Coordinates
    if pargs.jcoord is not None:
        coord = ltu.radec_to_coord(pargs.jcoord)
    else:
        coord = SkyCoord(ra=0., dec=0., unit='deg')

    # go
    if pargs.itype == 'dla':
        isys = DLASystem(coord, pargs.zabs, None, pargs.NHI)
    elif pargs.itype == 'lls':
        isys = LLSSystem(coord, pargs.zabs, None, NHI=pargs.NHI)
    else:
        raise IOError("Not prepared for this type of IGMSystem")

    # Write
    isys.write_json(pargs.outfile)

