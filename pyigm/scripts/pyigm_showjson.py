#!/usr/bin/env python

"""
Script to show the contents of a JSON file
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
    parser.add_argument("jsonfile", help="Name of JSON file")
    #parser.add_argument("--epoch", default=2000., type=float, help="Epoch [Not functional]")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):
    from linetools import utils as ltu

    pargs = parser()
    # Read
    jdict = ltu.loadjson(pargs.jsonfile)

    if 'class' not in jdict.keys():
        raise KeyError("This script only works with JSON files with named classes")
    if jdict['class'] == 'IGMSightline':
        from pyigm.igm.igmsightline import IGMSightline
        obj = IGMSightline.from_dict(jdict)
    else:
        raise IOError("Not prepared for this class: {:s}".format(jdict['class']))

    # name
    try:
        name = jdict['name']
    except KeyError:
        try:
            name = jdict['Name']
        except:
            name = 'None'
    print("Name of object: {:s}".format(name))

    # Generate table
    tbl = obj.build_table()
    if len(tbl) > 0:
        tbl.pprint(max_width=120)
    else:
        print("Table was empty..")

