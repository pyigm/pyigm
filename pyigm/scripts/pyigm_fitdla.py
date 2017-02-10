#!/usr/bin/env python

"""
Script to run igmguesses GUI
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from PyQt4 import QtGui
from pyigm.guis.igmguesses import IGMGuessesGui
from astropy import units as u
import pdb

try:
    ustr = unicode
except NameError:
    ustr = str


def parser(options=None):
    """ Parser for igmguesses
    Parameters
    ----------
    options

    Returns
    -------

    """

    import argparse
    parser = argparse.ArgumentParser(description='Parser for FitDLAGUI')
    parser.add_argument("in_file", type=str, help="Spectral file")
    parser.add_argument("zqso", type=float, help="Use QSO template with zqso")
    parser.add_argument("-out_file", type=str, help="Output LLS Fit file")
    parser.add_argument("-smooth", type=float, help="Smoothing (pixels)")
    parser.add_argument("-dla_fit_file", type=str, help="Input LLS Fit file")
    parser.add_argument("-conti_file", type=str, help="Input continuum spectrum")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):

    import sys
    from pyigm.guis.fitdla import XFitDLAGUI
    if args is None:
        pargs = parser(options=args)
    else: # better know what you are doing!
        pdb.set_trace()  # Not ready for this yet
        if isinstance(args[0],(Spectrum1D,tuple)):
            app = QtGui.QApplication(sys.argv)
            gui = XFitDLAGUI(args[0], **kwargs)
            gui.show()
            app.exec_()
            return
        else: # String parsing
            largs = ['1'] + [iargs for iargs in args]
            pargs = parser.parse_args(largs)


    # Run
    app = QtGui.QApplication(sys.argv)
    gui = XFitDLAGUI(pargs.in_file,outfil=pargs.out_file,smooth=pargs.smooth,
        dla_fit_file=pargs.dla_fit_file, zqso=pargs.zqso, conti_file=pargs.conti_file)
    gui.show()
    app.exec_()

if __name__ == '__main__':
    main()
