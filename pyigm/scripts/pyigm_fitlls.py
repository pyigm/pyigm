#!/usr/bin/env python

"""
Script to run fitlls GUI
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:
    ustr = unicode
except NameError:
    ustr = str


# Script to run XSpec from the command line or ipython
def parser(options=None):
    '''
    Parser for XFitLLSGUI

    Command line or from Python
    Examples:
      1.  python ~/xastropy/xastropy/xguis/spec_guis.py 1
      2.  spec_guis.run_fitlls(filename)
      3.  spec_guis.run_fitlls(spec1d)
    '''

    import argparse

    parser = argparse.ArgumentParser(description='Parser for FitLLSGUI (v1.0)')
    parser.add_argument("in_file", type=str, help="Spectral file")
    parser.add_argument("zqso", type=float, help="Use Telfer template with zqso")
    parser.add_argument("--out_file", type=str, help="Output LLS Fit file")
    parser.add_argument("--smooth", type=float, default=3., help="Smoothing (pixels)")
    parser.add_argument("--lls_fit_file", type=str, help="Input LLS Fit file")
    parser.add_argument("--norm", help="Use normalized spectrum", action="store_true")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(namespace=options)
    return args

def main(args=None):

    import sys
    from PyQt5.QtWidgets import QApplication
    from pyigm.guis.fitlls import XFitLLSGUI
    pargs = parser(options=args)

    # Run
    app = QApplication(sys.argv)
    gui = XFitLLSGUI(pargs.in_file, pargs.zqso, outfil=pargs.out_file,smooth=pargs.smooth,
                     lls_fit_file=pargs.lls_fit_file, norm=args.norm)
    gui.show()
    app.exec_()


# ################
if __name__ == "__main__":
    main()

