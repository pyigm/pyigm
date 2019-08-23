#!/usr/bin/env python

"""
Script to run igmguesses GUI
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

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
    parser = argparse.ArgumentParser(description='Parser for FitDLAGUI (v1.0)')
    parser.add_argument("in_file", type=str, help="Spectral file")
    parser.add_argument("zqso", type=float, help="Use QSO template with zqso")
    parser.add_argument("--out_file", type=str, help="Output DLA Fit file")
    parser.add_argument("--smooth", type=float, help="Smoothing (pixels)")
    parser.add_argument("--dla_fit_file", type=str, help="Input DLA Fit file")
    parser.add_argument("--conti_file", type=str, help="Input continuum spectrum")
    parser.add_argument("--zdla", type=float, help="Input DLA redshift")
    parser.add_argument("--NHI", type=float, help="Input DLA NHI")
    parser.add_argument("--norm", help="Use normalized spectrum", action="store_true")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(namespace=options)
    return args


def main(args=None):

    import sys
    from PyQt5.QtWidgets import QApplication
    from pyigm.guis.fitdla import XFitDLAGUI
    pargs = parser(options=args)

    if hasattr(args, 'norm'):
        inorm = args.norm
    else:
        inorm = False

    # Run
    app = QApplication(sys.argv)
    gui = XFitDLAGUI(pargs.in_file,outfil=pargs.out_file,smooth=pargs.smooth,
        dla_fit_file=pargs.dla_fit_file, zqso=pargs.zqso, conti_file=pargs.conti_file,
                     zdla=pargs.zdla, NHI=pargs.NHI, norm=inorm)
    gui.show()
    app.exec_()

if __name__ == '__main__':
    main()
