#!/usr/bin/env python

"""
Script to run igmguesses GUI
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from PyQt4 import QtGui
from pyigm.guis.igmguesses import IGMGuessesGui
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

    parser = argparse.ArgumentParser(description='Parser for IGMGuesses')
    parser.add_argument("in_file", type=str, help="Spectral file")
    parser.add_argument("-o", "--out_file", type=str, help="Output Guesses file [Default is IGM_model.json]")
    parser.add_argument("--fwhm", type=float, help="FWHM smoothing for fitting (pixels)")
    parser.add_argument("-p", "--previous_file", type=str, help="Input Guesses file [previously generated]")
    parser.add_argument("--n_max_tuple", type=int, help="Maximum number of transitions per ion species to display")
    parser.add_argument("--min_strength", type=float, help="Minimum strength for transitions to be displayed; choose values (0,14.7)")
    parser.add_argument("--min_ew", type=float, help="Minimum EW (in AA) for transitions to be stored within a component. This\
                                                    is useful to get rid of extremely weak transitions from the model")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    pargs = parser(options=args)
    import sys

    app = QtGui.QApplication(sys.argv)
    gui = IGMGuessesGui(pargs.in_file,
                        outfil=pargs.out_file,
                        fwhm=pargs.fwhm,
                        previous_file=pargs.previous_file,
                        n_max_tuple=pargs.n_max_tuple,
                        min_strength=pargs.min_strength,
                        min_ew=pargs.min_ew)
    gui.show()
    app.exec_()


if __name__ == '__main__':
    main()
