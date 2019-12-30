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

    parser = argparse.ArgumentParser(description='Parser for IGMGuesses')
    parser.add_argument("in_file", type=str, help="Spectral file")
    parser.add_argument("-o", "--out_file", type=str, help="Output Guesses .json file [Default is IGM_model.json]")
    parser.add_argument("--fwhm", type=float, help="FWHM Gaussian smoothing for fitting (pixels)")
    parser.add_argument("-p", "--previous_file", type=str, help="Input Guesses .json file [previously generated]")
    parser.add_argument("--n_max_tuple", type=int, help="Maximum number of transitions per ion species to display")
    parser.add_argument("--min_strength", type=float, help="Minimum strength for transitions to be displayed; choose values (0,14.7)")
    parser.add_argument("--min_ew", type=float, help="Minimum EW (in AA) for transitions to be stored within a component. This\
                                                    is useful to get rid of extremely weak transitions from the model")
    parser.add_argument("--vlim", type=float, help="Velocity limit (in km/s) for the display. This limit will apply to both sides")
    parser.add_argument("--external_model", type=str, help="An external model spectrum (.fits)")
    parser.add_argument("--scale", type=float, help="Scaling of screen [default=1.]")
    parser.add_argument("-redsh","--redsh",type=float,help="Redshift at which transitions are initially plotted")
    parser.add_argument("-wmin", "--wmin", type=float, help="Minimum wavelenght for which line transitions are calculated")
    parser.add_argument("-wmax", "--wmax", type=float, help="Maximum wavelenght for which line transitions are calculated")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(namespace=options)
    return args


def main(args=None):
    pargs = parser(options=args)
    import sys
    from pyigm.guis.igmguesses import IGMGuessesGui
    from PyQt5.QtWidgets import QApplication
    from astropy import units as u

    if pargs.vlim is not None:
        vlim_disp = [-1*pargs.vlim, 1.*pargs.vlim]*u.km/u.s
    else:
        vlim_disp = pargs.vlim
    if pargs.scale is not None:
        scale = pargs.scale
    else:
        scale = 1.


    app = QApplication(sys.argv)
    gui = IGMGuessesGui(pargs.in_file,
                        outfil=pargs.out_file,
                        fwhm=pargs.fwhm,
                        previous_file=pargs.previous_file,
                        n_max_tuple=pargs.n_max_tuple,
                        min_strength=pargs.min_strength,
                        min_ew=pargs.min_ew,
                        screen_scale=scale,
                        vlim_disp=vlim_disp,
                        external_model=pargs.external_model,
                        redsh=pargs.redsh,
                        spwvmin=pargs.wmin, spwvmax=pargs.wmax)
    gui.show()
    app.exec_()


if __name__ == '__main__':
    main()
