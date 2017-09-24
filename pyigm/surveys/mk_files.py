""" Script to generate output files from IGM Surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Execute with:  python mk_files -h


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Create one or more files for IGM Surveys.')
    parser.add_argument("--dla_lz", help="Create dla_lz_boot file", action="store_true")
    parser.add_argument("--all", help="Create dla_lz_boot file", action="store_true")
    parser.add_argument("--nproc", type=int, default=4, help="Number of processors")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):

    from pkg_resources import resource_filename
    from pyigm.surveys.dlasurvey import fit_atan_dla_lz

    pargs = parser()
    if pargs.dla_lz or pargs.all:
        """ Creates dla_lz_boot.fits
        """
        outfile = resource_filename('pyigm', 'data/DLA/dla_lz_boot.fits.gz')
        fit_atan_dla_lz(nstep=100, bootstrap=True, nboot=50000, nproc=pargs.nproc,
                        outfile=outfile)


if __name__ == '__main__':
    args = parser()
    main(args)
