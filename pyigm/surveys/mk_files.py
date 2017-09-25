""" Script to generate output files from IGM Surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Execute with:  python mk_files -h


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Create one or more files for IGM Surveys.')
    parser.add_argument("--dla_lz", help="Fit lz with atan model and create dla_lz_boot file", action="store_true")
    parser.add_argument("--dla_dpow", help="Fit f(N) of DLA with double power law", action="store_true")
    parser.add_argument("--dla_nenH", help="Add ne/nH into the fit file", action="store_true")
    parser.add_argument("--all", help="Create dla_lz_boot file", action="store_true")
    parser.add_argument("--nproc", type=int, default=4, help="Number of processors")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):


    from pkg_resources import resource_filename

    from linetools import utils as ltu

    from pyigm.surveys.analysis import fit_atan_dla_lz, fit_fN_dblpow
    from pyigm.surveys.dlasurvey import load_dla_surveys, update_dla_fits
    from pyigm.surveys.dlasurvey import DLASurvey


    pargs = parser()

    # DLA l(z) analysis [Prochaska & Neeleman 2017]
    if pargs.dla_lz or pargs.all:
        surveys = load_dla_surveys()
        boot_file = resource_filename('pyigm', 'data/DLA/dla_lz_boot.fits.gz')
        dfits, _ = fit_atan_dla_lz(surveys, nstep=100, bootstrap=False, nboot=50000, nproc=pargs.nproc,
                        boot_out=boot_file)
        # Write
        dfits['lz']['atan']['Ref'] = 'Prochaska & Neeleman 2017'
        update_dla_fits(dfits)

    # Fit double power law to f(N) of DLA [PW09 only]
    if pargs.dla_dpow or pargs.all:
        sdss_dr5 = DLASurvey.load_SDSS_DR5()
        dfits, best, Ndgrid, a3grid, a4grid, lik = fit_fN_dblpow(
            sdss_dr5.NHI, (-3., -1.1), (-6,-2), (21., 22.), nstep=100)
        # Write
        dfits['fN']['dpow']['Ref'] = 'PHW05'
        update_dla_fits(dfits)

    # DLA ne/nH
    if pargs.dla_nenH or pargs.all:
        dfits = {}
        dfits['nenH'] = {}
        dfits['nenH']['loglin'] = dict(bp=-2.592, m=-0.64)
        dfits['nenH']['loglin']['Ref'] = 'Neeleman+15; PN17'
        # Update
        update_dla_fits(dfits)

if __name__ == '__main__':
    args = parser()
    main(args)
