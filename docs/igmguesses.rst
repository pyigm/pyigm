.. highlight:: rest

.. _IGMGuesses:

**************
IGMGuesses GUI
**************

.. index:: IGMGuesses

[NT and JB will fill in here.]

Overview
========

IGMGuesses is a GUI for an straightforward identification of
absorption features in the spectra of background sources (e.g. QSOs,
GRBs, stars, galaxies, etc). IGMGuesses is not meant to provide
an optimal Voigt profile fitting of absorption features, but to
aid with the identification of absorption lines and provide reasonable
first guesses on the main parameters of absorption profiles (z, N, b);
these initial guesses should be then used for subsequent Voigt profile
fitting. Still, IGMGuesses does provide an easy handling of blends,
which can become a problem for high-redshift or high-S/N data.

Usage
=====

IGMGuesses should be called from terminal using the script `pyigm_igmguesses`::

    >pyigm_igmguesses in_file [-h] [-o OUT_FILE] [-p PREVIOUS_FILE] [--fwhm FWHM]
                        [--n_max_tuple N_MAX_TUPLE] [--min_strength MIN_STRENGTH]
                        [--min_ew MIN_EW]

Where `in_file` is the input spectrum. This file is usually a fits table, with
columns for wavelength, flux and flux uncertainty. IGMGuesses works on a
continuum-normalized spectrum, thus it has the extra requirement that the
`in_file` also contains an estimation of the continuum. An easy way to make
sure that `in_file` has a compatible format is by using the output file of
linetools' `lt_continuumfit` script (see https://github.com/linetools/linetools).

Optional arguments
++++++++++++++++++

    =============================================== ================================================================================
    Argument                                        Description
    =============================================== ================================================================================
    -h, --help                                      Show the help message and exit
    -o OUT_FILE, --out_file OUT_FILE                Output JSON file with absorption model [Default is called IGM_model.json]
    -p PREVIOUS_FILE, --previous_file PREVIOUS_FILE Input JSON file with absorption model [Previously generated]
    --fwhm FWHM                                     FWHM Gaussian smoothing for fitting (in pixels) [Default is 3]
    --n_max_tuple N_MAX_TUPLE                       Maximum number of transitions per ion species to display
    --min_strength MIN_STRENGTH                     Minimum strength for transitions to be displayed; choose values between (0,14.7)
    --min_ew MIN_EW                                 Minimum EW (in AA) for transitions to be shown/stored within a component.
                                                    This is useful to get rid of very weak transitions from the model
    =============================================== ================================================================================

The number of transitions displayed by default can be excessive for many
(especially low-redshift) spectra, so using the argument "--n_max_tuple 5"
can be a decent starting option.

Component definition
====================
We remind the user that IGMGuesses works on a "absorption component"
basis as given by linetools. Thus, absorption features seen in a spectrum
are the result of the superposition of single or multiple absorption
components. An absorption component is described by the tuple (ion, N, b, z),
where ion is the ion species (e.g. HI, CIV, CIII, SiII), N is the column density,
b is the Doppler parameter, z is the redshift. Be aware that this definition may be
different from the convention used in other softwares.

In order for IGMGuesses to estimate (N, b, z) one may also need to specify a
rest-frame velocity window associated to the component (but note that more generally
one could use different velocity windows for individual absorption lines in a
given component).


Line reliability
================

Because the identification of absorption lines in a given spectrum
is not always 100% certain, in IGMGuesses we have incorporated three
levels of reliability for a component identification, defined as follows.

Certain (a): These include ion components with multiple
transitions where at least two of them are available and visible
in the spectrum, showing the expected ratios and kinematic
structure between them. Common absorption seen at `z=0` can
also be considered reliable.

Possible (b): These can include ion components from single
(or multiple) transition ions, or those

Uncertain (c): These corresponds to those


Line identification algorithm
=============================

Although the users are free to use IGMGuesses as they please,
here we provide a simple algorithm for a systematic absorption
line identification which has empirically proven to be very
efficient.


1. Identify all absorption features available at redshift z = 0,
and assign them to the corresponding reliability category (see above).
