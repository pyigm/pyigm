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

Running IGMGuesses
==================

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

    =============================================== =============================================== ==============
    Argument                                        Description                                     Default
    =============================================== =============================================== ==============
    -h, --help                                      Show the help message and exit
    -o OUT_FILE, --out_file OUT_FILE                Output JSON file with absorption model          IGM_model.json
    -p PREVIOUS_FILE, --previous_file PREVIOUS_FILE Input JSON file with absorption model (if any)
    --fwhm FWHM                                     FWHM Gaussian smoothing for fitting (in pixels) 3
    --n_max_tuple N_MAX_TUPLE                       Maximum number of transitions per ion species   5
                                                    to display
    --min_strength MIN_STRENGTH                     Minimum strength for transitions to be          0
                                                    displayed; choose values between (0,14.7)
    --min_ew MIN_EW                                 Minimum EW (in AA) for transitions to be
                                                    shown/stored within a component. This is useful 0.005
                                                    for getting rid of very weak transitions from
                                                    the model
    =============================================== =============================================== ==============

The number of transitions available for some ions  can be excessive for many,
especially low-redshift spectra (e.g. CI, CI**), so using the default argument
"--n_max_tuple 5" is a decent starting option. Feel free however to try different
values depending on your scientific needs.


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
given component). This velocity window can also be used to figure out which components are blended
or unblended, in the expected subsequent Voigt profile fitting.


Line reliability
================

Because the identification of absorption lines in a given spectrum
is not always 100% certain, in IGMGuesses we have incorporated three
levels of reliability for a component identification, defined as follows.

- **Certain (label a)**: These include components with multiple
  transitions where at least two of them are available and visible
  in the spectrum, and showing the expected ratios and kinematic
  structure. Common absorption seen at `z=0` fall in this category,
  as well as strong HI showing multiple Lyman series transitions.

- **Possible (label b)**: These include components from single
  transition ions that are at the same redshift (within a reasonable
  velocity window) from another certain component (e.g. CIII at the
  same redshift than a certain HI). Another case where this category
  should apply is when we have components from ions with multiple
  transitions but that for some reason only 1 transition is clearly seen
  (e.g. due to heavy blends, poor S/N, wavelength coverage, etc). Examples of these
  could be weak HI where only HI Lya is visible, or a OVI component where one of
  the transition is blended with something else thus not certain.

- **Uncertain (label c)**: These correspond to those components that
  based on the user experience are likely to be an incorrect identification.

- **Unknown (not implemented yet)**: This category is for those absorption
  features that cannot be explained with current information.


Line identification algorithm
=============================

Although the users are free to use IGMGuesses as they please,
here we provide a simple algorithm for a systematic absorption
line identification which has empirically proven to be very
efficient.

- 1. Identify all absorption components available at redshift z = 0,
and assign them to the corresponding reliability category (see above).
Depending on the (RA, Dec) of the QSO also inspect dv close to known
structures (e.g. dv = -200 km/s for sightlines close to Andromeda galaxy).

- 2. Identify all absorption components available at redshift z = z_qso,
and assign them to the corresponding reliability category (see above).

- 3. Identify HI components showing at least two transitions (e.g. Ly-alpha
and Ly-beta, Ly-beta and Ly-gamma, etc), starting at z=z_qso until z=0, and
assign them to the 'certain' category. This classification includes the full
Lyman series transitions of the component available in the spectrum.

- 4. Identify all possible metal absorption components within a reasonable
rest-frame velocity window (dv) from each HI redshift found in the previous
step and assign them to the corresponding reliability category (see above).

- 5. Assume all the unidentified absorption features to be HI Lya starting from
z=z_qso down to z=0, and assign them to the 'possible' category. Then repeat
step 4.


Basic IGMGuesses usage
======================

Once IGMGuesses is launched from terminal, a GUI will appear with four
main panels, these are:

1. Velocity Windows: This is the main graphic display where different
transitions
