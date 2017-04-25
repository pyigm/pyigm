.. highlight:: rest

.. _IGMGuesses:

**************
IGMGuesses GUI
**************

.. index:: IGMGuesses

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

    ============================================== =============================================== ==============
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
    --vlim VLIM                                     Velocity limit (in km/s) for the display        500.
    --external_model                                Name of an external spectrum model fits file
    =============================================== =============================================== ==============

The number of transitions available for some ions  can be excessive for many,
especially low-redshift spectra (e.g. CI, CI**), so using the default argument
"--n_max_tuple 5" is a decent starting option. Feel free however to try different
values depending on your scientific needs. If an external model is given, you can toggle
displaying/hiding it using the keystroke 'E' (for external).


Component definition
====================

We remind the user that IGMGuesses works on a "absorption component"
basis as given by linetools AbsComponent object. Thus, absorption features
seen in a spectrum are the result of the superposition of single or multiple absorption
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

- **Certain (flag a)**: These include components with multiple
  transitions where at least two of them are available and visible
  in the spectrum, and showing the expected ratios and kinematic
  structure. Common absorption seen at `z=0` fall in this category,
  as well as strong HI showing multiple Lyman series transitions.
  Use 'P' to toggle on/off colorful display of components and
  this will appear in green.

- **Probable (flag b)**: These include components from single
  transition ions that are at the same redshift (within a reasonable
  velocity window) from another certain component (e.g. CIII at the
  same redshift than a certain HI). Another case where this category
  should apply is when we have components from ions with multiple
  transitions but that for some reason only 1 transition is clearly seen
  (e.g. due to heavy blends, poor S/N, wavelength coverage, etc). Examples of these
  could be weak HI where only HI Lya is visible, or a OVI component where one of
  the transition is blended with something else thus not certain. Use 'P' to
  toggle on/off colorful display of components and this will appear in blue.

- **Uncertain (flag c)**: These correspond to those components that
  based on the user experience are likely to be an incorrect identification.
  Hopefully components identified in this category will be later replaced by a
  better identification. These could include an unphysical narrow line, artifacts,
  etc. Use 'P' to toggle on/off colorful display of components and
  this will appear in red.

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
     structures (e.g. `dv = -200 km/s` for sightlines close to Andromeda galaxy).

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
     step (iv).


Basic IGMGuesses usage
======================

Panels description
------------------

Once IGMGuesses is launched from terminal, a GUI will appear with four
main panels, these are:

- **Velocity Windows** *(left)*: This is the main graphic display where different
  transitions of different ion are shown at a given redshift (see top).
  Transitions of the same ion are grouped by the same color
  sorted by strength. The model is plotted on top of the spectrum as a brown line.
  The spectrum +- 1 sigma uncertainty is plotted around the zero level in arbitrary
  units (re-scaled for clarity). Residuals are plotted around the zero level
  in the same units as the uncertainty, which helps to assess the statistical
  significance of absorption/emission features. You can add/remove columns or rows
  by using the keystrokes 'C'/'c' or 'K'/'k', respectively. You can also go to the
  previous/next "page" using the keystrokes '-'/'='. Other navigation/display
  options are available (by pressing '?' within IGMGuesses the full
  option list will be generated in the terminal)

- **Component Widget** *(top right)*: This is the widget that displays and
  controls the main parameter for modeling the currently selected absorption
  component, i.e. (N, b, z). You can slightly increase/decrease the current
  values using the keystrokes 'N'/'n', 'V'/'v' and '<'/'>', respectively. You
  can also modify them directly by hand just editing the respective values on
  the widget itself. After done with editing, make sure you press the "Update" button.
  Reliability flags must be chosen from the available list (see
  above for definitions). A string comment can also be entered.

- **Components List** *(middle right)*: This widget displays and allows to
  navigate between currently defined components. The name of the component
  includes the redshift and ion species. By selecting a component from this list,
  its parameters can be modified from the *Component Widget* but in order to move
  to the corresponding redshift one needs to manually navigate to it (by using the
  keystroke `^` or the Space Bar; see below)

- **Line List** *(bottom right)*: This widget displays the current parent LineList
  (see linetools LineList Class for further details)
  where ion transitions are selected from. By selecting/unselecting them you can
  control which transitions are displayed in the *Velocity Windows*. Built-in LineList
  can be loaded by keystrokes 'H' (HI Lyman series), 'T' (Common Strong IGM transitions),
  'F' (Full list of ISM known transitions). Of course, depending on redshift some
  transitions may or may not be available in the current spectrum; in order to select
  those available at the current redshift of interest from the current LineList
  you can use keystroke 'U' (update). By doing 'U' it will also restrict the
  subset of lines satisfying  `n_max_tuple` and `min_strength` as given in the
  initialization (see Table above).


Adding/Removing/Selecting/Editing Components
--------------------------------------------

- **Adding a component**: Click on the *Velocity Window* associated to the
  relevant ion transition and define the rest-frame velocity limits by pressing
  'A' twice (one for each limit). IGMGuesses will use the pixel information within
  those limits to fit a Voigt profile (convolved with a Gaussian of FWHM as given
  in the initialization; see above) in order to estimate the (N,b,z) parameters.
  It only uses the pixels in the selected transition of the given ion for
  guessing the model parameters. Thus, it is recommended to select a transition
  that is not saturated, blended or in a poor S/N spectral region, when possible.
  Still, a model of the component will be displayed encompassing *all* the transitions.
  Another convenient way to add a component is by using the single keystroke 'a',
  which tells IGMGuesses to use the same velocity window as the previously selected
  component.

  - **Removing a component**: Use 'D' to remove the *closest* (in velocity) component to
  the cursor position in a given velocity window panel. You can also remove a component
  regardless of whether is being displayed in the velocity window panels by selecting it
  from the *Components List* widget and then pressing the keystroke 'd'.

  - **Selecting a component**: Use 'S' to select a component from the velocity window
  panel at the cursor position (you may need to click on the white area of the velocity
  panel to make sure IGMGuesses recognizes you are selecting from the velocity panel).
  You can also select a component directly from the *Components List* widget.

  - **Editing a component**: Use the keystrokes N'/'n', 'V'/'v' and '<'/'>', to slightly
  increase/decrease the values of (N, b, z), respectively, for the currently selected
  component. Alternatively, you can modify these values by hand directly from the
  *Component Widget*. The lower/upper velocity limits of a given component can be
  modified by using the keystrokes '1'/'2', respectively. Use 'R' if you wish to refit the
  data using the current velocity limits.


Redshift navigation
-------------------

There are two ways that you can navigate in redshift space:

    - **Panning**: Use the Space Bar to set the redshift given by the the rest-frame velocity
    of the cursor.

    - **By hand**: Use '^' to select the redshift by hand. A pop-up widget should appear.


Model visualization
-------------------

The superposition of all absorption components identified provides the overall model
for the given spectrum. Such model is by default displayed as a brown line.
Other options for model visualization include:

    - **Colorful display**: use 'P' to toggle on/off colorful display in which each
    component of the model is plotted by different colors depending on their assigned
    reliability (see above).

    - **Show/hide labels**: use 'L' to toggle on/off the label of each identified component.
    The label is composed by the ion transition, its redshift, and its reliability flag
    (see above)

    - **Show/hide model**: use 'M' to toggle on/off the model being displayed.

    - **Show/hide external model**: use 'E' to toggle on/off the model being displayed.