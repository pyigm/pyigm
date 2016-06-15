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
aid with the identification of absorption lines and provide reasonable first guesses
on the main parameters of absorption profiles (z, N, b); these initial guesses should be then used for
subsequent Voigt profile fitting. Still, IGMGuesses should provide an easy
handling of blends, which can become a problem for high-redshift
or high-S/N data.

The input spectrum has to be normalized (e.g.
by using linetools' lt_continuumfit script).

Because the identification of absorption lines is not always
certain, in IGMGuesses we have incorporated three levels of
reliability for a line identification, defined as follows.

Line reliability
----------------

Certain (a): These include ion components with multiple
transitions where at least two of them are available and visible
in the spectrum, showing the expected ratios and kinematic
structure between them. Common absorption seen at `z=0` can also be considered reliable.

Possible (b): These can include ion components from single
(or multiple) transition ions, or those



Although the user is free to use IGMGuesses as their discretion,
here we provide a simple algorithm for a systematic absorption
line identification which has empirically proven to be very
efficient.




Line identification algorithm
-----------------------------

1. Identify all absorption features available at redshift z = 0.