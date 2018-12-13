.. highlight:: rest

.. _fitdla:

***********
XFitDLA GUI
***********

.. index:: IGMGuesses

Overview
========

XFitDLAGUI is a GUI for the by-hand fit to one or more DLAs
in an input spectrum.  One can modify both the continuum
and DLA properties and write the solution to disk as a JSON file.

Running XFitDLAGUI
==================

XFitDLAGUI should be called from the terminal using the
script `pyigm_fitdla`.  Here is the usage::

    usage: pyigm_fitdla [-h] [--out_file OUT_FILE] [--smooth SMOOTH]
                    [--dla_fit_file DLA_FIT_FILE] [--conti_file CONTI_FILE]
                    [--zdla ZDLA] [--NHI NHI]
                    in_file zqso

    Parser for FitDLAGUI (v1.0)

    positional arguments:
      in_file               Spectral file
      zqso                  Use QSO template with zqso

    optional arguments:
      -h, --help            show this help message and exit
      --out_file OUT_FILE   Output DLA Fit file
      --smooth SMOOTH       Smoothing (pixels)
      --dla_fit_file DLA_FIT_FILE
                            Input DLA Fit file
      --conti_file CONTI_FILE
                            Input continuum spectrum
                            Option 'from_spectrum': continuum data will be used from the original (in_file) spectrum.
      --zdla ZDLA           Input DLA redshift
      --NHI NHI             Input DLA NHI


There are 2 required arguments:
(1) `in_file` which is the input spectrum file;
(2) `zqso` which is an estimate of the quasar (or other source)
emission redshift.

Optional arguments
++++++++++++++++++

============================================== =============================================== ==============
Argument                                       Description                                     Default
============================================== =============================================== ==============
-h, --help                                     Show the help message and exit
-o OUT_FILE, --out_file OUT_FILE               Output JSON file with model                     DLA_fit.json
--smooth SMOOTH                                FWHM Gaussian smoothing for fitting (in pixels) 3
--dla_fit_file DLA_FIT_FILE                    Input .JSON file containing a previous DLA fit
--conti_file CONTI_FILE                        Input continuum spectrum file
--zdla                                         Input redshift for an initial DLA model
--NHI                                          Input NHI value for an initial DLA model
============================================== =============================================== ==============


Usage
=====

When the DLA gui launches, you should begin by left clicking
in the plot window.  This will activate plotting and key strokes
(see below for the help message).

From there, navigate to the DLA of interest, add a model, refine
it and the continuum, etc.

The Write button saves the current model to disk as a JSON file.

Keystrokes
==========

Here is a brief description of the key strokes that control
the DLA GUI (also displayed when launching the GUI)::

    i,o       : zoom in/out x limits
    I,O       : zoom in/out x limits (larger re-scale)
    Y         : zoom out y limits
    y         : guess y limits
    W         : Show original zooom
    t,b       : set y top/bottom limit
    l,r       : set left/right x limit
    a,m,d     : Add/modify/delete continuum knot
    A         : Add a new DLA
    g         : Move nearest Lyman line to cursor and reset z
    N/n       : Increase/decrease NHI
    V/v       : Increase/decrease bvalue
    Z/z       : Increase/decrease zabs
    D         : Delete DLA
    $         : Toggle displaying metal lines
    6,7,8,9   : Add forest lines
    ?         : Print these help notes
    Q         : Quit the GUI

