.. highlight:: rest

.. _fitlls:

***********
XFitLLS GUI
***********

.. index:: IGMGuesses

Overview
========

XFitLLSGUI is a GUI for the by-hand fit to one or more LLSs
in an input spectrum.  One can modify both the continuum
and LLS properties and write the solution to disk as a JSON file.

Running XFitLLSGUI
==================

XFitLLSGUI should be called from the terminal using the
script `pyigm_fitlls`.  Here is the usage::

    wolverine> pyigm_fitlls -h
    usage: pyigm_fitlls [-h] [--out_file OUT_FILE] [--smooth SMOOTH]
                        [--lls_fit_file LLS_FIT_FILE]
                        in_file zqso

    Parser for FitLLSGUI (v1.0)

    positional arguments:
      in_file               Spectral file
      zqso                  Use Telfer template with zqso

    optional arguments:
      -h, --help            show this help message and exit
      --out_file OUT_FILE   Output LLS Fit file
      --smooth SMOOTH       Smoothing (pixels)
      --lls_fit_file LLS_FIT_FILE
                            Input LLS Fit file


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
--lls_fit_file LLS_FIT_FILE                    Input .JSON file containing a previous DLA fit
============================================== =============================================== ==============


Usage
=====

When the LLS gui launches, you should begin by left clicking
in the plot window.  This will activate plotting and key strokes
(see below for the help message). Then press "#" to specify the model spec and to add a continuum.
You should then refine the continuum ("C", "1", "2").

From there, navigate to the LLS of interest, add a model ("A" or "F"),
refine it and the continuum, etc.

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
    C         : Set the continuum normalization to the cursor
    1,2       : Tilt the continuum
    A         : Add a new LLS at cursor
    F         : Add a new LLS automagically near the cursor (probably)
    g         : Move nearest Lyman line to cursor and reset z
    N/n       : Increase/decrease NHI
    V/v       : Increase/decrease bvalue
    Z/z       : Increase/decrease zabs
    D         : Delete LLS
    $         : Toggle displaying metal lines
    6,7,8,9   : Add forest lines
    ?         : Print these help notes
    Q         : Quit the GUI

