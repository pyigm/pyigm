.. highlight:: rest

**********
CUBA Class
**********

.. index:: CUBA

Notebooks
=========

.. toctree::
   :maxdepth: 1

   CUBA Examples <CUBA_examples>

Overview
========

This Class is designed to interface with the somewhat
cumbersome output of CUBA.  It may then be used to
perform EUVB calculations.


Instantiation
=============

The CUBA Class is instantiated with a CUBA file.
If one is not provided, a default file from within
the package is loaded.::

	from pyigm.euvb import cuba as pycuba
	cuba = pycuba.CUBA()



Attributes/Properties
=====================

========   ============== ============================================
Variable   Type           Description
========   ============== ============================================
z          ndarray        Redshifts of the evaluation
wave       Quantity array Wavelengths of the evaluation
Jnu        Quantity array Intensity of the EUVB
========   ============== ============================================


Plots
=====

One can generate a quick plot at a given redshift with::

   cuba.plot(z)

Methods
=======

zinterp
-------

The intensity may be evaluated at any redshift with this method::

   jnu = cuba.zinterp_jnu(z)

phi
---

The ionizing photon flux may be calculated with this method, again
at a given redshift::

   phi = cuba.phi(z)

Output
======

A future method will be provided to generate a file that may
be used with `Cloudy <http://www.nublado.org/>`_.
