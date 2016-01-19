.. highlight:: rest

*********
Continuum
*********

.. index:: contniuum

Notebooks
=========

.. toctree::
   :maxdepth: 1

Overview
========

Observational analysis of the IGM generally requires the
continuum normalization of the background source observed.
pyigm provides a few methods useful for such analysis.


Quasars
=======

It is common practice to adopt a quasar template as a first
guess for the continuum of these sources.  pyigm provides
two methods.

Telfer
------

One may generate a quasar continuum using the Telfer+05 template
(radio-quiet). ::

   from pyigm.continuum import quasar as pyicq
   telfer = pyicq.get_telfer_spec(3.)

The resultant object is a XSpectrum1D object of the
model continuum.  One may specify whether to include an
average IGM opacity model (see :doc:`fn` for a description
of the fN model).  Multiprocessing is employed for speed. ::

   telfer = pyicq.get_telfer_spec(3., igm=True, nproc=4)


WFC3
----

Individual quasars without apparent, significant Lyman limit absorption
that were observed with the HST/WFC3 instrument may provide *real*
continua.  These may be selected randomly from the full set or by
index. ::

   wfc3, _ = pyicq.wfc3_continuum(wfc3_indx=0, zqso=2.)
   wfc3, idx = pyicq.wfc3_continuum(zqso=2.5)

Again, the returned object is a XSpectrum1D object.