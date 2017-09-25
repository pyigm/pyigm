.. highlight:: rest

***************
DLASurvey Class
***************

Notebooks
=========

.. toctree::
   :maxdepth: 1

   DLA <DLASurvey_examples>

Overview
========

This Class, a child of :ref:`igmsurvey`, is designed to organize
and analyze a survey of DLA systems
(using DLASystem objects).

In general, a DLASurvey is a unique collection of
DLASystem objects.  It is specified by the number of
systems and the references.


Instantiation
=============

The DLASystem Class may be instantiated in a few ways.
The default sets the properties listed above::

    dlas = DLASurvey(ref='null')


One may instantiate and then add one or more DLASystem objects::

    coord = SkyCoord(ra=123.1143, dec=-12.4321, unit='deg')
    dlasys = DLASystem(coord, 1.244, [-300,300.]*u.km/u.s, 20.4)
    dlasys.name = 'Sys1'
    #
    coord2 = SkyCoord(ra=223.1143, dec=42.4321, unit='deg')
    dlasys2 = DLASystem(coord2, 1.744, [-300,300.]*u.km/u.s, 21.7)
    dlasys2.name = 'Sys2'
    # Add systems
    dlas.add_abs_sys(dlasys)
    dlas.add_abs_sys(dlasys2)

Attributes/Properties
=====================

========   ============== ============================================
Variable   Type           Description
========   ============== ============================================
nsys       int            Number of systems in the survey
ref        str            References for the survey
========   ============== ============================================


Datasets
========

We are striving to include all of the major DLA
survey data published by the community.

Here is a Table describing the various samples that may
be accessed.

.. _Peroux03: http://adsabs.harvard.edu/abs/2003MNRAS.346.1103P
.. _PW09: http://adsabs.harvard.edu/abs/2009ApJ...696.1543P
.. _G09: http://adsabs.harvard.edu/abs/2009A%26A...508..133G
.. _Neeleman+13: http://adsabs.harvard.edu/abs/2013ApJ...769...54N
.. _Crighton+15: http://adsabs.harvard.edu/abs/2015MNRAS.452..217C
.. _Neeleman+16: http://adsabs.harvard.edu/abs/2016ApJ...818..113N
.. _Sanchez+16: http://adsabs.harvard.edu/abs/2016MNRAS.456.4488S

========== =============================  =================== =====================================
Survey     Call                           Reference(s)              Description
========== =============================  =================== =====================================
P03        DLASurvey.load_P03()           `Peroux03`_         DLAs discovered from APM quasars
SDSS_DR5   DLASurvey.load_SDSS_DR5()      `PW09`_             DR5
G09        DLASurvey.load_G09()           `G09`_              ESI survey by Guimaraes
H100       DLASurvey.load_H100()          `Neeleman+13`_      100 unbiased HIRES spectra
GGG        DLASurvey.load_GGG()           `Crighton+15`_      DLAs from the GGG Survey
XQ-100     DLASurvey.load_XQ100()         `Sanchez+16`_       DLAs from the XQ-100 Survey
HST16      DLASurvey.load_HST16()         `Neeleman+16`_      Sample of DLAs drawn from HST spectra
========== =============================  =================== =====================================

Loading
+++++++

Here is an example of loading the H100 dataset::

    h100 = DLASurvey.load_H100()


Metallicities
+++++++++++++

If the DLA systems have measured metallicites
in the dataset, one may access them with the ZH
attribute::

    ZH = h100.ZH

Ionic column densities
++++++++++++++++++++++

Similarly, some datasets include ionic measurements.
These are loaded into each system but a Table of
measurements may be generated::

    SiII_clms = h100.ions((14,2))  # SiII

This astropy Table includes name (of the DLA),
flagN, logN, sig_logN, attributes.


Plots
=====

Binned Evaluations
==================

One can calculate standard DLA statistics in
bins of redshift, NHI, etc. with the following
methods.  Each requires that a sightline Table
exists, so that g(z) may be evaluated.

g(z)
++++

Provide the sightlines Table is filled and has keys
Z_START and Z_END, this method will generate a
selection function :math:`g(z)` curve::

   sdss = DLASurvey.load_SDSS_DR5()
   zeval, gz = sdss.calculate_gz()


f(N,X)
++++++

Calculate the NHI frequency distribituion in bins of NHI and z.  e.g., ::

    fN, fN_lo, fN_hi = sdss_stat.binned_fn([20.3, 20.5, 21., 21.5, 22.], [2, 2.5], log=True)

Setting log=True returns log10 values for f(N) and its error.

l(X)
++++

Calculate the incidence per unit dX in binned redshift
intervals.  Default is over all NHI values.  Here is an
example::

    lX, lX_lo, lX_hi = sdss_stat.binned_lox([2., 2.5, 3])

This calculates lX and its error in the intervals z=[2,2.5]
and z=[2.5,3.].

rhoHI
+++++

Similar to the last two methods but now for the HI mass density.
Here is an example::

    zbins = [2., 2.5, 3.]
    NHImnx = (20.3, 23.)
    rho, rho_lo, rho_hi = sdss_stat.binned_rhoHI(zbins, NHImnx)

rho will have units of Solar mass per Mpc^3.

Output
======
