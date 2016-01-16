.. highlight:: rest

***************
IGMSurvey Class
***************

.. index:: AbsSystem

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <IGMSurvey_examples>
   LLS <LLSSurvey_examples>
   DLA <DLASystem_examples>

Overview
========

This Class is designed to organize and analyze a survey of
absorption systems (defined as AbsSystem objects).

By definition, an IGMSurvey is a unique collection of
AbsSystem objects.  It is specified by the number of
systems and the references.


Instantiation
=============

The AbsSystem Class may be instantiated in a few ways.
The default sets the properties listed above::

    gensurvey = GenericIGMSurvey()


More commonly, one will instantiate with one or more IGMSystem objects::

    coord = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys = IGMSystem('MgII', coord, 1.244, [-300,300.]*u.km/u.s, NHI=16.)
    gensys.name = 'Sys1'
    #
    coord2 = SkyCoord(ra=223.1143*u.deg, dec=42.4321*u.deg)
    gensys2 = IGMSystem('MgII', coord2, 1.744, [-300,300.]*u.km/u.s, NHI=17.)
    gensys2.name = 'Sys2'
    #
    gensurvey.add_abs_sys(gensys1)
    gensurvey.add_abs_sys(gensys2)


Attributes/Properties
=====================

========   ============== ============================================
Variable   Type           Description
========   ============== ============================================
nsys       int            Number of systems in the survey
ref        str            References for the survey
========   ============== ============================================

Sub Classes
===========

LLS
+++

Subclass for an LLS survey.  There are many published surveys
that can be read in.  Several require access to the Internet
which will then generate a file on your drive for future use.
There is also a method to handle the .dat and .lst files used
by JXP.   See :doc:`LLSSurvey_examples` for more.

Here is a Table describing the various samples that may
be accessed.

.. _Prochaska+10: http://adsabs.harvard.edu/abs/2010ApJ...718..392P
.. _Omeara+11: http://adsabs.harvard.edu/abs/2011ApJS..195...16O
.. _Fumagalli+13: http://adsabs.harvard.edu/abs/2013ApJ...775...78F
.. _Prochaska+15: http://adsabs.harvard.edu/abs/2015ApJS..221....2P
.. _Fumagalli+16: http://adsabs.harvard.edu/abs/2016MNRAS.455.4100F

========== =============================  =================== ================================
Survey     Call                           Reference(s)              Description
========== =============================  =================== ================================
SDSS       LLSSurvey.load_SDSS_DR7()      `Prochaska+10`_     tau>2 survey
z2_HST     LLSSurvey.load_HST_ACS()       `Omeara+11`_        tau>2 with HST/ACS
           LLSSurvey.load_HST_WFC3()      `Omeara+11`_        tau>1 with HST/WFC3
z3_MagE    LLSSurvey.load_mage_z3()       `Fumagalli+13`_     tau>2 with Magellan/MagE
HD-LLS     LLSSurvey.load_HDLLS()         `Prochaska+15`_     Ionic column densities
                                          `Fumagalli+16`_     and metallicity PDF
Literature lls_literature.load_lls_lit()  See `Fumagalli+16`_ Ionic column densities
========== =============================  =================== ================================

Below are additional options for a few of these.

HD-LLS DR1
----------

The standard call loads the ionic column densities and
metallicity PDFs.
One call access the spectra with::

   hdlls = LLSSurvey.load_HDLLS(grab_spectra=True)

This will grab 154Mb of data from the internet, and place
them within pyigm/data/LLS/HD-LLS.


DLA
+++

Subclass for DLA survey.  Presently handles the .dat and .lst files used
by JXP.   See :doc:`DLASurvey_examples` for more.

Here is a Table describing the various samples that may
be accessed.

.. _PW09: http://adsabs.harvard.edu/abs/2009ApJ...696.1543P
.. _Neeleman+13: http://adsabs.harvard.edu/abs/2013ApJ...769...54N

========== =============================  =================== ================================
Survey     Call                           Reference(s)              Description
========== =============================  =================== ================================
SDSS_DR5   DLASurvey.load_SDSS_DR5()      `PW09`_             DR5
H100       DLASurvey.load_H100()          `Neeleman+13`_      100 unbiased HIRES spectra
========== =============================  =================== ================================

Plots
=====

Methods
=======

g(z)
++++

Provide the sightlines Table is filled and has keys
Z_START and Z_END, this method will generate a
selection function :math:`g(z)` curve::

   # LLS
   z3mage = LLSSurvey.load_mage_z3()
   zeval, gz = z3mage.calculate_gz()
   # DLA
   sdss = DLASurvey.load_SDSS_DR5()
   zeval, gz = sdss.calculate_gz()

Output
======
