.. highlight:: rest

***************
CGMSurvey Class
***************

.. _cgmsurvey:


.. _Johnson15: http://adsabs.harvard.edu/abs/2015MNRAS.449.3263J
.. _Burchett16: http://adsabs.harvard.edu/abs/2016ApJ...832..124B

CGMSurvey
=========

This Class is designed to organize and analyze a survey of
CGM systems.

By definition, an CGMSurvey is a unique collection of
CGM objects.  It is specified by the number of
systems and the references.

In general, it is used to load up datasets
from the literature, e.g. COS-Halos.

Johnson+15
==========

HI and OVI measurements by `Johnson15`_ for
galaxies at low redshift.  You can load them
up with::

    j15 = CGMAbsSurvey.load_J15()
    # OVI table
    OVI_tbl = j15.ion_tbl((8,6))

Note that we have not duplicated COS-Halos data here.
Only their IMACS and SDSS systems.

Burchett+16
============

HI and CIV absorption measurements for galaxies at z~0
measured by `Burchett16`_.  You can instantiate the
survey easily enough::

    b16 = CGMAbsSurvey.load_B16()
    # Generate a HI Table
    Lya_tbl = b16.ion_tbl((1,1))

