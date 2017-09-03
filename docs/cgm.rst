.. highlight:: rest

***********
CGM Classes
***********

.. index:: CGM

Notebooks
=========

.. toctree::
   :maxdepth: 1

   CGM Examples <CGM_examples>
   CGM P11 <CGM_P11>
   CGMModel <https://github.com/pyigm/pyigm/blob/master/docs/examples/CGM_Models.ipynb>
   COS-Halos <https://github.com/pyigm/pyigm/blob/master/docs/examples/COS_Halos_Examples.ipynb>


Overview
========

There are a series of classes designed for CGM analyses.
We describe each in turn.

.. _cgm-class:

CGM
===

This class is intended to provide a physical representation
of the CGM.  A Galaxy object is required for instantiation.
The key attributes are tabulated below.  This includes several
which may be controversial, e.g. the physical extent of the CGM,
and a description of the baryons in `phases'.
This is the least developed of the classes.

=========   ============== ============================================
Attribute   Type           Description
=========   ============== ============================================
galaxy      Galaxy object  Describes the galaxy hosting the CGM
rlim        Quantity       Physical extent of the CGM (e.g. 300 kpc)
phases      dict           Intended to organize the properties of various CGM phases, e.g. mass, metallicity
cgm_abs     list           List of CGMAbsSys classes
=========   ============== ============================================

Instantiation
-------------

The only requirement is a Galaxy object.::

   radec = (125*u.deg, 45.2*u.deg)
   gal = Galaxy(radec,z=0.3)
   cgm = CGM(gal)

CGMAbsSys
=========

This class enables an absorption-line analysis of a CGM.  This has been
the primary approach to CGM analysis to date.  This class requires
both a Galaxy object and :doc:`igmsys` object for instantiation.

=========   ================= ============================================
Attribute   Type              Description
=========   ================= ============================================
galaxy      Galaxy object     Describes the galaxy hosting the CGM
igm_sys     IGMSystem         IGM system object to describe absorption on the line-of-sight
cosmo       astropy.cosmology Cosmological model; defaults to WMAP9
rho         Quantity          Impact parameter from galaxy to sightline
=========   ================= ============================================

Instantiation
-------------

One must invoke CGMAbsSys with both a Galaxy and IGMSystem object.::

   radec_qso = (125*u.deg, 45.203*u.deg)
   igmsys = IGMSystem('CGM', radec_qso, gal.z, [-500,500]*u.km/u.s)
   #
   radec = (125*u.deg, 45.2*u.deg)
   gal = Galaxy(radec,z=0.3)
   #
   cgmabs = CGMAbsSys(gal,igmsys)

CGMAbsSurvey
============

This class organizes a survey of CGMAbsSys objects, e.g. COS-Halos.

=========   ================= ============================================
Attribute   Type              Description
=========   ================= ============================================
survey      str               Name of the survey, e.g. COS-Halos
ref         str               References for the survey
cgm_abs     list              list of CGMAbsSys objects
mask        bool array        Mask
=========   ================= ============================================

Instantiation
-------------

This class requires no input for instantiation.  But, it is expected
that one will fill the cgm_abs list with CGMAbsSys objects.

Properties and Methods
----------------------

nsys
++++

This property returns the number of CGMAbsSys objects in the survey (ignores mask).::

   nsys = cgmsurvey.nsys

getattr
+++++++

This is overloaded to return an array of properties from one of the
internal sets of objects in the survey.  The order of attribution is
CGMAbsSys objects, Galaxy objects, and then IGMSystem objects. ::

   rho_array = cgmsurvey.rho # Grabs rho from CGMAbsSys objects
   z_array = cgmsurvey.z  # Grabs galaxy redshifts
   coord = cgmsurvey.coord # Grabs galaxy coordinates

Miscelaneous
------------

There are a few methods related to CGM analysis available.

dN/dX
+++++

Calculate dN/dX given a cosmology and a paremterization of the
halos.  Code of interest is cgm.analysis.dndx_rvir::

   dNdX = cgm.analysis.dndx_rvir()


COS-Halos
---------

All of the measurements related to the COS-Halos survey
`Werk et al. (2011) <http://adsabs.harvard.edu/abs/2012ApJS..198....3W>`_,
`Tumlinson et al. (2013) <http://adsabs.harvard.edu/abs/2013ApJ...777...59T>`_,
`Werk et al. (2013) <http://adsabs.harvard.edu/abs/2013ApJS..204...17W>`_
`Prochaska et al. (2017) <http://adsabs.harvard.edu/abs/2017ApJ...837..169P>`_
are included in the data directory.  You can read in the dataset with::

   from pyigm.cgm import cos_halos as pch
   coshalos = pch.COSHalos()

See the Notebook in examples for further details.


CGMModel
========

Simple models designed to reproduce observations of the CGM
may be generated using the CGMModel class.  Current implementation
is based on the Mathews & Prochaska (2017) halo profile.

See the CGM_Models Notebook for examples.
