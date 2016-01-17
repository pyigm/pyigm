.. highlight:: rest

********************
IGMGalaxyField Class
********************

.. index:: igmfield

Notebooks
=========

.. toctree::
   :maxdepth: 1

   IGMGalaxyField <IGMGalaxyField_examples>

Overview
========

This Class is designed to enable analysis between galaxies and
absorbers in a chosen field.  In particular, cross-correlation
analysis.  Although, in many ways this is the starting point for
essentially any analysis of galaxies with the IGM.


Instantiation
=============

Instantiation only requires the field coordinates, presumably near the center.::

    field = ('PG1407+265',212.349634*u.deg,26.3058650*u.deg)
    lfield = igmf.IgmGalaxyField((field[1],field[2]), name=field[0], verbose=True)



Attributes/Properties
=====================

=========  ================= ============================================
Variable   Type              Description
=========  ================= ============================================
name       str               Given name for the field
cosmo      astropy.cosmology Given name for the field
igm        ??                A means of conveniently storing IGM system info
targets    astropy.Table     Table of target info
galaxies   astropy.Table     Table of galaxy info
observing  astropy.Table     Table of info on observing the galaxies
selection  ??                Object to enable calculation of the galaxy selection function
=========  ================= ============================================

Methods
=======

Impact Parameter
----------------

Calculate :math:`\rho`, the projected impact parameter from
a given object at a given redshift to the line-of-sight (LOS) coordinate. By defualt, the
projected impact parameter is calculated in co-moving coordinates and the LOS is assumed to be the
field coordinate.::

   rho = lfield.calc_rhoimpact(obj)

Observed Galaxies
-----------------

Generate a table of the observed galaxies in the field
within a given angular radius of the field coord.::

   targ, dates, idx = lfield.get_observed(5.*u.arcmin)

Unobserved Galaxies
-------------------

Generate a table of the unobserved galaxies in the field
within a given angular radius of the field coord.::

   need_targ = lfield.get_unobserved(5.*u.arcmin)

Associated Galaxies
-------------------

Generate a table of the galaxies "associated" to a given LOS.::

   close_gal, rho = lfield.get_associated_galaxies(0.13, R=300*u.kpc)

Mask Date
---------

Returns a list of the date(s) when a given mask was observed.::

   dates = lfield.get_mask_obsdate('PG1407_may_mid2')

Clean Duplicates
----------------

Returns a version of a table (e.g. targets, galaxies) without duplicates. The table
has to have columns for sky coordinates (e.g. RA, DEC) and the duplication criteria is
based on a angular tolerance (usually small; default is `tol = 1*u.arcsec`). Currently, the
duplication conflict is solved by only keeping the first entry but we expect other methods
will be available in the future.::

   targets = lfield.targets
   clean_targets = lfield.clean_duplicates(targets, tol=1*u.arcsec, method='first')

