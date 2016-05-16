.. highlight:: rest

.. _IGMSystem:

*********
IGMSystem
*********

.. index:: IGMSystem

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <IGMSystem_examples>
   LLS <LLSSystem_examples>
   DLA <DLASystem_examples>

Overview
========

The IGMSystem class is designed to organize and analyze an
absorption system.  This is generally constructed
of one or more (:ref:`AbsComponent`).
The base class is abstract, i.e. one must instantiate
one of its flavors (e.g. HILyman, MgII, LLS, DLA).

By definition, an IGMSystem is a unique collection of
absorption components.  It is specified by:

=============== ========   ============== ============================================
Property        Variable   Type           Description
=============== ========   ============== ============================================
Type            abs_type   str            Name of the type of system
RA, DEC         radec      tuple or coord RA,DEC in deg or astropy.coordinate
Redshift        z          float          absorption redshift
Velocity limits vlim       Quantity array -/+ velocity limits of the system
=============== ========   ============== ============================================


Instantiation
=============

The IGMSystem Class may be instantiated in a few ways.
The default sets the properties listed above::

	gensys = GenericIGMSystem((15.23*u.deg,-23.22*u.deg), 1.244, [-500,500]*u.km/u.s, NHI=16.)

More commonly, one will instantiate with one or more AbsComponent objects::

    # HI Lya, Lyb
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = radec
    # Finish
    HIsys = HISystem.from_components([abscomp])



Attributes
==========

Sub Classes
===========

Generic
-------

A catch-all subclass for IGMSystem.

HI
--

A subclass for HI systems.  All of the components are
required to have Zion=(1,1).

LLS
---

An IGMSystem class for LLS.  One is required to set the
NHI keyword upon instantiation.  See the :doc:`LLSSystem_examples`
notebook for more.

DLA
---

An IGMSystem class for DLA.  One is required to set the
NHI keyword upon instantiation and it must equal or
exceed 20.3.  See the :doc:`DLASystem_examples`
notebook for more.

Plots
=====

Methods
=======

Output
======
