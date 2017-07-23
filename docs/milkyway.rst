.. highlight:: rest

**************
Galaxy Classes
**************

.. index:: Galaxy

Notebooks
=========

.. toctree::
   :maxdepth: 1

Overview
========

There classes designed for analysis of gas in our Galaxy.
We describe each in turn.

GalaxyCGM
=========

This class enables analysis of the CGM of the Milky Way.
It is a child of :ref:`cgm-class`.

Thus far, the class only loads observational data.

=========   =============== ============================================
Data        Type            References
=========   =============== ============================================
OVII        Absorption line Fang+15
=========   =============== ============================================

Instantiation
-------------

Simply::

   mwcgm = GalaxyCGM()


