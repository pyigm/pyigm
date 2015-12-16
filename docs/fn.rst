.. highlight:: rest

**************
fN Model Class
**************

.. index:: fNmodel

Notebooks
=========

.. toctree::
   :maxdepth: 1

   fN Model <fN_examples>

Overview
========

This Class is designed to enable f(N) calculations in the IGM.
While the authors have reservations on the validity of f(N)
as a description of the IGM, it still offers good value for
a number of calculations.


The FNModel Class may take several forms as described below.

========= ======================================================= ============
Type      Description                                             Reference
========= ======================================================= ============
Hspline   Hermite spline which is constrained to monotonically    `Prochaska+14 <http://adsabs.harvard.edu/abs/2014MNRAS.438..476P>`_
          decrease
Gamma     Power-law + exponential                                 `Inoue+14 <http://adsabs.harvard.edu/abs/2014MNRAS.442.1805I>`_
========= ======================================================= ============

Instantiation
=============


Here are some examples of instantiation::

    from pyigm.fN.fnmodel import FNModel
    # P14 HSpline
    fN_P14 = FNModel('Hspline', zmnx=(2.,5.))
    # I14
    fN_I14 = FNModel('Gamma')


And the default model is currently the P14 formulation::

   fN_default = FNModel.default_model()


Attributes/Properties
=====================

========   ============== ============================================
Variable   Type           Description
========   ============== ============================================
zmnx       tuple          min/max redshift for using this model to evaluate f(N)
========   ============== ============================================

Methods
=======

l(X)
----

Calculate :math:`\ell(X)`, the incidence of absorption per absorption
path length :math:`dX` in a given interval of column density:

.. math::
   \ell(X) = \int_{N_{\rm min}}^{N_{\rm max}} f(N,X) \, dN

Easily evaluated with::

   lX = fN_default.calculate_lox(z, Nmin, Nmax)

tau_eff^LL
----------

Calculate the effective optical depth to Lyman limit photons

.. math::
   \tau_{\rm eff}^{\rm LL}(z_{912},z_\mathrm{q}) = \int_{z_{912}}^{z_q} \int_0^\infty f(N_{\rm HI},z) \, [1 - \exp[-N_{\rm HI} \sigma_{\rm ph}(z')]] \, dN_{\rm HI} dz'

Easily evaluated over across a redshift interval::

   from pyigm.fN import tau_eff as pyteff
   zval,teff_LL = pyteff.lyman_limit(fN_default, 0.5, 2.45)

MFP
---

Evaluate the mean free path to ionizing photons at a given redshift::

   mfp = fN_default.mfp(z)

rho HI
------

Evaluate the mass density in HI atoms at a specific redshift, integrating
over a :math:`N_{\rm HI}` interval:

.. math::
   \rho_{\rm HI} = \frac{m_p H_0}{c} \int_{N_{\rm min}}^{N_{\rm max}} f(N_{\rm HI},X) \, N_{\rm HI} \, dN_{\rm HI}

As such::

   rho_HI = fN_default.calculate_rhoHI(z, (Nmin, Nmax))

Most useful for DLA calculations.