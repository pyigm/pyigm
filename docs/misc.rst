.. highlight:: rest

*************
Miscellaneous
*************

.. index:: misc

Methods
=======

DM
--

This method calculates the Dispersion Measure through the IGM, e.g.
`Ioka+03 <http://adsabs.harvard.edu/abs/2003ApJ...598L..79I>`_.
Specifically, we calculate (numerically):

.. math::
   {\rm DM}_{\rm IGM} = \frac{3 c H_o \Omega_b}{8 \pi G m_p}
     \int_0^z \, \frac{(1+z) dz}{[ \Omega_m(1+z)^3 + \Omega_\Lambda]^{1/2}}

Trivially performed with::

   from pyigm.fN import tau_eff as pyteff
   DM = pyteff.DM(z)



