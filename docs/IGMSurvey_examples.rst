
Examples for IGMSurvey (v1.2)
=============================

.. code:: python

    # imports
    from astropy.coordinates import SkyCoord
    
    from pyigm.abssys.igmsys import IGMSystem
    from pyigm.abssys.igmsurvey import GenericIGMSurvey

Simple instantiation
--------------------

.. code:: python

    gensurvey = GenericIGMSurvey()
    gensurvey




.. parsed-literal::

    [IGMSurvey: nsys=0, type=Generic, ref=]



.. code:: python

    coord = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys1 = IGMSystem('MgII', coord, 1.244, [-300,300.]*u.km/u.s, NHI=16.)
    gensys1.name = 'Sys1'
    #
    coord2 = SkyCoord(ra=223.1143*u.deg, dec=42.4321*u.deg)
    gensys2 = IGMSystem('MgII', coord2, 1.744, [-300,300.]*u.km/u.s, NHI=17.)
    gensys2.name = 'Sys2'

.. code:: python

    gensurvey.add_abs_sys(gensys1)
    gensurvey.add_abs_sys(gensys2)

.. code:: python

    gensurvey.abs_sys()




.. parsed-literal::

    array([<IGMSystem: MgII 08:12:27.432 -12:25:55.56, 1.244, NHI=16, Z/H=0>,
           <IGMSystem: MgII 14:52:27.432 42:25:55.56, 1.744, NHI=17, Z/H=0>], dtype=object)



Parsing
-------

.. code:: python

    gensurvey.NHI




.. parsed-literal::

    array([ 16.,  17.])



.. code:: python

    gensurvey.zabs




.. parsed-literal::

    array([ 1.244,  1.744])


