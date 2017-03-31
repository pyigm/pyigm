*************
pyigm Scripts
*************

There are a number of scripts, many of which are GUIs,
provided with pyigm.  As regards the GUIs we warn
again that Mac users will need to set their matplotlib to
something other than MacOSX. See
`backends <http://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`__.

pyigm_igmguesses
----------------

GUI for examining IGM sightlines.
See `ref:igmguesses` for details.

pyigm_mkigmsys
--------------

Generate a simple JSON file for an IGM System.
Useful for things like the GUI lt_xabssys in `linetools`.

Here is the usage::

    usage: pyigm_mkigmsys [-h] [--jcoord JCOORD] [--zem ZEM] [--sigNHI SIGNHI]
                          [--vlim VLIM]
                          itype zabs NHI outfile

    Show contents of a JSON file, assuming one of several formats (v1.0).

    positional arguments:
      itype            Type of IGMSystem: dla, lls
      zabs             Absorption redshift
      NHI              log10 NHI value
      outfile          Name of JSON file to create

    optional arguments:
      -h, --help       show this help message and exit
      --jcoord JCOORD  Coordinates in JXXXXXXXX.X+XXXXXX.X format
      --zem ZEM        Emission redshift
      --sigNHI SIGNHI  Error in NHI
      --vlim VLIM      Velocity limits in format ###,###

