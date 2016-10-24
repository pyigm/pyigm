
Ingesting Prochaska+11 CGM dataset [v1.0]
=========================================

.. code:: python

    # imports
    from astropy.table import Table
    from astropy.coordinates import SkyCoord
    
    import pyigm
    from pyigm.cgm.cgmsurvey import CGMAbsSurvey

Summary
-------

.. code:: python

    ovi_file = pyigm.__path__[0]+'/data/CGM/P11/lowovidat.fits'

.. code:: python

    ovidat = Table.read(ovi_file)

.. code:: python

    ovidat[0:3]




.. raw:: html

    &lt;Table length=3&gt;
    <table id="table4599547088">
    <thead><tr><th>QSO</th><th>QSO_RA</th><th>QSO_DEC</th><th>QSO_ZEM</th><th>QSO_VMAG</th><th>QSO_UV</th><th>FLG_GAL</th><th>GAL_FIL</th><th>R_LIMIT</th><th>N_GAL [4]</th><th>COMPLETE [2,4]</th><th>FLG_FUSE</th><th>FUSE_EXP</th><th>FUSE_SNR</th><th>FLG_STIS</th><th>STIS_COMM [4]</th><th>FLG_GHRS</th><th>GHRS_COMM [4]</th></tr></thead>
    <thead><tr><th>str15</th><th>str15</th><th>str15</th><th>float32</th><th>float32</th><th>float32</th><th>int16</th><th>str58</th><th>float32</th><th>int32</th><th>int32</th><th>int16</th><th>float32</th><th>float32</th><th>int16</th><th>str60</th><th>int16</th><th>str60</th></tr></thead>
    <tr><td>Q0026+1259</td><td>00:29:13.8</td><td>+13:16:04</td><td>0.142</td><td>14.78</td><td>20.0</td><td>1</td><td>/u/xavier/LCO/OVI/FUSE/data/Q0026+1259/Q0026+1259_gal.fits</td><td>20.0</td><td>131 .. 47</td><td>22 .. 85</td><td>1</td><td>20000.0</td><td>7.0</td><td>0</td><td>..</td><td>1</td><td>G270M 5222s                                                  ..</td></tr>
    <tr><td>TONS180</td><td>00:57:20.0</td><td>-22:22:56</td><td>0.06198</td><td>16.6</td><td>30.0</td><td>1</td><td>/u/xavier/LCO/OVI/FUSE/data/TONS180/TONS180_gal.fits</td><td>19.7</td><td>110 .. 4</td><td>15 .. 92</td><td>1</td><td>132453.0</td><td>15.0</td><td>2</td><td>G140M 7000s 15.                                              ..</td><td>0</td><td>..</td></tr>
    <tr><td>TONS210</td><td>01:21:51.5</td><td>-28:20:57</td><td>0.116</td><td>14.7</td><td>70.0</td><td>1</td><td>/u/xavier/LCO/OVI/FUSE/data/TONS210/TONS210_gal.fits</td><td>20.0</td><td>71 .. 5</td><td>6 .. 87</td><td>1</td><td>56500.0</td><td>20.0</td><td>2</td><td>E140M 22000s                                                 ..</td><td>0</td><td>..</td></tr>
    </table>



.. code:: python

    qso_radec = SkyCoord(ra=ovidat['QSO_RA'], dec=ovidat['QSO_DEC'], unit=(u.hourangle, u.deg))

Dwarfs
------

.. code:: python

    cgm_dwarf_file = pyigm.__path__[0]+'/data/CGM/P11/dwarf_galabs_strct.fits'

.. code:: python

    cgm_dwarfs = Table.read(cgm_dwarf_file)

.. code:: python

    cgm_dwarfs[0:3]




.. raw:: html

    &lt;Table length=3&gt;
    <table id="table4599717648">
    <thead><tr><th>FIELD</th><th>ID</th><th>OBJ_ID</th><th>FLG_ANLY</th><th>FLG_SURVEY</th><th>OBJ_TYPE</th><th>MAG [10]</th><th>MAGERR [10]</th><th>FILTER [10]</th><th>IMG_FIL [10]</th><th>XYPIX [2]</th><th>RA</th><th>DEC</th><th>AREA</th><th>STARGAL</th><th>GAL_TYPE</th><th>GAL_COEFF [10]</th><th>Z</th><th>VCIRC</th><th>FSPEC_FIL [10]</th><th>DRA</th><th>DDEC</th></tr></thead>
    <thead><tr><th>str11</th><th>int32</th><th>str1</th><th>int16</th><th>int16</th><th>int16</th><th>float32</th><th>float32</th><th>str1</th><th>str27</th><th>float32</th><th>float64</th><th>float64</th><th>float32</th><th>float32</th><th>str5</th><th>float32</th><th>float64</th><th>float32</th><th>str29</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>Q0026+1259</td><td>1303</td><td>a</td><td>7</td><td>1</td><td>0</td><td>19.733 .. 0.0</td><td>0.125 .. 9.99</td><td>B ..</td><td>Images/Q0026+1259XB.fits    ..</td><td>1014.76 .. 1034.29</td><td>7.28886499155</td><td>13.2746021629</td><td>10.43</td><td>0.18</td><td>Late</td><td>0.331349 .. 0.0</td><td>0.0329451337457</td><td>0.0</td><td>..</td><td>45.6990947951</td><td>0.0191390593786</td></tr>
    <tr><td>TONS180</td><td>2295</td><td>a</td><td>7</td><td>1</td><td>0</td><td>18.923 .. 1.0</td><td>0.088 .. 0.05</td><td>B ..</td><td>Images/TONS180XB.fits       ..</td><td>1318.89 .. 607.18</td><td>14.2667432785</td><td>-22.44755991</td><td>10.92</td><td>0.19</td><td>Late</td><td>-0.0115093 .. 0.0</td><td>0.0233643911779</td><td>0.0</td><td>..</td><td>154.07390626</td><td>0.0207054292147</td></tr>
    <tr><td>PKS0405-123</td><td>90033</td><td>a</td><td>7</td><td>1</td><td>0</td><td>0.0 .. 0.44</td><td>0.0 .. 0.06</td><td>..</td><td>..</td><td>61.9512 .. -12.1839</td><td>61.9512481689</td><td>-12.1838884354</td><td>0.0</td><td>0.0</td><td>Late</td><td>0.0 .. 0.0</td><td>0.167</td><td>1.0</td><td>..</td><td>97.7760719383</td><td>0.0877822129587</td></tr>
    </table>



Funny columns -- Renaming
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    cgm_dwarfs.rename_column('DRA', 'rho(kpc)')
    cgm_dwarfs.rename_column('DDEC', 'L(L*)')

MAG has HI
^^^^^^^^^^

::

    2 == Ref
    3 == z
    4 == Lya EW
    5 == sigEW Lya
    8 == NHI
    9 == sigNHI

.. code:: python

    cgm_dwarfs[1]['MAG']




.. parsed-literal::

    array([  1.89230003e+01,   1.80701466e+01,   5.00000000e+00,
             2.34109219e-02,   2.22000000e+02,   2.90000000e+01,
             5.60000000e+01,   4.00000000e+00,   1.38000002e+01,
             1.00000000e+00], dtype=float32)



MAGERR has OVI
^^^^^^^^^^^^^^

::

    2 == Ref
    3 == z
    4 == EW 1031
    5 == sigEW 1031
    8 == NOVI    
    9 == sigNOVI

.. code:: python

    cgm_dwarfs[1]['MAGERR']




.. parsed-literal::

    array([  8.79999995e-02,   6.49999976e-02,   4.00000000e+00,
             2.33999994e-02,   4.30000000e+01,   1.50000000e+01,
             3.00000000e+01,   0.00000000e+00,   1.34799995e+01,
             5.00000007e-02], dtype=float32)



Refs
^^^^

.. code:: python

    refdict = {1: 'tripp08', 2: 'tc08a', 3: 'ds08', 4: 'dsr+06', 5: 'pss04', 6: 'cm09', 9: 'p+11'}

Ingest
------

::

    python ingest_lit.py

Read
----

.. code:: python

    p11_tarfile = pyigm.__path__[0]+'/data/CGM/P11/P11_sys.tar'

.. code:: python

    p11 = CGMAbsSurvey.from_tarball(p11_tarfile, chk_lowz=False)


.. parsed-literal::

    WARNING: UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: 10**-1 nm. [astropy.units.format.utils]
    /Users/xavier/local/Python/linetools/linetools/lists/linelist.py:374: RuntimeWarning: divide by zero encountered in log10
      self._data['log(w*f)'] = np.log10(qm_strength)
    /Users/xavier/anaconda/lib/python2.7/site-packages/numpy/ma/core.py:824: RuntimeWarning: invalid value encountered in less_equal
      return umath.less_equal(x, self.critical_value)


.. parsed-literal::

    Loading abundances from Asplund2009
    Abundances are relative by number on a logarithmic scale with H=12
    Skipping a likely folder: CGM_JSON


.. parsed-literal::

    /Users/xavier/local/Python/linetools/linetools/isgm/abssystem.py:288: UserWarning: Input AbsComponent with Zion=(8, 6) does not match AbsSystem rules. Not appending
      warnings.warn('Input AbsComponent with Zion={} does not match AbsSystem rules. Not appending'.format(abscomp.Zion))
    /Users/xavier/local/Python/linetools/linetools/isgm/abssystem.py:294: UserWarning: Failed velocity overlap
      warnings.warn('Failed velocity overlap')
    /Users/xavier/local/Python/linetools/linetools/isgm/abssystem.py:288: UserWarning: Input AbsComponent with Zion=(1, 1) does not match AbsSystem rules. Not appending
      warnings.warn('Input AbsComponent with Zion={} does not match AbsSystem rules. Not appending'.format(abscomp.Zion))


.. code:: python

    p11




.. parsed-literal::

    <CGM_Survey:  nsys=54, ref=>
    <IGMSystem: IGMSystem 00:29:13.8 13:16:04, 0.0329451, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 00:29:13.8 13:16:04, 0.039311, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 00:57:20 -22:22:56, 0.0233644, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 00:57:20 -22:22:56, 0.045619, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 03:11:55.2 -76:51:51, 0.202643, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 03:11:55.2 -76:51:51, 0.0593531, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.0964516, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.297609, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.352, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.153212, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.203022, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.36124, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.167, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 04:07:48.4 -12:11:37, 0.167043, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 10:07:26.1 12:48:56, 0.0296659, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 10:07:26.1 12:48:56, 0.0092207, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 10:31:54.3 -14:16:51, 0.0508329, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 11:19:08.7 21:19:18, 0.165951, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 11:19:08.7 21:19:18, 0.0600208, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 11:19:08.7 21:19:18, 0.0593876, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 11:19:08.7 21:19:18, 0.13829, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:14:17.7 14:03:13, 0.06438, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:14:17.7 14:03:13, 0.0646084, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:14:17.7 14:03:13, 0.0519877, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:14:17.7 14:03:13, 0.0511341, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:19:20.9 06:38:38, 0.124102, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:19:20.9 06:38:38, 0.0131789, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:19:20.9 06:38:38, 0.00666757, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:19:20.9 06:38:38, 0.0080957, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:29:06.7 02:03:09, 0.00620912, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:33:25.8 09:31:23, 0.0118122, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:33:25.8 09:31:23, 0.125431, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:33:25.8 09:31:23, 0.0905207, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 12:33:25.8 09:31:23, 0.206801, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 13:05:33 -10:33:19, 0.191709, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 13:05:33 -10:33:19, 0.0935802, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 13:05:33 -10:33:19, 0.0425621, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 13:05:33 -10:33:19, 0.145303, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 13:09:47 08:19:49, 0.127621, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 13:09:47 08:19:49, 0.0337313, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 14:29:06.4 01:17:06, 0.0299413, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 14:29:06.4 01:17:06, 0.0281113, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 15:55:43 11:11:24, 0.0150682, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 15:55:43 11:11:24, 0.0395311, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 15:55:43 11:11:24, 0.0420751, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 15:55:43 11:11:24, 0.0416811, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:55:01.5 -09:22:25, 0.13262, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:55:01.5 -09:22:25, 0.0503651, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:55:01.5 -09:22:25, 0.155532, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:55:01.5 -09:22:25, 0.0807993, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:55:01.5 -09:22:25, 0.0788036, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:55:01.5 -09:22:25, 0.0809754, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:55:01.5 -09:22:25, 0.0516533, NHI=0, Z/H=0>
    <IGMSystem: IGMSystem 21:58:51.8 -30:13:30, 0.0169383, NHI=0, Z/H=0>



.. code:: python

    p11.rho




.. math::

    [46.32007,~231.08717,~153.91919,~281.05024,~34.785816,~241.08236,~271.56657,~263.85542,~170.30404,~193.11905,~277.85574,~228.92416,~101.21597,~122.3351,~181.21741,~78.670472,~306.79467,~159.32049,~133.01256,~221.58367,~137.97841,~151.94839,~72.602911,~182.7867,~138.10197,~94.401827,~106.07147,~92.907301,~37.14687,~85.547486,~41.258645,~260.27253,~120.93128,~243.45857,~225.98147,~69.865838,~227.78876,~88.463791,~99.493909,~274.85217,~162.60742,~309.02131,~176.70156,~282.59707,~293.98086,~184.38542,~228.32783,~305.40178,~266.93147,~292.71387,~236.98361,~35.061178,~271.77462,~108.99879] \; \mathrm{kpc}



