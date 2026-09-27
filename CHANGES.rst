2026 (unreleased)
-----------------

Packaging
~~~~~~~~~

- **Python 3.11 is now the minimum supported version.**  Python 2.7 through
  3.10 are no longer supported.
- Packaging migrated from a hand-edited ``setup.py`` to PEP 621
  ``pyproject.toml``.  ``setup.py``, ``setup.cfg`` and ``.travis.yml`` are gone.
- **Runtime dependencies are now declared.**  Previously ``pip install pyigm``
  installed no dependencies at all.  pyigm now requires ``numpy>=2.0``,
  ``scipy>=1.10``, ``astropy>=5.3``, ``matplotlib>=3.7``, ``h5py>=3.8``,
  ``PyYAML>=6.0`` and ``linetools>=0.3.2``.  Optional extras: ``test``,
  ``mcmc``, ``gui``, ``docs``.
- The version number is now derived from git tags via ``setuptools_scm``
  instead of being hard-coded to ``0.1.dev0``.
- ``pyigm.__version__`` is now available.
- The ``bin/`` wrapper scripts were replaced by console-script entry points.
  The six ``pyigm_*`` commands are unchanged from a user's point of view.
  (``bin/pyigm_mtlmcmc`` had been broken regardless: it called a ``parser()``
  function that no longer exists in that module.)

Bug fixes
~~~~~~~~~

- ``imp`` (removed from the standard library in Python 3.12) and
  ``pkg_resources`` (no longer shipped by setuptools) are no longer used.
  Both are replaced by :mod:`pyigm.resources`, built on
  ``importlib.resources``.  Without this, 19 of 35 test modules could not even
  be imported.
- ``astropy.cosmology.core`` (removed in astropy 8) replaced with the public
  ``astropy.cosmology.FlatLambdaCDM`` at seven sites.
- ``scipy.integrate.simps`` (removed in scipy 1.14) replaced with ``simpson``
  in :mod:`pyigm.euvb.cuba`.
- ``np.trapz`` (removed in numpy 2.0) replaced with ``np.trapezoid`` in
  :mod:`pyigm.clustering.randist`.
- ``h5py`` ``Dataset.value`` (removed in h5py 3.0) replaced with ``ds[()]`` in
  the COS-Halos and metallicity readers.
- ``tau_eff.lyman_ew`` no longer references the renamed data file
  ``EW_SPLINE_b24.yml``.  That file became ``EW_SPLINE_b24_2017.yml`` in
  November 2022 and the function has raised ``FileNotFoundError`` ever since.
- ``tau_eff.lyman_ew`` accepts ``wrest`` as either a Quantity or a bare
  Angstrom array.  Previously a Quantity was swallowed by a bare ``except:``
  that dropped the caller into ``pdb`` -- which, under multiprocessing, hung
  the worker forever.
- ``FNModel`` is picklable again, so ``tau_eff`` can be evaluated in parallel.
  scipy 1.18's ``PchipInterpolator`` holds an unpicklable module reference.
- ``IGMSurvey.calculate_gz`` no longer seeds ``np.arange`` with a float32
  table value, which cost ~4e-4 in redshift by the 4000th grid point.
- ``GalaxyCGM.load_coolgas`` now recognizes masked ``logN`` entries in the
  Richter+17 table.  ``np.isnan()`` of a masked element is itself masked, so
  those rows were silently being given ``flag_N=1`` and a column density of
  zero.
- ``IGMSightline.write_to_igmguesses`` no longer imports the Qt GUI module.
  Writing a JSON file required PyQt5 to be installed.
- :mod:`pyigm.metallicity.generate_pdf` no longer imports the defunct
  ``xastropy``, and no longer calls the removed ``np.random.random_integers``.
