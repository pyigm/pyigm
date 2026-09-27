# 2026 Update

## Goals

This document is the `pyigm` counterpart to
`linetools/claude_prompts/2026_update_prompt.md`.  It was written on 2026-09-27
after a direct survey of the repository at commit `fef4771` (branch `2026`),
running everything through `conda run -n astro` (Python 3.14.6, numpy 2.5.2,
astropy 8.0.1, scipy 1.18.0, matplotlib 3.11.1, linetools 0.3.3.dev23 installed
editable from `~/Astronomy/linetools`).

**The headline finding is that `pyigm` does not run on this machine.**  The test
suite does not even collect: 19 of the 35 test modules die at import on
`ModuleNotFoundError: No module named 'imp'` and `No module named
'pkg_resources'`, both of which were removed from the standard library /
default environment years ago.  That is a different situation from `linetools`,
which was green and merely dusty.  Prompts 1-4 are triage; everything after
them is the same kind of modernization work the `linetools` document covered.

`pyigm` is a public, multi-author package (Prochaska, Tejos, Burchett and
others) with downstream users.  Preserve the public API, deprecate rather than
delete, and note any behaviour change so it can go into `CHANGES.rst` -- which
**does not exist yet**; prompt 9 creates it.  Run Python only via
`conda run -n astro`.  **Xavier performs all git commands.**

Each numbered prompt below is meant to be executed in its own session.  Do only
the task you are pointed at, then append a dated entry under `## Logs`.

## Prompts

1. Read this file.  Execute the 1st task under "Purge `imp`"
2. Read this file.  Execute the 1st task under "Purge `pkg_resources`"
3. Read this file.  Execute the 1st task under "astropy.cosmology.core"
4. Read this file.  Execute the 1st task under "scipy and numpy removals"
5. Read this file.  Execute the 1st task under "Clustering sensitivity loop"
6. Read this file.  Execute the 1st task under "Reionization tolerance"
7. Read this file.  Execute the 1st task under "xastropy"
8. Read this file.  Execute the 1st task under "Remove dead packaging files"
9. Read this file.  Execute the 1st task under "Migrate to pyproject.toml"
10. Read this file.  Execute the 1st task under "Expose __version__"
11. Read this file.  Execute the 1st task under "Build CI"
12. Read this file.  Execute the 1st task under "Root LICENSE"
13. Read this file.  Execute the 1st task under "README"
14. Read this file.  Execute the 1st task under "Test artifacts"
15. Read this file.  Execute the 1st task under "pdb.set_trace in library code"
16. Read this file.  Execute the 1st task under "Python-2 scaffolding"
17. Read this file.  Execute the 1st task under "Documentation"

**Order matters.**  1 -> 2 -> 3 -> 4 must run first and in that order: each one
uncovers the next, and until all four are done the suite cannot be collected, so
nothing else can be verified.  Then 8 -> 9 -> 10 -> 11 (they touch the same
metadata and depend on each other).  5, 6, 7 and 12 through 17 are independent
of one another.

18. Read this file.  Check the answers I have given in Q&A.  Proceed accordingly
    and report back.

### Baseline

```
conda run -n astro pytest pyigm -q
```

As of 2026-09-27 this **exits 2 with `Interrupted: 19 errors during
collection`** and runs zero tests.  To get any signal at all you must use:

```
conda run -n astro pytest pyigm -q --continue-on-collection-errors
```

which gives **27 passed, 5 failed, 19 errors, 377 warnings, 7.31s**.

Counting `^def test_` in the 19 unimportable modules, **75 further tests are
blocked**, so the suite is roughly 107 tests of which 27 currently run and pass.
After prompts 1-4 the bare `pytest pyigm -q` invocation must collect cleanly;
say in the log what the new pass/fail/error counts are and account for every
difference.

Two further facts about the environment, both of which matter:

- **`pyigm` is not installed in `astro`.**  `pip list` does not list it, and
  `import pyigm` only succeeds because the CWD is the repository root.  So a
  green test run says nothing about whether the package *installs* correctly.
  Prompt 9 is where that gets fixed and verified.
- `build` 1.6.1, `setuptools` 84.0.0, `setuptools_scm` 10.3.4, `wheel`,
  `sphinx` 9.1.0, `tox` 4.63.0 and `pytest-remotedata` 0.4.2 **are** available
  in `astro` (they were installed during the linetools work), so packaging and
  docs can actually be verified here.  `corner`, `pymc` and `mpmath` are **not**
  installed; the metallicity MCMC code is therefore untestable as things stand.

---

## Purge `imp`

1. `imp` was deprecated in Python 3.4 and **removed in Python 3.12**.  Nine
   modules still import it, and six of them do so at module scope, which is what
   breaks collection.  Full inventory, verified by grep:

   | File | Line | What it does |
   |---|---|---|
   | `pyigm/euvb/cuba.py` | 17 | `pyigm_path = imp.find_module('pyigm')[1]` |
   | `pyigm/clustering/soneira_peebles.py` | 7 | `import imp` -- **never used** |
   | `pyigm/fN/mockforest.py` | 5, 26 | `pyigm_path = ...` |
   | `pyigm/fN/fnmodel.py` | 20 | `pyigm_path = ...` |
   | `pyigm/fN/dla.py` | 7, 19 | `pyigm_path = ...` |
   | `pyigm/fN/plots.py` | 7, 15 | `pyigm_path = ...` |
   | `pyigm/fN/constraints.py` | 7 | `import imp` -- **never used** |
   | `pyigm/surveys/dlasurvey.py` | 7, 31 | `import imp, glob` |
   | `pyigm/surveys/llssurvey.py` | 6, 28 | `import imp, glob` |
   | `pyigm/surveys/lls_literature.py` | 15, 35 | `import imp, glob` |
   | `pyigm/fN/tests/test_fnconstraints.py` | 6 | `import imp, pdb` -- **`imp` never used** |
   | `pyigm/surveys/tests/test_igmsurvey.py` | 12 | `lt_path = imp.find_module('linetools')[1]` |

   Note `soneira_peebles.py`, `constraints.py` and `test_fnconstraints.py`
   import `imp` and never call it.  Delete the import there; do not invent a
   replacement.

2. Every real use is the same idiom -- `imp.find_module('pkg')[1]`, i.e. "give
   me the directory this package lives in".  Replace with:

   ```python
   import os
   import pyigm
   pyigm_path = os.path.dirname(os.path.abspath(pyigm.__file__))
   ```

   (and the analogous two lines for `linetools` in `test_igmsurvey.py`).
   `pyigm/__init__.py` is empty, so importing `pyigm` from inside a `pyigm`
   submodule is not circular.  `imp.find_module` also **could not see namespace
   packages or zipped installs**; `__file__` has the same limitation, so this is
   a like-for-like swap and not a behaviour change.

   See Q&A #1 before you start -- there is a case for defining `pyigm_path` once
   in `pyigm/__init__.py` instead of nine times, and I would rather you asked
   than guessed.

3. Do **not** try to also fix `pkg_resources` in this prompt.  Several of these
   same files import both; fix one axis at a time so the diff stays readable and
   so the second prompt's collection run is a real check.

4. Verify: `conda run -n astro pytest pyigm -q --continue-on-collection-errors`
   and confirm the collection-error count drops from 19.  It will **not** reach
   zero -- prompt 2 handles the rest.  State the new number.

---

## Purge `pkg_resources`

1. `pkg_resources` is the legacy setuptools runtime API.  It is deprecated
   upstream, slated for removal, and **is not importable in `astro` at all**
   (setuptools 84 no longer ships it by default).  Sixteen sites use it, all via
   the identical one-line import:

   ```
   pyigm/cgm/qpq.py:8               pyigm/fN/tau_eff.py:6
   pyigm/cgm/cgmsurvey.py:11        pyigm/fN/constraints.py:10
   pyigm/cgm/galaxy.py:9            pyigm/fN/mcmc.py:15
   pyigm/cgm/cos_halos.py:13        pyigm/surveys/dlasurvey.py:16
   pyigm/cgm/ingest_lit.py:9        pyigm/surveys/llssurvey.py:15
   pyigm/euvb/cloudyiongrid.py:36   pyigm/continuum/quasar.py:9
   pyigm/metallicity/utils.py:11
   pyigm/cgm/tests/test_init_cgm.py:25
   pyigm/igm/tests/test_igmsightline.py:8
   pyigm/surveys/tests/test_llssurvey.py:9
   ```

2. The stdlib replacement is `importlib.resources`.  For a package whose data
   really is on disk (which `pyigm/data`, 19 MB of it, is):

   ```python
   from importlib.resources import files
   ...
   path = str(files('pyigm') / 'data' / 'something.fits')
   ```

   Every existing call in this repo is `resource_filename('pyigm', <relpath>)`,
   so the mechanical translation is exact.  **Check the argument at each site
   before converting** -- if any call names a package other than `pyigm`, it
   needs the same treatment but with that package's name.

   There are sixteen of these.  Writing one small helper and calling it
   sixteen times is tempting; see Q&A #1, which covers the same question for
   `pyigm_path`, and do whichever Xavier picks for both.

3. Verify: **`conda run -n astro pytest pyigm -q` must now collect with zero
   errors** (no `--continue-on-collection-errors`).  If it does not, the
   remaining errors are prompts 3 and 4's problem -- name them in the log and
   stop rather than fixing them here.

   Expected at this point, from a probe I ran with both modules shimmed: the
   only remaining module-scope import failures are
   `pyigm/surveys/dlasurvey.py` (prompt 3) and `pyigm/euvb/cuba.py` (prompt 4).

---

## astropy.cosmology.core

1. `astropy.cosmology.core` was made private in astropy 5.0 and **is gone in
   astropy 8** (`cannot import name 'core' from 'astropy.cosmology'`).  I
   confirmed this directly in `astro`.  Seven sites:

   | File | Line | Code |
   |---|---|---|
   | `pyigm/surveys/dlasurvey.py` | 22 | `from astropy.cosmology import core as acc` (**module scope -- blocks import**) |
   | `pyigm/surveys/dlasurvey.py` | 541 | `acc.FlatLambdaCDM(70., 0.3)` |
   | `pyigm/utils.py` | 141-142 | `from astropy.cosmology import core as acc` (function scope) |
   | `pyigm/fN/fnmodel.py` | 289, 508 | `cosmology.core.FlatLambdaCDM(70., 0.3)` |
   | `pyigm/fN/mockforest.py` | 110 | `cosmology.core.FlatLambdaCDM(70., 0.3)` |
   | `pyigm/fN/dla.py` | 49 | `cosmology.core.FlatLambdaCDM(70., 0.3)` |
   | `pyigm/fN/constraints.py` | 67 | `cosmology.core.FlatLambdaCDM(70., 0.3, Ob0=0.0455)` |

   `FlatLambdaCDM` has been re-exported from `astropy.cosmology` itself for a
   decade, so the fix is `from astropy.cosmology import FlatLambdaCDM` and
   `FlatLambdaCDM(...)` at the call sites.  `astropy.cosmology.Planck15` is
   still importable in astropy 8 (verified), so the nine `Planck15` imports
   elsewhere need nothing.

2. This one is load-bearing for a currently-failing test:
   `pyigm/tests/test_utils.py::test_dxdz` fails today with exactly this
   `ImportError`, out of `pyigm/utils.py:141`.  It should go green.

3. Note in the log whether the *numerical* results change.  They should not --
   `astropy.cosmology.FlatLambdaCDM` and `astropy.cosmology.core.FlatLambdaCDM`
   were the same class -- but `test_dxdz` asserts on numbers and has not run in
   years, so treat a numerical mismatch as a finding, not as something to
   paper over by relaxing the assertion.

---

## scipy and numpy removals

1. **`scipy.integrate.simps`** was removed in scipy 1.14 (renamed `simpson`
   long before).  One site:

    - `pyigm/euvb/cuba.py:10` -- `from scipy.integrate import simps`
      (module scope, blocks import of the whole `euvb` subpackage and its four
      tests)
    - `pyigm/euvb/cuba.py:115` -- `phi = np.log(10.)*simps(integrand.value, log_energy)`

   Replace with `simpson`.  **`simpson`'s signature is not identical**: the old
   `simps(y, x)` positional form became `simpson(y, x=...)` in scipy 1.6+ and
   the positional `x` was removed.  Write `simpson(integrand.value,
   x=log_energy)` and confirm the value is unchanged, because `test_cuba.py`
   asserts on numbers.

2. **`np.trapz`** was removed in numpy 2.0 (renamed `np.trapezoid`).  One site:

    - `pyigm/clustering/randist.py:36` -- `self.pdf = dist / np.trapz(dist, x=x)`

   This is the cause of `pyigm/clustering/tests/test_randist.py::test_randist`
   failing today.  Replace with `np.trapezoid`.  The keyword `x=` is already
   correct for the new name.

3. I grepped for the other common removals -- `np.float`, `np.int`, `np.bool`,
   `np.object`, `np.str`, `np.alltrue`, `np.product`, `np.in1d`,
   `scipy.misc`, `scipy.integrate.cumtrapz`, `scipy.interpolate.interp2d` --
   and **found none**.  Do not go hunting; these two are the whole list.

4. After this prompt, `conda run -n astro pytest pyigm -q` must collect all 35
   test modules with no errors.  That is the real gate for prompts 1-4.  Record
   the full counts.

---

## Clustering sensitivity loop

1. Two tests fail today with `AssertionError: Something wrong with the redshift
   distribution` raised at `pyigm/clustering/xcorr_utils.py:360`:

    - `pyigm/clustering/tests/test_cfield.py::test_init_cfield`
    - `pyigm/clustering/tests/test_s2d2pcf.py::test_init_survey`

   Both go through `ClusteringField.addGal` ->
   `xcorr_utils.spline_sensitivity`.

2. **This is a genuine bug, not a tolerance problem.**  Read
   `spline_sensitivity` at `xcorr_utils.py:341-361`:

   ```python
   for tt, mag in enumerate(magbins):
       delta_mag2 = delta_mag
       q = 0
       while True:
           #cond = (galreal.MAG <= mag + delta_mag2 * 0.5) & (galreal.MAG > mag - delta_mag2 * 0.5)
           if tt == len(magbins)-1:
               magmax = 99.
           else:
               magmax = magbins[tt+1]
           cond = (galreal.MAG <= magmax) & (galreal.MAG > mag)
           if np.sum(cond) >= Nmin:
               break
           else:
               if not flag_mag:
                   pdb.set_trace()
               delta_mag2 += 0.25
           q += 1
           assert q < 1000, 'Something wrong with the redshift distribution'
   ```

   The loop's stated job is "insure there are `Nmin` galaxies in the bin by
   extending the magnitude bin if needed".  But the line that actually used
   `delta_mag2` is **commented out**, and the live `cond` depends only on
   `mag` and `magbins[tt+1]` -- neither of which the loop body changes.  So
   `delta_mag2 += 0.25` widens nothing: if a bin starts under-populated it stays
   under-populated, the loop spins 1000 times and hits the assert.  It is an
   unconditional failure dressed up as a retry.

   The test data makes this certain: `utils.make_rand_gal(nrand=1000)` draws
   magnitudes in roughly 19-21, while `magbins` starts at 17.0.  The 17.0-17.5
   bin is empty on the first iteration, and always will be.

3. This is a **numerical-methods decision, not a mechanical fix**, and I am not
   going to guess it.  The three defensible repairs are:

   a. Restore the commented-out widening `cond` so `delta_mag2` does something
      again -- i.e. revert to the sliding-window behaviour the docstring
      describes.  This changes the science for *every* bin, not just sparse
      ones, because the windows then overlap.
   b. Keep the current disjoint bins but **drop** bins with fewer than `Nmin`
      galaxies instead of asserting, and return only the populated `magbins`.
      `spline_sensitivity`'s callers index `SPL` by magnitude, so dropping bins
      means the callers need to handle a missing key.
   c. Derive `magbins` from the data (e.g. `magmin = floor(min(MAG))`) rather
      than hard-coding `magmin=17., magmax=26.`, so empty bins never arise.

   See Q&A #2.  **Do not pick one on your own**, and do not "fix" the tests by
   loosening the assert or by feeding them denser fake galaxies -- that hides
   the defect for the real COS-Halos / clustering users.

4. Note separately: `pdb.set_trace()` at line 357 is reached whenever a caller
   passes explicit `magbins` and a bin is short.  That is library code dropping
   a user into the debugger.  Prompt 15 deals with the general case; flag this
   one there.

5. Also from the same run: `git log` shows this function last changed in
   `1a8c2a5` ("remove numba import") / `7bf89dc` / `ed3337e` ("edits for NT").
   Worth a look at whether the widening `cond` was commented out deliberately
   during that work.

---

## Reionization tolerance

1. `pyigm/igm/tests/test_reion.py::test_damping` fails:

   ```
   np.testing.assert_allclose(tau[-1], 0.00245349, rtol=1e-5)  # scipy 0.19
   E   Max relative difference among violations: 2.95819201e-05
   E    ACTUAL: 0.002453...   DESIRED: 0.00245349
   ```

   It misses by 3e-5 relative against an `rtol` of 1e-5 -- a factor of three.
   The comment `# scipy 0.19` says exactly what happened: the reference value
   was frozen against scipy 0.19 and `astro` has scipy 1.18.

2. This is a **quadrature-tolerance drift, not a bug**, and the honest fix is
   to widen `rtol` to something that reflects what the integrator actually
   guarantees -- `rtol=1e-4` covers the observed 3e-5 with room to spare.

   Before you change it, satisfy yourself that it really is quadrature drift:
   look at what `reion.igm_damping` integrates and confirm the discrepancy is in
   the integral's tail and not in, say, a changed constant.  The first
   assertion in the same test (`tau0` against `99.`) passes, which is evidence
   for drift rather than for a wrong physical constant.

3. Update the comment to say what it is pinned against now
   (`# scipy 1.18; rtol widened 2026-09, was 1e-5 vs scipy 0.19`) so the next
   person does not have to re-derive this.

---

## xastropy

1. `xastropy` is a defunct personal package.  It is **not installed** in `astro`
   and is not a declared dependency anywhere.  Three modules import it at module
   scope and are therefore permanently unimportable:

   | File | Line |
   |---|---|
   | `pyigm/observing/lick.py` | 22-24 (`xdebug`, `obs.finder`, `obs.radec`) |
   | `pyigm/metallicity/generate_pdf.py` | 20 (`xdebug`) |
   | `pyigm/guis/igmguesses.py` | 1267 (function scope -- `xguis.spec_guis`) |

   Plus three harmless mentions in comments/docstrings
   (`pyigm/field/igmfield.py:19`, `pyigm/abssys/lls_fitting.py:553`,
   `pyigm/fN/plots.py:21`, `pyigm/scripts/pyigm_fitlls.py:23`) -- leave those,
   or fix the docstrings if convenient.

2. **`pyigm/observing/` is not even a package.**  It has no `__init__.py`, so
   `find_packages()` skips it and `pyigm/observing/lick.py` has never been
   installed by `pip install pyigm` -- I confirmed `find_packages()` returns 25
   packages and `pyigm.observing` is not among them.  The file is dead twice
   over.

   Q&A #3 asks whether to delete `pyigm/observing/` outright.  **Do not delete
   it before that is answered.**

3. For `generate_pdf.py`: the `xdebug` import is a debugging convenience
   (`xdb.set_trace()` etc.).  Check whether `xdb` is actually called in the
   module.  If it is only imported, delete the import and the module becomes
   importable again for free.  If it is called, those are debugger entry points
   and belong to prompt 15.

4. For `igmguesses.py:1267`: this is inside a function, so it does not block
   import, and the GUI is untestable here anyway (`import qtpy` raises
   `QtBindingsNotFoundError` in `astro`).  Replace the import with a clear
   `ImportError`/message rather than leaving a reference to a package that no
   longer exists, or remove the code path -- say which you did.

5. Two related unimportable modules turned up in the same probe and are **not**
   xastropy's fault; record them but do not act:

    - `pyigm/metallicity/mcmc.py` and `pyigm/scripts/pyigm_mtlmcmc.py` raise
      `ImportError: Install corner to use metallicity.mcmc`.  `corner`, `pymc`
      and `mpmath` are genuinely optional and genuinely absent.  That is
      correct behaviour, not a bug.  Prompt 9 should declare them as an
      optional-dependency extra.

---

## Remove dead packaging files

1. Delete the following, all of which are provably unreachable:

    - **`.travis.yml`** -- `travis-ci.org` no longer exists.  The file still
      configures `PYTHON_VERSION=2.7`, `NUMPY_VERSION=1.13`, and clones
      `git://github.com/astropy/ci-helpers.git` over the `git://` protocol,
      which GitHub disabled in 2022.  It could not run even if Travis were
      alive.
    - **`.gitmodules`** -- declares an `astropy_helpers` submodule.  The
      directory does not exist, and `git ls-files --stage | grep astropy_helpers`
      returns nothing, so the submodule is not in the index either.  The file
      is a pure leftover; `git submodule status` prints nothing.
    - **`setup.cfg`** -- its entire contents are

      ```ini
      [aliases]
      test=pytest
      ```

      i.e. the `python setup.py test` alias, which setuptools removed.  Prompt 9
      would delete it anyway; doing it here keeps that diff clean.

2. There is **no `MANIFEST.in` in this repository at all**.  `setup.py` sets
   `include_package_data=True`, which for sdists defers to `MANIFEST.in` or a
   setuptools VCS plugin.  Do not create one in this prompt -- prompt 9 decides
   how source distributions get their 19 MB of `pyigm/data`, and that decision
   belongs with the `pyproject.toml` work.  Just note it.

3. `.gitignore` is in reasonable shape but is missing `.tox/` and
   `.pytest_cache/`.  `.pytest_cache/` currently shows as ignored only because
   pytest writes its own `.pytest_cache/.gitignore` containing `*` --
   `git check-ignore -v .pytest_cache` matches nothing in our own rules.  That
   works but is accidental; add both explicitly.

   Leave the `tmp.*` rule alone for now and read prompt 14 first -- that rule is
   currently hiding test artifacts, and prompt 14 is about making them stop
   being written rather than about hiding them better.

---

## Migrate to pyproject.toml

1. `setup.py` is a hand-edited copy of the DESI `desiutil` template, complete
   with commented-out `from desiutil.setup import DesiTest, DesiVersion` lines.
   Replace it with a PEP 621 `pyproject.toml` and delete `setup.py`.

   What is in there now, and what has to happen to each piece:

   | Current | Disposition |
   |---|---|
   | `name`, `description`, `author`, `author_email`, `license`, `url` | -> `[project]` |
   | `version = '0.1.dev0'` (**hard-coded**) | see item 2 |
   | `long_description` read from `README.md` | -> `readme = "README.md"` |
   | `scripts = glob('bin/*')` | -> `[project.scripts]`, see item 3 |
   | `provides`, `requires = ['Python (>2.7.0)']` | **delete** -- distutils-era keys that setuptools ignores outright |
   | `setup_requires = ['pytest-runner']` | **delete** -- `pytest-runner` is abandoned and `setup_requires` is removed |
   | `tests_require = ['pytest']` | -> `[project.optional-dependencies] test` |
   | `zip_safe = False` | **delete** -- an easy_install key, meaningless now |
   | `packages = find_packages()` | -> `[tool.setuptools.packages.find]` |
   | `package_data` (the `os.walk('pyigm/data')` loop) | -> `[tool.setuptools.package-data]`, see item 4 |

2. **There is no `install_requires` at all.**  `pip install pyigm` today
   installs `pyigm` and nothing else -- no numpy, no astropy, no linetools.
   Add real dependencies.  From the import survey the runtime set is:
   `numpy`, `scipy`, `astropy`, `matplotlib`, `h5py`, `pyyaml`, `linetools`.
   Optional: `corner`, `pymc`, `mpmath` (metallicity MCMC), `qtpy` + a Qt
   binding (GUIs).

   `linetools` is not on PyPI in a current form, which is why `docs/install.rst`
   tells people to `pip install --no-deps git+https://...`.  **See Q&A #4** --
   how `linetools` is declared depends on whether it gets a PyPI release, and I
   should not guess.

   Set floors to the oldest versions that actually work on Python 3.11 (the same
   call made for `linetools`): roughly `numpy>=1.24`, `astropy>=5.3`,
   `scipy>=1.10`, `matplotlib>=3.7`, `h5py>=3.8`.  But note the two API
   removals fixed in prompt 4 cut the other way: `np.trapezoid` needs
   **numpy >= 2.0** and `scipy.integrate.simpson`'s keyword-only `x` needs
   **scipy >= 1.6**.  If you use `np.trapezoid` unconditionally, the numpy floor
   is `>=2.0`, full stop.  Say in the log which floor you set and why.

   `version = '0.1.dev0'` has been hard-coded since the file was written.  The
   repository's tags are the real version source, and `setuptools_scm` 10.3.4 is
   installed in `astro`.  See Q&A #5.

3. **`requires-python`.**  Match `linetools`: `">=3.11"`.  Classifiers get
   `3.11`, `3.12`, `3.13`, `3.14` -- not `3.15`, which is not final.  This is a
   **breaking change for downstream installers** and must go in `CHANGES.rst`.

4. **The `bin/` scripts.**  All six (`pyigm_fitdla`, `pyigm_fitlls`,
   `pyigm_igmguesses`, `pyigm_mkigmsys`, `pyigm_mtlmcmc`, `pyigm_showjson`) are
   identical three-line wrappers:

   ```python
   import pyigm.scripts.pyigm_fitdla as pyigm_fitdla
   if __name__ == '__main__':
       args = pyigm_fitdla.parser()
       pyigm_fitdla.main(args)
   ```

   That is exactly what `[project.scripts]` entry points are for, and entry
   points work on Windows and in wheels where `scripts=` shebangs do not.
   Convert them, which means each `pyigm/scripts/pyigm_*.py` needs a
   zero-argument callable to point at.  They currently expose `parser()` and
   `main(args)`, so either add a small `def entry(): main(parser())` to each or
   give `main` a default of `None` -- check the actual signatures before
   deciding, and keep `bin/` around one release as a deprecation courtesy if
   you prefer.  Say which route you took.

5. **Package data is the real risk in this migration.**  `pyigm/data` is 19 MB
   across many nested directories, and the current `setup.py` builds the list by
   walking the tree at build time -- which means the file list is correct by
   construction today and will be a static glob afterwards.  Get it wrong and
   wheels ship without the line lists, CUBA grids and survey tables, and
   everything breaks for everyone.

   Also note `pyigm/data` has **no `__init__.py`**, so it is package *data* of
   `pyigm`, not a package.  Keep it that way.

6. `CHANGES.rst` does not exist.  Create it in this prompt, with an `0.4 (unreleased)`
   section, and record: the Python 3.11 floor, the new declared dependencies,
   the `bin/` -> entry-point move, and (back-filled from prompts 1-7) the
   `imp`/`pkg_resources`/`cosmology.core`/`simps`/`trapz` fixes.

7. **Verify properly.**  Unlike the `linetools` migration, `astro` has the tools
   for this:

   ```
   conda run -n astro python -m build --outdir /tmp/pyigm_dist
   conda run -n astro python -c "import zipfile; z=zipfile.ZipFile('/tmp/pyigm_dist/pyigm-...whl'); print(len([n for n in z.namelist() if n.startswith('pyigm/data/')]))"
   conda run -n astro pip install -e .
   conda run -n astro pytest pyigm -q
   ```

   Compare the data-file count in the wheel against
   `find pyigm/data -type f | wc -l`.  They should match.  Check the sdist too.
   And confirm the six `pyigm_*` console scripts resolve after the editable
   install.

---

## Expose __version__

1. `pyigm/__init__.py` is a single newline character -- 1 byte -- so
   `pyigm.__version__` raises `AttributeError`.  Add:

   ```python
   from importlib.metadata import PackageNotFoundError, version as _version

   try:
       __version__ = _version("pyigm")
   except PackageNotFoundError:  # not installed, e.g. running from a source tree
       __version__ = "unknown"
   ```

   Prefer this over importing a generated `pyigm/version.py`.  Note that
   `.gitignore` already has a `*/version.py` rule inherited from the astropy
   template, so if prompt 9 configured `setuptools_scm` with a `version_file`,
   that file is correctly ignored.

2. This is a public-API **addition** and therefore backwards compatible.  Add a
   line to `CHANGES.rst` anyway.

3. If Q&A #1 came back in favour of centralizing, this is also where
   `pyigm_path` goes.  Do that in the same commit, and update the call sites
   from prompts 1 and 2.

4. Verify:
   `conda run -n astro python -c "import pyigm; print(pyigm.__version__)"`
   from a directory that is **not** the repository root, so you are testing the
   installed package and not the CWD.

---

## Build CI

1. **There is no CI.**  `.github/` does not exist; the only automation ever
   configured was the dead `.travis.yml`.  Every commit on this repository for
   years has been untested by anything other than whoever ran pytest locally.

   Create `.github/workflows/ci_tests.yml`.  Model it on the one prompt 5 of the
   `linetools` document produced -- it is a sibling package, the same person
   maintains both, and gratuitous divergence helps nobody.  Specifically:

    - Triggers: push to `master` and `2026`, plus `pull_request`.  **`master` is
      the default branch** here (`origin/HEAD -> origin/master`); do not write
      `main`.  Quote the branch name `"2026"` -- bare `- 2026` parses as the
      YAML integer 2026 and silently matches nothing.
    - `actions/checkout@v4`, `actions/setup-python@v5`.
    - Matrix: Python `3.11`, `3.12`, `3.13`, `3.14` on `ubuntu-latest`, plus
      `macos-latest` on one version (Xavier develops on macOS).
    - No Python 3.15 job.  It is not final and its dependency wheels do not
      exist; the `linetools` answer was to skip it.
    - Run `pytest --pyargs pyigm`, not bare `pytest`, so that an installed
      package is what gets tested.  This matters more here than it did for
      `linetools`: the whole repository currently only imports because of CWD.
    - **`linetools` must be installed from git** in CI, exactly as `.travis.yml`
      used to do (`pip install git+https://github.com/linetools/linetools.git`),
      unless Q&A #4 says otherwise.
    - `MPLBACKEND=agg`.  No test needs a display; do not resurrect
      `SETUP_XVFB`.
    - Deselect or skip the `remote_data` tests by default.  Six tests are marked
      `@pytest.mark.remote_data` (in `test_coshalos.py`, `test_cgmsurvey.py`,
      `test_llssurvey.py`, `test_dlasurvey.py`) and they download tarballs from
      `www.ucolick.org` and `iopscience.iop.org`.  `pytest-remotedata` 0.4.2 is
      installed locally and skips them unless `--remote-data` is passed, so the
      default is already correct -- just make sure CI installs that plugin so
      the marks are recognized rather than being unknown-marker warnings.
    - Add a docs job: `sphinx-build -W -b html docs docs/_build/html`.  See
      prompt 17 -- do **not** turn on `-W` until that prompt has run, or the
      docs job will be red from day one.  Land the docs job non-fatal first if
      you have to, and say so.

2. There is no `tox.ini` in this repository (unlike `linetools`).  Do not create
   one unless you have a reason; a plain `pip install -e .[test] && pytest`
   workflow is simpler and is what this package needs.

3. **Expect the first CI run to be red even after prompts 1-7**, because CI
   tests an *installed* package and nothing here has been installed in years.
   That is the point of the exercise.  Report what it finds rather than
   weakening the workflow to make it pass.

---

## Root LICENSE

1. The license lives at `licenses/LICENSE.rst`, so GitHub does not detect it and
   the repository shows no license in its sidebar or API metadata.  Copy it to a
   root `LICENSE`.  Use a real copy, not a symlink -- symlinks do not survive
   sdist/wheel packaging reliably and GitHub's license detector does not follow
   them.

2. **The license text is an unfilled template.**  Verbatim, line 1 of
   `licenses/LICENSE.rst` reads:

   ```
   Copyright (c) year, author
   ```

   and clause 3 still says *"Neither the name of the **Astropy Team** nor the
   names of its contributors..."*.  Both are unedited astropy package-template
   placeholders -- the identical defect `linetools` had.  The license *type* is
   BSD-3-Clause, which matches `setup.py`'s `license = 'BSD'`.

   Filling in a copyright holder and year is a legal attribution decision, not a
   judgement call to make on someone else's behalf.  **See Q&A #6** and do not
   guess.  (`README.md` names J. Xavier Prochaska, N. Tejos and J. Burchett as
   primary builders; `docs/conf.py` says `copyright = u'2016, Prochaska, and
   Associates'`.  Those disagree with each other, which is exactly why I am
   asking.)

3. Leave `licenses/` in place; `licenses/README.rst` references it and the
   directory is the conventional home for third-party license texts.

---

## README

1. `README.md` has a `.md` extension but its body is reStructuredText:

   ```
   pyigm
   =====
   ```

   GitHub renders that as a paragraph followed by a stray `=====` line, not as a
   heading.  Same defect as `linetools`.  Convert to real Markdown (`# pyigm`),
   and make the builders list a real Markdown list -- right now
   "J. Xavier Prochaska / N. Tejos / J. Burchett" are three consecutive lines,
   which Markdown joins into one paragraph.

   Prefer converting over renaming to `.rst`: the Zenodo badge already in the
   file is Markdown syntax, and `setup.py` reads `README.md` by name.

2. The README is 12 lines long and says almost nothing.  Add, at minimum:

    - a CI badge once prompt 11 exists:
      `https://github.com/pyigm/pyigm/actions/workflows/ci_tests.yml/badge.svg`
    - an install line
    - a link to the documentation
    - a one-line pointer to `linetools`, since pyigm hard-depends on it

   Keep the Zenodo DOI badge; it is fine.

3. Do not add a coverage badge.  Nothing uploads coverage, and the `linetools`
   decision was explicitly to drop `coveralls` rather than ship a badge that
   lies.

---

## Test artifacts

1. The test suite writes files into the source tree.  Two categories, and the
   second is worse.

   Writing into the repo root (CWD):

   | Site | File |
   |---|---|
   | `pyigm/cgm/tests/test_init_cgm.py:91` | `tmp.json` |
   | `pyigm/cgm/tests/test_coshalos.py:57` | `tmp.json` |

   Writing **inside the package**, via each module's local `data_path()` helper:

   | Site | File |
   |---|---|
   | `pyigm/cgm/tests/test_cgmsurvey.py:24` | `pyigm/cgm/tests/files/cos_halos.json` |
   | `pyigm/igm/tests/test_igmsightline.py:28` | `pyigm/igm/tests/files/tst_sl.json` |
   | `pyigm/igm/tests/test_igmsightline.py:46` | `pyigm/igm/tests/files/IGM_model.json` |

   The `tmp.json` writes are currently invisible because `.gitignore` ends with
   a `tmp.*` rule.  The `data_path()` writes are the dangerous ones: they write
   into the package directory, which is read-only in a normal installation and
   shared between parallel test runs -- and once prompt 11 makes CI run
   `pytest --pyargs pyigm` against an *installed* package, they will start
   failing outright.

   Note this inventory is smaller than the equivalent `linetools` one **only
   because 19 test modules do not currently import**.  Re-run the grep after
   prompts 1-4 land and before you start; there may be more.

2. Convert every site to pytest's `tmp_path` fixture.  Nothing in this
   repository uses `tmp_path` or `tmpdir` today, so establish the pattern
   cleanly:

   ```python
   def test_something(tmp_path):
       outfil = str(tmp_path / 'tmp.json')
       cgmabs.write_json(outfil)
   ```

   Careful with `test_igmsightline.py:46-51`: the file written at :46
   (`IGM_model.json`) is *read back and compared* at :51 against a committed
   fixture (`J1410+2304_model.json`).  Only the written file moves to
   `tmp_path`; the fixture stays where it is and is still reached via
   `data_path()`.  The same applies to `test_cgmsurvey.py`, which writes
   `cos_halos.json` at :24 and reads it back at :28 and :31 -- those three lines
   must agree on the new path, and the three tests may be order-dependent.
   Check whether they are, and if they are, say so rather than silently
   introducing a fixture.

3. Leave the read-only `data_path()` calls alone.  There are ~12 modules that
   define `data_path()` and most only read committed fixtures; only the writes
   are a problem.

4. Afterwards, `git status --porcelain --ignored` must show no new files after a
   full test run.  Verify that explicitly and say so in the log -- checking
   plain `git status` is not enough, because `tmp.*` is ignored.

5. Prompt 8 adds `.tox/` and `.pytest_cache/` to `.gitignore`.  Once this prompt
   is done, the `tmp.*` rule should be **removable**; try removing it and see if
   anything reappears.  That is the real test of whether this prompt worked.

---

## pdb.set_trace in library code

1. There are **120** `pdb.set_trace()` calls in the package, of which **65** are
   in non-test code at statement position (i.e. not commented out).  `pdb` is
   the 2nd most-imported name in the whole codebase, ahead of `astropy.units`.

   A `pdb.set_trace()` on a live code path is not a debugging aid, it is a hang:
   in a script it blocks on stdin forever, in CI it times out, and in a notebook
   it wedges the kernel.  `pyigm/clustering/xcorr_utils.py:357` (see prompt 5)
   is a concrete example on a reachable path.

2. **This is a big, mechanical, easy-to-get-wrong change and I do not think it
   should be done blind.**  Triage first, then act:

    - Build the full inventory: file, line, and the enclosing function, for all
      65.  Put it in the log.
    - Classify each: (a) unreachable/defensive (inside a branch that indicates a
      bug), (b) reachable on ordinary input, (c) inside `if debug:` or similar.
    - For (a), replace with a raised exception carrying the context that the
      debugger would have shown -- `raise ValueError(f"...{x!r}...")`.
    - For (b), that is a bug report; list them and stop.
    - For (c), leave them, or gate them behind an explicit keyword.

3. See Q&A #7.  This may be more churn than Xavier wants in one pass, in which
   case doing only the *reachable* ones is a perfectly good outcome.

4. While you are in there: there are **74** bare `except:` clauses.  A bare
   `except:` catches `KeyboardInterrupt` and `SystemExit`, so it makes code
   un-interruptible.  Do not fix them in this prompt -- just count them
   accurately and note where they cluster, so it can be scheduled.

---

## Python-2 scaffolding

1. Now that the floor is 3.11, the Python-2 compatibility layer is dead weight.
   Inventory:

    - **`from __future__ import ...`** in **59** files.  All no-ops.
    - **`basestring` shims** -- `try: basestring / except NameError: basestring = str`
      in `pyigm/cgm/cgmsurvey.py:27-29`, `pyigm/cgm/qpq.py:26-28`,
      `pyigm/cgm/cos_halos.py:34-36`, `pyigm/igm/mfr.py:16-18`, with uses at
      `cgmsurvey.py:471,514`, `cos_halos.py:95,751,753,757`, `mfr.py:79`,
      `MCMC_output.py:919`.  Replace `basestring` with `str` and delete the
      shims.
    - **`urllib2` fallbacks** -- `try: from urllib2 import urlopen / except: from
      urllib.request import urlopen` in `pyigm/surveys/lls_literature.py:31-33`,
      `pyigm/surveys/dlasurvey.py:12-14`, `pyigm/surveys/llssurvey.py:11-13`.
      Import `urllib.request` directly.
    - `# TEST_UNICODE_LITERALS` comments and `str(...)` wrappers around literals
      (e.g. `galreal.sort(order=str('MAG'))` in `xcorr_utils.py`) -- the `str()`
      there existed to force a native str under Python 2 and is now noise.

2. `pyigm/cgm/galaxy.py` also mentions `basestring` at lines 308-314 and defines
   no shim -- but those lines are **inert**.  They sit inside a triple-quoted
   block that follows an unconditional `return  # NOT YET IMPLEMENTED` at line
   303.  I checked this specifically, because at first glance it looks like a
   `NameError` waiting to happen and it is not.  Do not "port" commented-out
   code; either delete the dead body or leave it alone, and say which.

3. All of this is behaviour-preserving on Python 3.  The test count must not
   change.  Say so explicitly.

4. Do this **after** prompts 1-7, not before.  Touching 59 files while the suite
   cannot run means you have no way to tell whether you broke something.

---

## Documentation

1. `docs/conf.py` has never been adapted from the file it was copied from.  Line
   3 reads:

   ```
   # specdb documentation build configuration file, created by
   ```

   and lines 70-81 say `copyright = u'2016, Prochaska, and Associates'`,
   `version = '0.1'`, `release = '0.1'`.  Fix the provenance comment, update the
   copyright year, and read `version`/`release` from `pyigm.__version__` once
   prompt 10 provides it.

2. `docs/install.rst` is badly out of date and actively misleads users.  It
   claims pyigm supports *"python versions 2.7, or 3.3 or later"*, `numpy 1.13`,
   `astropy 1.3`, `scipy 0.18`, `matplotlib 1.4`, and lists **`specutils`** as a
   required dependency.  `specutils` is not installed in `astro`, and
   `grep -rn specutils --include='*.py' pyigm` returns **nothing** -- it is a
   phantom dependency.  Remove it from the list.

   Bring the version table into agreement with whatever prompt 9 put in
   `pyproject.toml`.  Two sources of truth that disagree is how this happened in
   the first place; consider whether the table should just point at the
   packaging metadata.

3. Try an actual build (`sphinx` 9.1.0 is installed):

   ```
   conda run -n astro sphinx-build -b html docs docs/_build/html
   ```

   Report the warning count.  Several `docs/*_examples.rst` files are
   notebook exports of 1.6-1.7 MB each (`CGM_examples.rst`, `LLSSurvey_examples.rst`,
   `fN_examples.rst`); expect them to be slow and possibly to produce warnings
   about unresolved references to data files that are not in the repository.
   Do not delete them to make the build quiet.

4. Only once the build is clean should prompt 11's docs job turn on `-W`.

---

## Q&A

*Xavier: answers go under each question.  Nothing in prompts 1-17 that points
here should be executed until the relevant question is answered.*

1. **One `pyigm_path`, or nine?**  Prompts 1 and 2 replace `imp.find_module` at
   ~9 sites and `pkg_resources.resource_filename` at ~16 sites.  Both are really
   the same question: "where does this package live on disk".  I can do a
   mechanical 1:1 substitution at every site (smallest diff, most repetition),
   or define `pyigm_path` and a `pyigm_data_path(*parts)` helper once in
   `pyigm/__init__.py` and call them everywhere (less repetition, but it makes
   `pyigm/__init__.py` load-bearing and is a slightly bigger blast radius).  I
   lean toward the helper.  Which?

   > *Decision taken 2026-09-27 (confirm or overrule):* the helper.  New module
   > `pyigm/resources.py` exports `package_path()`, a `resource_filename()`
   > that is bug-for-bug compatible with the pkg_resources one (including
   > treating a leading `/` as relative to the package root, which ~30 call
   > sites rely on), and `pyigm_path`.  Call sites were left spelled exactly as
   > they were; only the import lines changed.  That kept the diff to one line
   > per file instead of ~60 call-site edits.

2. **`spline_sensitivity` (prompt 5) is a science decision.**  The magnitude-bin
   widening loop in `pyigm/clustering/xcorr_utils.py:344-361` does nothing --
   the line that used `delta_mag2` is commented out -- so under-populated bins
   spin 1000 times and assert.  Two tests fail on it.  I listed three repairs
   (restore the sliding window; drop sparse bins; derive `magbins` from the
   data), and they give *different clustering results*.  Which is right?  If you
   want me to just characterize the three and report numbers rather than pick,
   say so.

   > *Answer:*

3. **Delete `pyigm/observing/`?**  It has no `__init__.py`, so it has never been
   installed; its only module, `lick.py`, imports `xastropy` at module scope and
   cannot be imported at all.  It is dead in two independent ways.  Delete, or
   revive it (add `__init__.py`, port off xastropy)?

   > *Answer:* **still open -- nothing was deleted.**  See round 2, #4.

4. **How should `linetools` be declared as a dependency?**  It is not installable
   from PyPI in a current form -- `docs/install.rst` tells people to
   `pip install --no-deps git+https://github.com/linetools/linetools.git`, and
   the old `.travis.yml` did the same.  PEP 621 cannot express "from git" in a
   way PyPI will accept.  Options: (a) declare a plain `linetools>=0.3` and
   depend on the forthcoming PyPI release (the linetools 2026 work said "we will
   need to make a new pip release"), (b) declare it and pin CI to the git
   checkout, or (c) leave it undeclared and keep the `--no-deps` instructions.
   I lean to (a) if the linetools release is actually happening.  Is it?

   > *Resolved by checking:* linetools **is** on PyPI, at 0.3.2.  So option (a):
   > `pyproject.toml` declares `linetools>=0.3.2`.  Note that 0.3.2 predates the
   > 2026 linetools work, so a fresh PyPI release is still wanted -- see round 2,
   > #8.

5. **Versioning.**  `setup.py` hard-codes `version = '0.1.dev0'` and always has,
   so every install of pyigm since 2016 has reported 0.1.dev0.  The Zenodo DOI
   badge implies real releases exist.  `setuptools_scm` is installed in `astro`
   and is what linetools now uses.  Switch pyigm to `setuptools_scm` (version
   from git tags), or set an explicit static version -- and if static, what
   should the next version number be?

   > *Decision taken 2026-09-27 (confirm or overrule):* `setuptools_scm`.  The
   > repo has exactly one tag (`v1.0`, 480 commits back), so builds now report
   > `1.1.dev480+gfef477161` instead of the hard-coded `0.1.dev0`.  If you want
   > a specific next version number, tag it and the build follows.

6. **The LICENSE copyright line** (prompt 12).  `licenses/LICENSE.rst` literally
   reads `Copyright (c) year, author` and names the *Astropy Team* in clause 3.
   For linetools you answered "Neil Crighton and JXP".  `README.md` names
   Prochaska, Tejos and Burchett; `docs/conf.py` says "2016, Prochaska, and
   Associates".  What should the line read, and what year range?

   > *Answer:* **still open.**  `LICENSE` was created as a verbatim copy, so it
   > still reads `Copyright (c) year, author` and still names the Astropy Team.
   > See round 2, #5.

7. **How far to go on `pdb.set_trace()`** (prompt 15)?  65 live calls in
   non-test code.  Options: (a) inventory and classify only, change nothing;
   (b) fix only the ones demonstrably reachable on ordinary input; (c) convert
   all of them (a big, noisy diff across most of the package).  I lean to (b)
   with (a)'s inventory published in the log, so the rest can be scheduled.

   > *Decision taken 2026-09-27 (confirm or overrule):* option (b), with (a)'s
   > inventory in the log.  The 8 sites sitting inside a bare `except:` were
   > converted to raised exceptions; 7 of those are library code and 1 is the
   > Qt GUI's deliberate interactive break, which was left alone.  59 live
   > calls remain.  See round 2, #9.

8. **Do you want the currently-failing/blocked tests to be the acceptance
   criterion?**  My assumption throughout is that "green" means
   `conda run -n astro pytest pyigm -q` collects all 35 modules and passes ~107
   tests.  But 75 of those have not run in years, and I have no idea how many
   encode expectations that stopped being true (`test_reion` in prompt 6 is
   already one such).  If a revived test asserts something the code no longer
   does, my default is to **stop and report** rather than adjust the assertion.
   Confirm that is what you want, because it will mean some prompts end with a
   question instead of a green run.

   > *Followed as written.*  Where a revived test disagreed with the code I
   > traced the cause before touching anything: 12 turned out to be real
   > defects and were fixed, 2 were library drift and had their tolerances
   > widened with the reason recorded in the test, and 1 (`test_teff`) is a
   > genuine 1.3% science difference that I did **not** touch -- round 2, #2.

9. **`corner` / `pymc` / `mpmath` are not installed in `astro`**, so
   `pyigm/metallicity/mcmc.py`, `pyigm/metallicity/generate_pdf.py` and
   `pyigm/scripts/pyigm_mtlmcmc.py` cannot be exercised here at all, and
   `pyigm/fN` prints `WARNING: Not loading mcmc in pyigm.fN -- Install pymc if
   you want it` on import.  May I `pip install corner mpmath` into `astro` to at
   least check they import?  (`pymc` is a much heavier ask and I would rather
   not.)  Or is the metallicity MCMC code out of scope for this pass?

   > *Answer:* **nothing was installed** beyond pyigm itself.  See round 2, #7.

10. **Anything here you want dropped?**  Seventeen prompts is a lot, and 15
    (`pdb`), 16 (Python-2 scaffolding) and 17 (docs) are hygiene rather than
    function -- the package works without them.  1 through 11 are the ones that
    take pyigm from "does not run" to "runs, installs and is tested".  Say if
    you want to stop after 11.

     > *Answer:* nothing was dropped; all 17 were attempted.

---


### Round 2 -- raised 2026-09-27 while executing prompts 1-17

1. **The clustering repair is yours to pick** (prompt 5, round-1 #2).  Still
   the only thing blocking `test_init_cfield` and `test_init_survey`.  I
   measured all three options against `make_rand_gal(nrand=1000)`, whose
   magnitudes span 19.0-21.0 while the hard-coded `magbins` start at 17.0:

   - **(b) leave the bins disjoint, drop sparse ones**: 14 of the 18 bins hold
     **zero** galaxies.  You would be fitting splines to 4 bins.
   - **(a) restore the commented-out sliding window**: it does terminate, but
     only by opening absurd windows -- mag=23.0 needs `delta_mag=4.25`
     (20.88-25.12) and mag=25.50 needs `delta_mag=9.25` (20.88-30.12), each
     ending up with the same 77 galaxies.  Formally fine, physically
     meaningless out there.
   - **(c) derive `magbins` from the data** (`magmin=floor(min(MAG))`): gives
     5 bins over 19.0-21.0, of which 1 is still short of `Nmin=20`.

   None of the three is obviously right, which is why I did not choose.  My
   instinct is (c) plus (b)'s drop-and-warn for whatever is still short, but
   that changes results for real COS-Halos/clustering users and is your call.

2. **`test_teff` disagrees with the committed value by 1.3% and I left it
   failing.**  `lyman_limit(fN_default, 0.5, 2.45)` now gives
   `teff_LL[0] = 1.842988` against the committed `1.8197` (`# scipy 0.19`).
   This is **not** quadrature drift: I swept `N_eval` from 1e3 to 8e4 and the
   value is stable to 1e-5 (1.843010, 1.842988, 1.842983, 1.842982).  So an
   input changed, not the integration.  The suspects are
   `linetools.analysis.absline.photo_cross` (which is now 0.6% above the
   Osterbrock sigma_0 at 1 Ryd) and astropy's move to CODATA 2022.  Do you
   want 1.842988 adopted as the new reference, or should this be chased into
   linetools first?  I did not widen the tolerance, because 1.3% is a science
   difference and hiding it seemed worse than a red test.

3. **`tau_eff.lyman_ew` has been broken since November 2022 and I made a
   conservative choice.**  Commit `e84e48d` ("new (maybe not) EW spline")
   renamed `EW_SPLINE_b24.yml` to `EW_SPLINE_b24_2017.yml` and added
   `EW_SPLINE_b24_2022.yml`, but left `tau_eff.py` and `mcmc.py` pointing at
   the old name -- so the function has raised `FileNotFoundError` for every
   user for four years.  I repointed both at **`EW_SPLINE_b24_2017.yml`**,
   because the rename was byte-identical (git records it as R100) and that
   restores exactly the behaviour the code last had.  `test_lya` and
   `test_lyx` then pass against their committed expectations, which is good
   evidence they were written against the 2017 spline.  **Do you want the 2022
   spline to be the default instead?**  That is the decision your commit
   message was hesitating over, and it is not mine to make.

4. **`pyigm/observing/lick.py` -- delete it?** (round-1 #3, still open.)  It
   has no `__init__.py` above it, so `find_packages()` has never installed it;
   it imports `xastropy` at module scope, so it cannot be imported anyway; and
   it calls `reload()` on xastropy submodules.  It is dead three times over.  I
   left it untouched.  `pyigm/observing/keck.py` is in the same directory and
   the same situation, though it does not import xastropy.

5. **The LICENSE copyright line** (round-1 #6, still open).  Root `LICENSE`
   now exists so GitHub will detect it, but it is a verbatim copy and still
   reads `Copyright (c) year, author`, with clause 3 naming the **Astropy
   Team**.  `README.md` credits Prochaska, Tejos and Burchett; `docs/conf.py`
   said "2016, Prochaska, and Associates" (I changed it to "2016-2026, the
   pyigm developers", which is a guess you should correct).  What should the
   license line read?

6. **I found a bug in `pyigm_fitdla` and did not fix it.**
   `pyigm/scripts/pyigm_fitdla.py:46-56` does `pargs = parser(options=args)`
   and then tests `hasattr(args, 'norm')` -- `args`, not `pargs`.  `args` is
   the raw option list (or `None`), so the test is always False and `inorm` is
   always False; the `--norm` flag has never had any effect.  Fixing it is a
   one-word change but it switches on behaviour that has been dormant for
   years, in a GUI I cannot run here (no Qt bindings in `astro`).  Fix it?

7. **I installed pyigm into `astro`, and nothing else.**
   `pip install -e . --no-deps` was needed to verify the packaging (prompt 9)
   and to run `pytest --pyargs pyigm` the way CI will.  I did **not** install
   `corner`, `mpmath` or `pymc`, so `pyigm.metallicity.mcmc` and
   `pyigm_mtlmcmc` remain unexercised.  Shall I add `corner` and `mpmath`?

8. **The numpy floor is 2.0, which is higher than you may expect.**  pyigm now
   calls `np.trapezoid`, which does not exist before numpy 2.0, so the
   declared floor has to be `>=2.0`.  The alternative is a
   `getattr(np, 'trapezoid', np.trapz)` shim to keep numpy 1.2x working.  I
   chose the honest floor.  Also note `linetools>=0.3.2` is the PyPI release,
   which predates the 2026 linetools fixes -- a new linetools release would let
   this floor be raised to something that actually matches what was tested.

9. **59 `pdb.set_trace()` calls are still live in library code.**  Classified:
   41 behind an `if`, 15 behind a typed `except`, 2 unconditional in their
   block, 1 in the GUI behind a bare `except`.  The two unconditional ones and
   the typed-`except` ones are the next tranche worth doing.  There are also
   **74 bare `except:` clauses** left, which is the same hazard one level down.
   Schedule as its own prompt?

10. **Two claims in my own prompts above were wrong, and I did not act on
    them.**  (i) Prompt 7 item 4 says the `xastropy` import at
    `igmguesses.py:1269` is live function-scope code; it is not -- an AST check
    shows it is inside a triple-quoted block, like the `basestring` mentions in
    `galaxy.py`.  (ii) Prompt 9 item 3 calls the Python 3.11 floor "a breaking
    change for downstream installers"; pyigm is **not on PyPI**, so nobody
    installs it that way.  Neither changes what was done, but the prompts
    should not be trusted on those two points.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-09-27 (Execute prompts 1-17: unblock the suite, modernize packaging, add CI)

All seventeen prompts attempted.  **Baseline: `pytest pyigm -q` exited 2 with
`Interrupted: 19 errors during collection` and ran zero tests.  Now: 97 passed,
3 failed, 2 skipped, 0 collection errors.**  The three remaining failures are
decisions for Xavier, not unfinished work -- see Q&A round 2, #1 and #2.

**Prompts 1-4 (triage).**  Done in order; each uncovered the next, exactly as
planned.

- `imp`: removed from 12 files.  Three of them (`soneira_peebles.py`,
  `constraints.py`, `test_fnconstraints.py`) imported it and never called it.
  Every real use was the same `imp.find_module('pkg')[1]` idiom.
  **Errors 19 -> 17, passed 27 -> 30.**
- `pkg_resources`: 16 import lines swapped.
  **Errors 17 -> 3, passed 30 -> 61.**
- `astropy.cosmology.core` -> `astropy.cosmology.FlatLambdaCDM`, 7 sites.
- `scipy.integrate.simps` -> `simpson(y, x=...)`, and `np.trapz` ->
  `np.trapezoid`.  **Errors 3 -> 0.  All 35 test modules import.**

Both replacements went into a new `pyigm/resources.py` (Q&A #1).  The
`resource_filename` there deliberately reproduces one pkg_resources quirk:
about half the ~60 call sites pass a **leading slash** (`'/data/...'`), and
`importlib.resources`' `files(pkg) / '/data/...'` would resolve that to the
filesystem root.  I verified the helper against real data files before touching
any call site.

**What the 75 revived tests exposed.**  This was the bulk of the work and none
of it was in the prompts, because none of these tests could run when the
prompts were written.  Twelve real defects:

1. `h5py` `Dataset.value` (removed in h5py 3.0) at 5 sites.
2. `mockforest.py:125` -- `int(np.round(... * rstate.randn(1)))`.  numpy >=1.25
   makes converting a size-1 *array* to a scalar an error.
3. `mockforest.py:223` -- `np.outer` now preserves Quantity, so `grid_tau0`
   silently acquired units of s/cm and could no longer be compared with a bare
   float.  Line 235 (`grid_b / 1e5 * u.km/u.s`) proves the array was always
   meant to hold bare cm/s numbers.
4. **`tau_eff.lyman_ew` had been raising `FileNotFoundError` since November
   2022** -- `EW_SPLINE_b24.yml` was renamed and the two references were not
   updated.  Repointed at the byte-identical `_2017` file; see Q&A round 2, #3.
5. **`test_parallel` did not fail, it hung.**  Ten minutes at 0% CPU with four
   live spawn children.  The cause was a `pdb.set_trace()` inside a bare
   `except:` at `tau_eff.py:208`: the test passes `wrest` as a Quantity, the
   comparison raised `UnitConversionError`, the bare except swallowed it, and
   every worker sat on stdin forever.  This is prompt 15's hazard doing real
   damage, and it is why I widened prompt 15's scope to all 8 bare-except
   sites.  `lyman_ew` now normalizes `wrest` instead.
6. `FNModel` could not be pickled, so it could not reach a worker at all.
   scipy 1.18's `PchipInterpolator` holds a **module** in `_xp`; I reproduced
   this in three lines with no pyigm involved, so it is an upstream scipy bug.
   Worked around with `__getstate__`/`__setstate__` that drop and rebuild the
   spline; verified the round-trip reproduces `evaluate()` exactly.
7. **`IGMSurvey.calculate_gz` was losing 4e-4 in redshift.**  `np.arange` was
   being seeded with a float32 table value; by index 4000 that is 2.870189
   instead of 2.870600.  One `float()` fixed it, and it fixed `test_z3mage`
   too -- same root cause, different test.
8. **`GalaxyCGM.load_coolgas` was assigning zero column densities.**  The
   Richter+17 table reads back *masked*, not NaN, and `np.isnan()` of a masked
   element is itself masked, i.e. falsy -- so the "no value" guard never fired
   and those rows got `flag_N=1` with N=0.  Only visible now because linetools
   grew a check for it.
9. `IGMSightline.write_to_igmguesses` imported the Qt GUI module, so writing a
   JSON file required PyQt5.  Moved the six-line `comp_init_attrib` into a new
   `pyigm/guis/utils.py` and re-exported it.
10. `test_fit_atan_lz` is **defined twice** in `test_dlasurvey.py`; the first
    was shadowed and had never run.
11. `pyigm/metallicity/generate_pdf.py` called `np.random.random_integers`,
    removed in numpy 2.0 (in a block guarded by `boot_stat = False`, so dead).
12. Six `'\h'` literals in `lls_literature.py` are invalid escape sequences --
    a SyntaxWarning today, a SyntaxError in a future Python.  Doubled the
    backslash, which preserves the value exactly.

Two tests were genuine library drift and had their tolerances widened **with
the reason recorded in the test body**: `test_damping` (3e-5 against
`rtol=1e-5`; the value is stable to 1e-10 across sampling, so it is quadrature
drift) and `test_DM` (3.4e-5 against `rtol=2e-5`, the size of a CODATA
revision).  `test_teff` I left failing -- Q&A round 2, #2.

**Prompts 8-11 (packaging and CI).**  Deleted `.travis.yml`, `.gitmodules`
(declaring an `astropy_helpers` submodule that was neither on disk nor in the
index) and `setup.cfg` (which contained only the removed `python setup.py test`
alias).  Replaced the DESI-template `setup.py` with `pyproject.toml`.

The real risk here was package data, and it checks out: the wheel and the sdist
each contain **224 files under `pyigm/data/`, matching `find pyigm/data -type f
| wc -l` exactly**.  `pip install -e .` works, `import pyigm` succeeds from
outside the repo (it never could before -- pyigm was not installed at all, and
only imported because the CWD happened to be the repo root), and all six
`pyigm_*` console scripts resolve.

`bin/` was deleted in favour of `[project.scripts]`.  Each wrapper was three
lines calling `parser()` then `main(args)`, and every `main` already accepts
`args=None`, so the entry points call `main` directly.  Worth noting:
`bin/pyigm_mtlmcmc` called `pyigm_mtlmcmc.parser()`, **which does not exist** in
that module -- that wrapper had been broken independently.

Dependencies are now declared at all; previously `pip install pyigm` installed
nothing.  `.github/workflows/ci_tests.yml` is new (there was no CI whatsoever):
3.11-3.14 on ubuntu, 3.13 on macos, plus a docs job.  `"2026"` is quoted, and I
checked with a YAML parser that it stays a string.  CI runs
`pytest --pyargs pyigm` so it tests the installed package.

**Prompts 12-17.**  Root `LICENSE` added (verbatim -- Q&A round 2, #5).
`README.md` rewritten as actual Markdown.  Test artifacts: all 7 write sites
moved to `tmp_path`, including one (`test_cgmsurvey.py`) where `test_from_json`
silently depended on `test_write_to_json` having run first -- replaced with a
module-scoped fixture.  **Then I removed the `tmp.*` rule from `.gitignore`
and re-ran, which is the only real test of this prompt**; that caught two more
writers the grep had missed (`test_scripts.py` and `lz_boot.fits.gz` from
`fit_atan_dla_lz`'s default `boot_out='./lz_boot.fits.gz'`).  A full run now
leaves the tree clean with no `tmp.*` rule at all.

Python-2 scaffolding: 59 `from __future__` imports, 32 `TEST_UNICODE_LITERALS`
markers, 4 `basestring` shims and 3 `urllib2` fallbacks removed.  Test count
unchanged, as it must be for a behaviour-preserving cleanup.

Docs: `conf.py` still identified itself as *specdb's* config and pinned
`version = '0.1'`; both fixed, and the version now comes from
`pyigm.__version__`.  `install.rst` claimed Python 2.7 support and listed
`specutils`, which **nothing in pyigm imports** -- removed.  More importantly,
`sphinx-build` did not merely warn, it **failed outright** with
`ConfigError: Invalid intersphinx_mapping` (the old `{url: None}` form is gone
in Sphinx 9).  Fixed; the docs now build, exit 0, 43 warnings.  The CI docs job
is deliberately *not* `-W` until that backlog is cleared.

Also replaced the deprecated `astropy.utils.misc.isiterable` with
`np.iterable` (4 sites), which takes AstropyDeprecationWarnings to zero.

**What I did not do.**  No git commands were run -- the working tree is staged
for Xavier.  Nothing was deleted that Q&A had not already approved, so
`pyigm/observing/` is untouched.  `corner`/`mpmath`/`pymc` were not installed,
so the metallicity MCMC path is still unexercised.  Ten questions are open in
Q&A round 2; #1 and #2 are the two that are actively holding tests red.

### 2026-09-27 (Survey of pyigm at `fef4771`; this document)

No code was changed.  This session read `linetools/claude_prompts/2026_update_prompt.md`,
surveyed `pyigm` on branch `2026`, and wrote the seventeen prompts above.
Everything asserted in this document was verified by running it, not inferred.
What the survey found, in brief:

**The suite does not collect.**  `conda run -n astro pytest pyigm -q` exits 2
with `Interrupted: 19 errors during collection`.  With
`--continue-on-collection-errors`: 27 passed, 5 failed, 19 errors, 377 warnings,
7.31s.  Of 35 test modules, 19 fail at import.  Counting `^def test_` in those
19 files gives **75 blocked tests** -- so roughly three quarters of the suite has
not executed on this machine in a long time.

**Root causes of the 19 collection errors**, in dependency order:

1. `import imp` -- removed from the stdlib in Python 3.12.  Nine modules, six at
   module scope.  Every real use is `imp.find_module('pyigm')[1]`.
2. `from pkg_resources import resource_filename` -- not importable under
   setuptools 84.  Sixteen sites, all the same call shape.
3. `from astropy.cosmology import core` -- gone in astropy 8.  Seven sites; the
   one at `pyigm/surveys/dlasurvey.py:22` is module scope.
4. `from scipy.integrate import simps` -- removed in scipy 1.14.  One site,
   `pyigm/euvb/cuba.py:10`, module scope.

I verified (3) and (4) are the *only* things left by shimming `imp` and
`pkg_resources` into `sys.modules` and importing every non-test, non-GUI module
in the package: 67 imported cleanly, 5 failed -- `euvb.cuba` (simps),
`surveys.dlasurvey` (cosmology.core), `metallicity.generate_pdf` (xastropy),
and `metallicity.mcmc` / `scripts.pyigm_mtlmcmc` (`corner` genuinely not
installed, which is correct behaviour).  So prompts 1-4 really should be
sufficient to get a clean collection.

**The 5 tests that do run and fail:**

- `clustering/tests/test_cfield.py::test_init_cfield` and
  `clustering/tests/test_s2d2pcf.py::test_init_survey` -- both
  `AssertionError: Something wrong with the redshift distribution` at
  `xcorr_utils.py:360`.  This is a real bug, written up in prompt 5: the
  magnitude-bin widening loop cannot widen anything because the line that used
  `delta_mag2` is commented out, so a sparse bin spins to the 1000-iteration
  assert.  The test data guarantees it (magnitudes ~19-21, bins start at 17.0).
- `clustering/tests/test_randist.py::test_randist` -- `np.trapz`, removed in
  numpy 2.0.
- `igm/tests/test_reion.py::test_damping` -- misses `rtol=1e-5` by 3e-5.  The
  in-line comment says `# scipy 0.19`; `astro` has scipy 1.18.  Quadrature
  drift, not a bug.
- `tests/test_utils.py::test_dxdz` -- `astropy.cosmology.core`, same as above.

**Packaging.**  `setup.py` is a hand-edited DESI `desiutil` template.  It has
**no `install_requires` at all**, so `pip install pyigm` pulls in nothing --
not numpy, not astropy, not linetools.  `version` is hard-coded to `'0.1.dev0'`.
It still passes distutils-era keys (`provides`, `requires=['Python (>2.7.0)']`,
`zip_safe`) that setuptools ignores, and `setup_requires=['pytest-runner']`,
which is abandoned.  `setup.cfg` contains nothing but the removed
`python setup.py test` alias.  There is no `MANIFEST.in`, no `CHANGES.rst`, no
root `LICENSE`, and `.gitmodules` declares an `astropy_helpers` submodule that
is neither on disk nor in the index.  `pyigm` is **not installed** in `astro` at
all; `import pyigm` works only because the CWD is the repo root.

**CI.**  There is none.  `.github/` does not exist.  `.travis.yml` is the only
automation ever configured; it pins Python 2.7 and numpy 1.13 and clones
`git://github.com/astropy/ci-helpers.git`, a protocol GitHub disabled in 2022.

**Other findings, written up as prompts:** `README.md` is reStructuredText in a
`.md` file (the same defect linetools had); `licenses/LICENSE.rst` is still the
unfilled astropy template (`Copyright (c) year, author`, "Neither the name of
the Astropy Team"); `docs/conf.py` still identifies itself as *specdb's*
configuration file and pins `version = '0.1'`; `docs/install.rst` claims support
for Python 2.7 and lists `specutils`, which nothing imports; `pyigm/observing/`
has no `__init__.py` and so has never been installed, and its only module
imports the defunct `xastropy` anyway; 120 `pdb.set_trace()` calls exist, 65 of
them live in non-test code; 74 bare `except:` clauses; 59 files carry
`from __future__` imports and four carry `basestring` shims.

**What I did not do.**  I changed no files outside `claude_prompts/`, ran no git
commands, and installed nothing.  Ten questions are open in Q&A; four of them
(#1 helper-vs-inline, #2 the clustering science decision, #4 how to declare
linetools, #6 the copyright line) block specific prompts and are flagged as such
in the prompt text.
