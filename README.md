# pyigm

[![CI Tests](https://github.com/pyigm/pyigm/actions/workflows/ci_tests.yml/badge.svg)](https://github.com/pyigm/pyigm/actions/workflows/ci_tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1045480.svg)](https://doi.org/10.5281/zenodo.1045480)

A package for the analysis of the Intergalactic Medium and the Circumgalactic
Medium.

## Primary builders

- J. Xavier Prochaska
- N. Tejos
- J. Burchett

## Installation

pyigm requires Python 3.11 or later.

```
pip install git+https://github.com/pyigm/pyigm.git
```

For a development checkout:

```
git clone https://github.com/pyigm/pyigm.git
cd pyigm
pip install -e ".[test]"
pytest --pyargs pyigm
```

Optional extras: `mcmc` (metallicity MCMC), `gui` (the Qt GUIs), `docs`.

## Documentation

<https://pyigm.readthedocs.io/>

## Related

pyigm is built on [linetools](https://github.com/linetools/linetools), which
provides the line lists, spectra and absorption-system machinery it depends on.
