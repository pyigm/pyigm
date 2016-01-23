.. highlight:: rest

************
Installation
************

Dependencies
============

pyigm depends on these packages:

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.10 or later
* `astropy <http://www.astropy.org>`_ version 1.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.16 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `specutils <https://github.com/astropy/specutils>`_ version 0.2 or later
* `linetools <https://github.com/linetools/linetools>`_ version 0.1 or later
* `h5py <http://www.h5py.org/>`_ version 2.5 or later

We strongly recommend that you use `Anaconda
<https://www.continuum.io/downloads>`_ to install them. With Anaconda
you can check for the presence and versions of the dependencies with::

  conda list "^python$|numpy|astropy$|scipy$|matplotlib|pyyaml|specutils"

If you're missing any, install them with (for example)::

  conda install astropy matplotlib

If their versions are too old, update them with (for example)::

  conda update astropy

Specutils can't be installed with conda; instead it needs to be
installed using `pip <https://pip.pypa.io/en/latest/>`_::
  
  pip install --no-deps specutils

This is also true for linetools.
If you aren't using Anaconda then all of the dependencies can also be
installed with pip.

The following packages are required only for MCMC analysis of
metallicity PDFs (e.g. Fumagalli+16):

* `corner <https://github.com/dfm/corner.py>`_ v0.3 or later
* `mpmath <http://www.mpmath.org/>`_ version 0.19 or later
* `emcee <http://http://dan.iel.fm/emcee/current/>`_ version 2.1 or later

Installing pyigm
================

If you plan to play around with the code and possibly contribute
changes, then follow the instructions in the section below,
:ref:`installsource`. Otherwise simply use::

    pip install --no-deps git+https://github.com/pyigm/pyigm.git

and you're done!


.. _installsource:

Installing pyigm from Source
============================

*I just want to play with the code*
-----------------------------------

Install the development version like this::

    git clone https://github.com/pyigm/pyigm.git
    cd pyigm
    python setup.py develop

Now you can easily make tweaks to the code, which are immediately
applied to your installed version (you'll have to reload the relevant
modules to see those changes in an existing Python session, though).

*I want to make a code contribution to pyigm*
---------------------------------------------

Fantastic! In that case, follow the `Astropy developer guidelines
<http://docs.astropy.org/en/stable/development/workflow/development_workflow.html>`_,
replacing every instance of `astropy` in those instructions with
`pyigm`. This will install a 'fork' of pyigm that you can use
to submit code changes to the main repository.


Running Tests
=============

If you install pyigm from source, then you can run tests to see
whether your installation works correctly. From the source directory
run::

    python setup.py test

This takes a couple of minutes to run. If you notice any failures,
we'd love you to report them on the `pyigm issue tracker
<http://github.com/pyigm/pyigm/issues>`_.


Building Documentation
======================

Only do this if you're a developer! If you want build the
documentation, you also need to install Sphinx (version 1.3+) and
astropy_helpers::

  conda install sphinx
  pip install astropy-helpers

If you want the generate inheritance diagrams in the docs then you
also need to install graphviz (`MacOSX
<http://www.graphviz.org/Download_macos.php>`_, `Ubuntu
<http://www.graphviz.org/Download_linux_ubuntu.php>`_), but this isn't
required. Once you've installed the dependencies, change to the
`/docs` directory under the source directory and run::

  make html

The documentation should now be in _build/html.
