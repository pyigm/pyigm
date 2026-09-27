""" Locating files that ship inside installed packages

This module replaces two APIs that pyigm relied on for a decade and that no
longer exist in the environments it now supports:

* ``imp.find_module('pyigm')[1]`` -- ``imp`` was removed in Python 3.12
* ``pkg_resources.resource_filename`` -- no longer shipped by setuptools

Both are re-expressed here in terms of :mod:`importlib.resources`.
"""

import os

from importlib.resources import files


def package_path(package):
    """ Absolute path of the directory an installed package lives in

    Parameters
    ----------
    package : str
      Importable package name, e.g. 'pyigm'

    Returns
    -------
    path : str
      Directory holding the package, with no trailing separator.  This is the
      value ``imp.find_module(package)[1]`` used to return.
    """
    return str(files(package))


def resource_filename(package, resource):
    """ Absolute path of a file shipped inside a package

    Drop-in replacement for the removed ``pkg_resources.resource_filename``.
    As with the original, a leading '/' in `resource` is relative to the
    package root and not to the root of the filesystem.

    Parameters
    ----------
    package : str
      Importable package name, e.g. 'pyigm'
    resource : str
      Path of the file relative to the package directory, with '/' separators

    Returns
    -------
    path : str
    """
    relative = resource.replace('/', os.sep).lstrip(os.sep)
    return os.path.join(package_path(package), relative)


# Path of the pyigm package itself.  Many modules import this by name.
pyigm_path = package_path('pyigm')
