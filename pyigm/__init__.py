""" pyigm: analysis of the Intergalactic and Circumgalactic Medium """

from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__ = _version("pyigm")
except PackageNotFoundError:  # not installed, e.g. running from a source tree
    __version__ = "unknown"
