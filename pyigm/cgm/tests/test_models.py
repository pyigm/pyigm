# Module to run tests on GalaxyCGM

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from ..models import ModifiedNFW


def test_modified_NFW():
    mNFW = ModifiedNFW()

