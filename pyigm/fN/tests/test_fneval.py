# Module to run tests on FNModel

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from pyigm.fN import dla


def test_lX_dla():
    # z<2
    lX = dla.lX(1.5)
    assert isinstance(lX, float)
    assert np.isclose(lX, 0.)
    # sensible values
    z = np.arange(2.5, 4.5, 0.1)
    lX = dla.lX(z)
    assert len(lX) == len(z)
    np.testing.assert_allclose(lX[0], 0.06)





