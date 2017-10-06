# Module to run tests on scripts

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np
import json

from .. import pyigm_mkigmsys
from pyigm.abssys.dla import DLASystem


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_mkigmsys():
    outfile = 'tmp.json'
    pargs = pyigm_mkigmsys.parser(['dla','3.0',outfile,'--NHI=20.5', '--zem=4.','--vlim=-232,300'])
    pyigm_mkigmsys.main(args=pargs)
    # Read as JSON
    with open(outfile,'r') as f:
        jdict = json.load(f)
    assert np.isclose(jdict['zabs'], 3.)
    # Now load as a system
    dla = DLASystem.from_dict(jdict)
    assert np.isclose(dla.NHI, 20.5)
