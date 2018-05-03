"""
Module to run tests on pyigm.cgm.qpq
"""


import numpy as np
from pyigm.cgm.qpq import QPQ5,QPQ6,QPQ7,QPQ8
import pdb

def test_qpq():
    # Load classes
    qpq5 = QPQ5().survey_data
    qpq6 = QPQ6().survey_data
    qpq7 = QPQ7().survey_data
    qpq8 = QPQ8().survey_data
    assert len(qpq5.cgm_abs) == 73
    assert len(qpq6.cgm_abs) == 646
    assert len(qpq7.cgm_abs) == 427
    assert len(qpq8.cgm_abs) == 14

    # EW : QPQ6
    ew6 = qpq6.trans_tbl('HI 1215')['EW']
    i = np.where(ew6 > 0)[0]
    assert len(i) == 461
    flgs = qpq6.trans_tbl('HI 1215')['flag_EW']
    i = np.where(flgs == 1)[0]
    assert len(i) == 454

    ## EW : QPQ7
    ew7 = qpq7.trans_tbl('MgII 2803')['EW']
    i = np.where(ew7 > 0)[0]
    assert len(i) == 101
    flgs = qpq7.trans_tbl('MgII 2803')['flag_EW']
    i = np.where(flgs == 1)[0]
    assert len(i) == 20

    ## EW : QPQ8              # need to modify trans_tbl in order to use this
    #ew8 = qpq8.trans_tbl('AlII 1670')['EW']
    #i = np.where(ew8 > 0)[0]
    #assert len(i) == 19
    #flgs = qpq8.trans_tbl('AlII 1670')['flag_EW']
    #i = np.where(flgs == 1)[0]
    #assert len(i) == 16

