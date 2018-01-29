"""  Module for the QPQ survey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
from pkg_resources import resource_filename
from astropy.table import Table, Column
from pyigm.cgm.cgmsurvey import CGMAbsSurvey
import pyigm

try:
    basestring
except NameError:  # For Python 3
    basestring = str


def load_qpq(v):
    """ Reads table with data

    Parameters
    ----------
    v : int
      dataset (6,7,8)

    Returns
    -------
    qpqdata : Table with data
    """

    if v == 8:
        q8file = resource_filename('pyigm', 'data/CGM/QPQ/qpq8_all_measured.dat')
        qpqdata = Table.read(q8file, format='ascii')

    if v == 7:
        q7file = resource_filename('pyigm', 'data/CGM/QPQ/qpq7_pairs.fits.gz')
        qpqdata = Table.read(q7file)

    if v == 6:
        q6file = resource_filename('pyigm', 'data/CGM/QPQ/qpq6_final_cut_2015nov16.fits')
        qpqdata = Table.read(q6file)

    if v not in [6,7,8]:
        print('Please choose 6, 7, or 8')
        return

    return qpqdata


def qpq7tab(itrans,non0=False):
    """ Read table with QPQ7 transitions

    Parameters
    ----------
    itrans : str
      transition, e.g. 'HI 1215
    non0 : bool
      return table with only non-zero EWs

    Returns
    -------
    tab7 : Table with transition itrans
    """

    q7 = load_qpq(7)

    if itrans == 'HI 1215':
        tab7 = q7['QSO', 'R_PHYS', 'Z_FG', 'QSO_BG', 'FLG_EWLYA', 'EWLYA', 'SIG_EWLYA', 'FLG_NHI', 'NHI', 'SIG_NHI']
        tab7.rename_column('EWLYA', 'EW')
        tab7.rename_column('SIG_EWLYA', 'sig_EW')
        tab7.rename_column('FLG_EWLYA', 'flag_EW')
        tab7.rename_column('SIG_NHI', 'sig_NHI')
        tab7.rename_column('FLG_NHI', 'flag_NHI')
        tab7.rename_column('Z_FG', 'z')


        if non0:
            tab7 = tab7[ tab7['EW']>0 ]
            print('Loading data for {:n} of {:n} QSOs'.format(len(tab7),len(q7)))

    else:
        trans = ['OI 1302', 'SiII 1304', 'CII 1334', 'NII 1083',
                 'SiIII 1206', 'SiIV 1393', 'SiIV 1402', 'SiII 1526',
                 'CIV 1548', 'CIV 1550', 'NV 1238', 'MgII 2796',
                 'MgII 2803']
        jj = [j for j in range(len(trans)) if (itrans == trans[j])]
        if len(jj) == 0:
            print('No data for ',itrans)
            return
        else:
            tab7 = q7['QSO', 'R_PHYS', 'Z_FG', 'QSO_BG']
            tab7.rename_column('Z_FG', 'z')
            tab7['flag_EW'] = q7['FLG_METAL_EW'][:, jj[0]]
            tab7['WREST'] = q7['METAL_WREST'][:, jj[0]]
            tab7['EW_RAW'] = q7['RAW_METAL_EW'][:, jj[0]]
            tab7['EW'] = q7['METAL_EW'][:, jj[0]]
            tab7['sig_EW'] = q7['METAL_SIGEW'][:, jj[0]]
            tab7['S2N'] = q7['METAL_S2N'][:, jj[0]]
            tab7['VCEN'] = q7['METAL_VCEN'][:, jj[0]]
            tab7['flag_EYE'] = q7['FLG_METAL_EYE'][:, jj[0]]
            if non0:
                tab7 = tab7[tab7['EW'] > 0]
                print('Loading data for {:n} of {:n} QSOs'.format(len(tab7), len(q7)))

    return tab7


def qpq6tab(itrans,non0=False):
    """ Read table with QPQ6 transitions

    Parameters
    ----------
    itrans : str
      transition, e.g. 'HI 1215
    non0 : bool
      return table with only non-zero EWs

    Returns
    -------
    tab6 : Table with transition itrans
    """

    q6 = load_qpq(6)

    if itrans == 'HI 1215':
        tab6 = q6['QSO', 'R_PHYS', 'Z_FG', 'QSO_BG', 'FLG_EWLYA', 'EWLYA', 'SIG_EWLYA', 'FLG_NHI', 'NHI', 'SIG_NHI']
        tab6.rename_column('EWLYA', 'EW')
        tab6.rename_column('SIG_EWLYA', 'sig_EW')
        tab6.rename_column('FLG_EWLYA', 'flag_EW')
        tab6.rename_column('SIG_NHI', 'sig_NHI')
        tab6.rename_column('FLG_NHI', 'flag_NHI')
        tab6.rename_column('Z_FG', 'z')
        if non0:
            tab6 = tab6[ tab6['EW']>0 ]
            print('Loading data for {:n} of {:n} QSOs'.format(len(tab6),len(q6)))

        return tab6

    else:
        print('No data for {:s} in QPQ6'.format(itrans))




def qpq8tab(ion):
    """ Read table with QPQ8 transitions

    Parameters
    ----------
    ion : str
      ion, e.g. 'HI'; or wavelength of a transition, e.g. '1215'

    Returns
    -------
    tab8 : Table with transition itrans
    """

    q8 = load_qpq(8)
    q7 = load_qpq(6)
    colnames = q8.colnames
    if ion.isdigit():
        colew = 'E'
    else:
        colew = 'c'
    jj = [jname for jname in colnames if (((ion+colew) == jname[0:len(ion)+1]) | (jname in ['Pair','subsys']) )]
    tab8 = q8[jj]

    # match with QPQ7?
    ##
    to_match = False
    if to_match:
        match78 = []
        R7 = []
        z7 = []
        qsos7 = np.asarray(q7['QSO'])
        qsos7b = []
        for iqso in qsos7:
            if iqso[-1] == ' ':
                qsos7b.append(iqso[-11:-1])
            else:
                qsos7b.append(iqso[-10:])

        for iqso in q8['Pair']:
            jj = [j for j in range(len(qsos7b)) if (iqso == qsos7b[j])]
            if len(jj) != 1:
                print('Not unique matches')
                pdb.set_trace()
            else:
                jj = jj[0]
                match78.append(jj)
                R7.append(q7[jj]['R_PHYS'])
                z7.append(q7[jj]['Z_FG'])

        tab8['R'] = R7
        tab8['z'] = z7


    return tab8


