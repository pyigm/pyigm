"""  Module for the QPQ survey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
from pkg_resources import resource_filename

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent
from linetools import utils as ltu

from pyigm.cgm.cgmsurvey import CGMAbsSurvey
import pyigm
from pyigm.field.galaxy import Galaxy
from .cgm import CGMAbsSys
from pyigm.abssys.igmsys import IGMSystem

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

    if v == 8.5:  # additional information about QPQ8
        q8file = resource_filename('pyigm', 'data/CGM/QPQ/qpq8_pairs.fits')
        qpqdata = Table.read(q8file)

    if v == 7:
        q7file = resource_filename('pyigm', 'data/CGM/QPQ/qpq7_pairs.fits.gz')
        qpqdata = Table.read(q7file)

    if v == 6:
        q6file = resource_filename('pyigm', 'data/CGM/QPQ/qpq6_final_cut_2015nov16.fits')
        qpqdata = Table.read(q6file)

    if v not in [6,7,8,8.5]:
        print('Please choose 6, 7, 8, or 8.5 for QPQ6, QPQ7, QPQ8, additional information about QPQ8')
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

    # match with QPQ7? # J1420+1603 not found in QPQ 7 or 6
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



class QPQ6(CGMAbsSurvey):
    """Inherits CGM Abs Survey

        Parameters:
        -----------
        cdir : str, optional
          Path to the data
        """

    def __init__(self, cdir=None, load=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'QPQ6'
        self.ref = 'Prochaska+13'
        #
        if cdir is None:
            self.cdir = pyigm.__path__[0] + '/data/CGM/QPQ/'
        else:
            self.cdir = cdir
        #self.data_file = pyigm.__path__[0] + '/data/CGM/QPQ/'+'qpq6_final_cut_2015nov16.fits'
        self.data_file = self.cdir + 'qpq6_final_cut_2015nov16.fits'

        # Init?
        if load:
            self.load_data(**kwargs)


    def load_data(self, **kwargs):
        #
        #q6file = resource_filename('pyigm', 'data/CGM/QPQ/qpq6_final_cut_2015nov16.fits')
        q6file = self.data_file
        qpqdata = Table.read(q6file)
        nmax = len(qpqdata)   # max number of QSOs
        for i in range(nmax):
            # Instantiate the galaxy
            gal = Galaxy((qpqdata['RAD'][i],qpqdata['DECD'][i]), z=qpqdata['Z_FG'][i])
            gal.L_BOL = qpqdata['L_BOL'][i]
            gal.L_912 = qpqdata['L_912'][i]
            gal.G_UV = qpqdata['G_UV'][i]
            gal.flg_BOSS = qpqdata['FLG_BOSS'][i]
            gal.zsig = qpqdata['Z_FSIG'][i]*u.km/u.s

            # Instantiate the IGM System
            igm_sys = IGMSystem((qpqdata['RAD_BG'][i], qpqdata['DECD_BG'][i]),
                                qpqdata['Z_FG'][i], [-5500, 5500.] * u.km / u.s,
                                abs_type='CGM')   ## if velocity range lower - does not load all abslines
            igm_sys.zem = qpqdata['Z_BG'][i]
            igm_sys.NHI = qpqdata['NHI'][i]
            igm_sys.sig_NHI = qpqdata['SIG_NHI'][i]
            igm_sys.flag_NHI = qpqdata['FLG_NHI'][i]
            igm_sys.s2n_lya = qpqdata['S2N_LYA'][i]
            igm_sys.flg_othick = qpqdata['FLG_OTHICK'][i]
            igm_sys.z_lya = qpqdata['Z_LYA'][i]

            iname = qpqdata['QSO'][i]
            # Instantiate
            rho = qpqdata['R_PHYS'][i]*u.kpc
            cgabs = CGMAbsSys(gal, igm_sys, name=iname, rho = rho, **kwargs)
            aline = AbsLine(1215.67*u.AA, closest=True, z=igm_sys.zabs)
            aline.attrib['EW'] = qpqdata['EWLYA'][i]* u.AA   # Rest EW
            aline.attrib['sig_EW'] = qpqdata['SIG_EWLYA'][i] * u.AA
            if aline.attrib['EW'] > 3. * aline.attrib['sig_EW']:
                aline.attrib['flag_EW'] = 1
            else:
                aline.attrib['flag_EW'] = 3

            aline.attrib['coord'] = igm_sys.coord
            aline.limits._wvlim = qpqdata['WVMNX'][i]*u.AA
            dv = ltu.rel_vel(aline.limits._wvlim,aline.wrest*(1+qpqdata['Z_FG'][i]))
            aline.limits._vlim = dv

            abslines = []
            abslines.append(aline)
            ###
            comp = AbsComponent.from_abslines(abslines, chk_vel=False)

            # add ang_sep
            qsocoord = SkyCoord(ra=qpqdata['RAD'][i], dec=qpqdata['DECD'][i], unit='deg')
            bgcoord = SkyCoord(ra=qpqdata['RAD_BG'][i], dec=qpqdata['DECD_BG'][i], unit='deg')
            cgabs.ang_sep = qsocoord.separation(bgcoord).to('arcsec')

            cgabs.igm_sys.add_component(comp)
            self.cgm_abs.append(cgabs)

            #if i == 11:
            #    pdb.set_trace()

class QPQ7(CGMAbsSurvey):
    """Inherits CGM Abs Survey

        Parameters:
        -----------
        cdir : str, optional
          Path to the data
        """

    def __init__(self, cdir=None, load=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'QPQ7'
        self.ref = 'Prochaska+14'
        #
        if cdir is None:
            self.cdir = pyigm.__path__[0] + '/data/CGM/QPQ/'
        else:
            self.cdir = cdir
        #self.data_file = pyigm.__path__[0] + '/data/CGM/QPQ/'+'qpq7_pairs.fits.gz'
        self.data_file = self.cdir + 'qpq7_pairs.fits.gz'

        # Init?
        if load:
            self.load_data(**kwargs)


    def load_data(self, **kwargs):
        #
        #q7file = resource_filename('pyigm', 'data/CGM/QPQ/qpq7_pairs.fits.gz')
        q7file = self.data_file
        qpqdata = Table.read(q7file)
        nmax = len(qpqdata)   # max number of QSOs
        for i in range(nmax):
            # Instantiate the galaxy
            gal = Galaxy((qpqdata['RAD'][i],qpqdata['DECD'][i]), z=qpqdata['Z_FG'][i])
            gal.L_BOL = qpqdata['L_BOL'][i]
            gal.L_912 = qpqdata['L_912'][i]
            gal.G_UV = qpqdata['G_UV'][i]
            gal.flg_BOSS = qpqdata['FLG_BOSS'][i]
            gal.zsig = qpqdata['Z_FSIG'][i]*u.km/u.s

            # Instantiate the IGM System
            igm_sys = IGMSystem((qpqdata['RAD_BG'][i], qpqdata['DECD_BG'][i]),
                                qpqdata['Z_FG'][i], [-5500, 5500.] * u.km / u.s,
                                abs_type='CGM')
            igm_sys.zem = qpqdata['Z_BG'][i]
            igm_sys.NHI = qpqdata['NHI'][i]
            igm_sys.sig_NHI = qpqdata['SIG_NHI'][i]
            igm_sys.flag_NHI = qpqdata['FLG_NHI'][i]
            igm_sys.s2n_lya = qpqdata['S2N_LYA'][i]
            igm_sys.flg_othick = qpqdata['FLG_OTHICK'][i]
            igm_sys.z_lya = qpqdata['Z_LYA'][i]

            iname = qpqdata['QSO'][i]
            # Instantiate
            rho = qpqdata['R_PHYS'][i]*u.kpc
            cgabs = CGMAbsSys(gal, igm_sys, name=iname, rho = rho, **kwargs)
            aline = AbsLine(1215.67*u.AA, closest=True, z=igm_sys.zabs)
            aline.attrib['EW'] = qpqdata['EWLYA'][i]* u.AA   # Rest EW
            aline.attrib['sig_EW'] = qpqdata['SIG_EWLYA'][i] * u.AA
            if aline.attrib['EW'] > 3. * aline.attrib['sig_EW']:
                aline.attrib['flag_EW'] = 1
            else:
                aline.attrib['flag_EW'] = 3

            aline.attrib['coord'] = igm_sys.coord
            #aline.limits._wvlim = qpqdata['WVMNX'][i]*u.AA   ##   (no data in QPQ7 file)
            #dv = ltu.rel_vel(aline.limits._wvlim, aline.wrest * (1 + qpqdata['Z_FG'][i]))
            #aline.limits._vlim = dv

            abslines = []
            abslines.append(aline)
            comp = AbsComponent.from_abslines(abslines, chk_vel=False)
            cgabs.igm_sys.add_component(comp)

            # add metal lines
            for j in range(100):
                if qpqdata[i]['FLG_METAL_EW'][j] >0:
                    wave0 = qpqdata[i]['METAL_WREST'][j]
                    iline = AbsLine(wave0*u.AA, closest=True, z=igm_sys.zabs)
                    iline.attrib['EW'] = qpqdata['METAL_EW'][i][j]* u.AA   # Rest EW
                    iline.attrib['sig_EW'] = qpqdata['METAL_SIGEW'][i][j] * u.AA
                    iline.attrib['flag_EW'] = qpqdata['FLG_METAL_EW'][i][j]
                    iline.analy['flg_eye'] = qpqdata['FLG_METAL_EYE'][i][j]
                    iline.attrib['coord'] = igm_sys.coord
                    abslines = []
                    abslines.append(iline)
                    comp = AbsComponent.from_abslines(abslines, chk_vel=False)
                    cgabs.igm_sys.add_component(comp)

            # add ang_sep
            qsocoord = SkyCoord(ra=qpqdata['RAD'][i], dec=qpqdata['DECD'][i], unit='deg')
            bgcoord = SkyCoord(ra=qpqdata['RAD_BG'][i], dec=qpqdata['DECD_BG'][i], unit='deg')
            cgabs.ang_sep = qsocoord.separation(bgcoord).to('arcsec')


            self.cgm_abs.append(cgabs)





class QPQ8(CGMAbsSurvey):
    """Inherits CGM Abs Survey
        Contains information about different components

        Parameters:
        -----------
        cdir : str, optional
          Path to the data
        """

    def __init__(self, cdir=None, load=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'QPQ8'
        self.ref = 'Lau+16'
        #
        if cdir is None:
            self.cdir = pyigm.__path__[0] + '/data/CGM/QPQ/'
        else:
            self.cdir = cdir
        #self.data_file = pyigm.__path__[0] + '/data/CGM/QPQ/'+'qpq8_pairs.fits'
        self.data_file = self.cdir +'qpq8_all_measured.dat'

        # Init?
        if load:
            self.load_data(**kwargs)

    def load_data(self, **kwargs):
        #
        #q8file = resource_filename('pyigm', 'data/CGM/QPQ/qpq8_all_measured.dat')
        q8file = self.data_file
        qpqdata = Table.read(q8file,format='ascii')
        nmax = len(qpqdata)   # max number of QSOs
        q8filecoord = resource_filename('pyigm', 'data/CGM/QPQ/qpq8_pairs.fits')
        qpqdatacoord = Table.read(q8filecoord)

        # match names with qpqdatacoord
        qnames = []
        for i in range(len(qpqdatacoord)):
            qname = qpqdatacoord['QSO'][i].strip()
            qnames.append(qname[-10:])

        for i in range(nmax):
            # match names with qpqdatacoord
            j = np.where(np.asarray(qnames) == qpqdata['Pair'][i])[0][0]
            # Instantiate the galaxy
            gal = Galaxy((qpqdatacoord['RAD'][j],qpqdatacoord['DECD'][j]), z=qpqdatacoord['Z_FG'][j])
            gal.L_BOL = qpqdatacoord['L_BOL'][j]
            gal.L_912 = qpqdatacoord['L_912'][j]
            gal.G_UV = qpqdatacoord['G_UV'][j]
            gal.zsig = qpqdatacoord['Z_FSIG'][j]*u.km/u.s

            # Instantiate the IGM System
            igm_sys = IGMSystem((qpqdatacoord['RAD_BG'][j], qpqdatacoord['DECD_BG'][j]),
                                qpqdatacoord['Z_FG'][j], [-5500, 5500.] * u.km / u.s,
                                abs_type='CGM')
            # Redshifts: QSO emission redshifts
            igm_sys.zem = qpqdatacoord['Z_BG'][j]
            igm_sys.NHI = qpqdata['HIcol'][i]
            igm_sys.sig_NHI = [qpqdata['HIcolhierr'][i],qpqdata['HIcolloerr'][i]]
            igm_sys.s2n_lya = qpqdatacoord['S2N_LYA'][j]

            iname = qpqdata['Pair'][i]+'_'+qpqdata['subsys'][i]
            # Instantiate
            rho = qpqdatacoord['R_PHYS'][j]*u.kpc
            cgabs = CGMAbsSys(gal, igm_sys, name=iname, rho = rho, **kwargs)


            # add metal lines
            #if False:
            #    for jj in range(100):
            #        flgew = 0
            #        if qpqdatacoord['METAL_EW'][j][jj] > 0:
            #            flgew = 1
            #            if qpqdatacoord['METAL_EW'][j][jj] > 3*qpqdatacoord['METAL_SIGEW'][j][jj]:
            #                flgew = 3
            #        if flgew > 0:
            #            wave0 = qpqdatacoord[j]['METAL_WREST'][jj]
            #            iline = AbsLine(wave0*u.AA, closest=True, z=igm_sys.zabs)
            #            iline.attrib['EW'] = qpqdatacoord['METAL_EW'][j][jj]* u.AA
            #            iline.attrib['sig_EW'] = qpqdatacoord['METAL_SIGEW'][j][jj] * u.AA
            #            iline.attrib['flag_EW'] = flgew
            #            iline.attrib['coord'] = igm_sys.coord
            #            abslines = []
            #            abslines.append(iline)
            #            comp = AbsComponent.from_abslines(abslines, chk_vel=False)
            #            cgabs.igm_sys.add_component(comp)

            ### add metal lines
            ### not included CII*, SiII*
            lines = [['CII 1334'], ## ['CII* 1335'],
                     ['CIV 1548','CIV 1550'],['NI 1134', 'NI 1199'],['NII 1083'],['NV 1238','NV 1242'],
                     ['OI 1302'],['OVI 1037'],['MgI 2852'],['MgII 2796','MgII 2803'],['AlII 1670'],['AlIII 1854','AlIII 1862'],
                     ['SiII 1190','SiII 1193','SiII 1304','SiII 1260','SiII 1526','SiII 1808'], ## ['SiII* 1264'],
                     ['SiIII 1206'],
                     ['SiIV 1393','SiIV 1402'],['FeII 1608','FeII 2344','FeII 2374','FeII 2382','FeII 2586','FeII 2600'],
                     ['FeIII 1122']]

            for icmp in range(len(lines)):
                abslines = []
                for ii in range(len(lines[icmp])):
                    wave0 = float(lines[icmp][ii].split(' ')[1])
                    ewstr = str(lines[icmp][ii].split(' ')[1]) + 'EW'
                    ewerrstr = str(lines[icmp][ii].split(' ')[1]) + 'EWerr'
                    if ewstr == '1808EW':
                        ewstr = '1808E'
                    if ewerrstr == '1122EWerr':
                        ewerrstr = '122EWerr'
                    if qpqdata[ewstr][i] != '/':
                        iline = AbsLine(wave0 * u.AA, closest=True, z=igm_sys.zabs)
                        ## EW
                        iline.attrib['EW'] = float(qpqdata[ewstr][i]) * u.AA  # Rest EW
                        iline.attrib['sig_EW'] = float(qpqdata[ewerrstr][i]) * u.AA
                        flgew = 1
                        if iline.attrib['EW'] < 3.*iline.attrib['sig_EW']:
                            flgew=3
                        iline.attrib['flag_EW'] = flgew
                        ## column densities
                        colstr = str(lines[icmp][ii].split(' ')[0]) + 'col'
                        colerrstr = str(lines[icmp][ii].split(' ')[0]) + 'colerr'
                        iline.attrib['logN'] = qpqdata[colstr][i]
                        iline.attrib['sig_logN'] = qpqdata[colerrstr][i]
                        ##
                        iline.attrib['coord'] = igm_sys.coord
                        abslines.append(iline)
                if len(abslines) > 0:
                    comp = AbsComponent.from_abslines(abslines, chk_vel=False)
                    cgabs.igm_sys.add_component(comp)


            # add ang_sep
            qsocoord = SkyCoord(ra=qpqdatacoord['RAD'][j], dec=qpqdatacoord['DECD'][j], unit='deg')
            bgcoord = SkyCoord(ra=qpqdatacoord['RAD_BG'][j], dec=qpqdatacoord['DECD_BG'][j], unit='deg')
            cgabs.ang_sep = qsocoord.separation(bgcoord).to('arcsec')

            cgabs.vlim = [qpqdata['v_lobound'][i],qpqdata['v_upbound'][i]]*u.km/u.s

            self.cgm_abs.append(cgabs)
