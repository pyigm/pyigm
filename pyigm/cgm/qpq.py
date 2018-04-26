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
from linetools.lists.linelist import LineList

from pyigm.cgm.cgmsurvey import CGMAbsSurvey
import pyigm
from pyigm.field.galaxy import Galaxy
from .cgm import CGMAbsSys
from pyigm.abssys.igmsys import IGMSystem

try:
    basestring
except NameError:  # For Python 3
    basestring = str

ism = LineList('ISM')

def load_qpq(v):
    """ Reads table with data

    Parameters
    ----------
    v : int
      dataset (5,6,7,8)

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

    if v == 5:
        q7file = resource_filename('pyigm', 'data/CGM/QPQ/qpq7_pairs.fits.gz')
        qpqdata7 = Table.read(q7file)
        ii = np.where(qpqdata7['R_PHYS'] < 300.)
        qpqdata = qpqdata7[ii[0]]

    if v not in [5,6,7,8,8.5]:
        print('Please choose 5, 6, 7, 8, or 8.5 for QPQ6, QPQ7, QPQ8, additional information about QPQ8')
        return

    return qpqdata



class QPQ6(CGMAbsSurvey):
    """Inherits CGM Abs Survey

        Parameters:
        -----------
        cdir : str, optional
          Path to the data
        """

    def __init__(self, cdir=None, nmax=None, load=True, from_dict=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'QPQ6'
        self.ref = 'Prochaska+13'
        #
        if cdir is None:
            self.cdir = pyigm.__path__[0] + '/data/CGM/QPQ/'
        else:
            self.cdir = cdir
        if from_dict:
            self.data_file = self.cdir + 'jsons/qpq6.json'
        else:
            self.data_file = self.cdir + 'qpq6_final_cut_2015nov16.fits'

        self.nmax = nmax
        self.from_dict = from_dict

        # Init?
        if load:
            self.load_data(**kwargs)


    def load_data(self, **kwargs):
        #
        q6file = self.data_file
        if self.from_dict:
            qpq6dict = CGMAbsSurvey.from_json(q6file)
            ism = LineList('ISM')
            qpq6dict.build_systems_from_dict(llist=ism)
            self.survey_data = qpq6dict
            #self.cgm_abs = qpq6dict.cgm_abs

        else:

            qpqdata = Table.read(q6file)
            if self.nmax is not None:
                nmax = self.nmax
            else:
                nmax = len(qpqdata)
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

    def __init__(self, cdir=None, nmax=None, load=True, from_dict=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'QPQ7'
        self.ref = 'Prochaska+14'
        #
        if cdir is None:
            self.cdir = pyigm.__path__[0] + '/data/CGM/QPQ/'
        else:
            self.cdir = cdir
        #self.data_file = pyigm.__path__[0] + '/data/CGM/QPQ/'+'qpq7_pairs.fits.gz'

        if from_dict:
            self.data_file = self.cdir + 'jsons/qpq7.json'
        else:
            self.data_file = self.cdir + 'qpq7_pairs.fits.gz'

        self.nmax = nmax
        self.from_dict = from_dict


        # Init?
        if load:
            self.load_data(**kwargs)


    def load_data(self, **kwargs):
        #
        #q7file = resource_filename('pyigm', 'data/CGM/QPQ/qpq7_pairs.fits.gz')
        q7file = self.data_file

        if self.from_dict:
            qpq7dict = CGMAbsSurvey.from_json(q7file)
            ism = LineList('ISM')
            qpq7dict.build_systems_from_dict(llist=ism)
            self.survey_data = qpq7dict
            #self.cgm_abs = qpq7dict.cgm_abs

        else:
            qpqdata = Table.read(q7file)
            if self.nmax is not None:
                nmax = self.nmax
            else:
                nmax = len(qpqdata)
            #nmax = len(qpqdata)   # max number of QSOs
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



class QPQ5(CGMAbsSurvey):
    """Inherits CGM Abs Survey

        Parameters:
        -----------
        cdir : str, optional
          Path to the data
        """

    def __init__(self, load=True, from_dict=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'QPQ7'
        self.ref = 'Prochaska+13'
        #
        self.cdir = pyigm.__path__[0] + '/data/CGM/QPQ/'

        if from_dict:
            self.data_file = self.cdir + 'jsons/qpq5.json'
        else:
            self.data_file = self.cdir + 'qpq7_pairs.fits.gz'

        self.from_dict = from_dict


        # Init?
        if load:
            self.load_data(**kwargs)


    def load_data(self, **kwargs):
        #

        if self.from_dict:
            q5file = self.data_file
            qpq5dict = CGMAbsSurvey.from_json(q5file)
            ism = LineList('ISM')
            qpq5dict.build_systems_from_dict(llist=ism)
            self.survey_data = qpq5dict
            #self.cgm_abs = qpq5dict.cgm_abs

        else:
            qpqdata = load_qpq(5)
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
                aline = AbsLine(1215.67*u.AA, closest=True, z=igm_sys.zabs, linelist=ism)
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
                        iline = AbsLine(wave0*u.AA, closest=True, z=igm_sys.zabs, linelist=ism)
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

    def __init__(self, cdir=None, nmax=None, load=True, from_dict=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'QPQ8'
        self.ref = 'Lau+16'
        #
        if cdir is None:
            self.cdir = pyigm.__path__[0] + '/data/CGM/QPQ/'
        else:
            self.cdir = cdir
        #self.data_file = pyigm.__path__[0] + '/data/CGM/QPQ/'+'qpq8_pairs.fits'

        if from_dict:
            self.data_file = self.cdir + 'jsons/qpq8.json'
        else:
            self.data_file = self.cdir + 'qpq8_all_measured.dat'

        self.nmax = nmax
        self.from_dict = from_dict


        # Init?
        if load:
            self.load_data(**kwargs)

    def load_data(self, **kwargs):
        #
        #q8file = resource_filename('pyigm', 'data/CGM/QPQ/qpq8_all_measured.dat')
        q8file = self.data_file

        if self.from_dict:
            q8file = self.data_file
            qpq8dict = CGMAbsSurvey.from_json(q8file)
            ism = LineList('ISM')
            qpq8dict.build_systems_from_dict(llist=ism)
            self.survey_data = qpq8dict
            #self.cgm_abs = qpq8dict.cgm_abs


        else:

            qpqdata = Table.read(q8file,format='ascii')
            #nmax = len(qpqdata)   # max number of QSOs
            q8filecoord = resource_filename('pyigm', 'data/CGM/QPQ/qpq8_pairs.fits')
            qpqdatacoord = Table.read(q8filecoord)
            if self.nmax is not None:
                nmax = self.nmax
            else:
                nmax = len(qpqdatacoord) #(qpqdata)

            # match names with qpqdatacoord
            qnames = []
            for i in range(len(qpqdatacoord)):
                qname = qpqdatacoord['QSO'][i].strip()
                qnames.append(qname[-10:])
            qnames2 = []
            for i in range(len(qpqdata)):
                qname = qpqdata['Pair'][i]
                qnames2.append(qname)


            for j in range(nmax):   # i,j
                # match names with qpqdatacoord
                i = np.where(np.asarray(qnames2) == qnames[j])[0]
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

                iname = qpqdata['Pair'][i][0] #+'_'+qpqdata['subsys'][i]
                # Instantiate
                rho = qpqdatacoord['R_PHYS'][j]*u.kpc
                cgabs = CGMAbsSys(gal, igm_sys, name=iname, rho = rho, **kwargs)


                ### add metal lines
                ### not included CII*, SiII*
                lines = [['CII 1334'], ## ['CII* 1335'],
                         ['CIV 1548','CIV 1550'],['NI 1134', 'NI 1199'],['NII 1083'],['NV 1238','NV 1242'],
                         ['OI 1302'],['OVI 1037'],['MgI 2852'],['MgII 2796','MgII 2803'],['AlII 1670'],['AlIII 1854','AlIII 1862'],
                         ['SiII 1190','SiII 1193','SiII 1304','SiII 1260','SiII 1526','SiII 1808'], ## ['SiII* 1264'],
                         ['SiIII 1206'],
                         ['SiIV 1393','SiIV 1402'],['FeII 1608','FeII 2344','FeII 2374','FeII 2382','FeII 2586','FeII 2600'],
                         ['FeIII 1122']]

                for kk in i:

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
                            if qpqdata[ewstr][kk] != '/':
                                # find z
                                v0 = 0.5*(qpqdata['v_lobound'][kk]+ qpqdata['v_upbound'][kk]) * u.km / u.s
                                dv = v0
                                zref = igm_sys.zabs
                                z_cmp = ltu.z_from_dv(dv, zref)

                                ## vlim
                                v1 = qpqdata['v_lobound'][kk] * u.km / u.s
                                z1 = ltu.z_from_dv(v1, zref)
                                v1_cmp = ltu.dv_from_z(z1,z_cmp)

                                v2 = qpqdata['v_upbound'][kk] * u.km / u.s
                                z2 = ltu.z_from_dv(v2, zref)
                                v2_cmp = ltu.dv_from_z(z2, z_cmp)

                                # iline
                                iline = AbsLine(wave0 * u.AA, closest=True, z=z_cmp)
                                iline.attrib['coord'] = igm_sys.coord

                                ## EW
                                iline.attrib['EW'] = float(qpqdata[ewstr][kk]) * u.AA  # Rest EW
                                iline.attrib['sig_EW'] = float(qpqdata[ewerrstr][kk]) * u.AA
                                flgew = 1
                                if iline.attrib['EW'] < 3.*iline.attrib['sig_EW']:
                                    flgew=3
                                iline.attrib['flag_EW'] = flgew

                                ## column densities
                                colstr = str(lines[icmp][ii].split(' ')[0]) + 'col'
                                colerrstr = str(lines[icmp][ii].split(' ')[0]) + 'colerr'
                                iline.attrib['logN'] = qpqdata[colstr][kk]
                                iline.attrib['sig_logN'] = qpqdata[colerrstr][kk]

                                abslines.append(iline)

                        if len(abslines) > 0:
                            comp = AbsComponent.from_abslines(abslines, chk_vel=False)
                            comp.limits._vlim = [v1_cmp.value, v2_cmp.value] * u.km / u.s
                            cgabs.igm_sys.add_component(comp)


                # add ang_sep
                qsocoord = SkyCoord(ra=qpqdatacoord['RAD'][j], dec=qpqdatacoord['DECD'][j], unit='deg')
                bgcoord = SkyCoord(ra=qpqdatacoord['RAD_BG'][j], dec=qpqdatacoord['DECD_BG'][j], unit='deg')
                cgabs.ang_sep = qsocoord.separation(bgcoord).to('arcsec')

                self.cgm_abs.append(cgabs)

