""" Class for LLS Surveys
"""
import numpy as np
import imp, glob
import pdb
import urllib2
import h5py


from astropy.table import QTable, Column, Table
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools import utils as ltu

from pyigm.surveys.igmsurvey import IGMSurvey
from pyigm.metallicity.pdf import MetallicityPDF

pyigm_path = imp.find_module('pyigm')[1]

class LLSSurvey(IGMSurvey):
    """
    An LLS Survey class
    """

    @classmethod
    def load_HST_ACS(cls):
        """ Load the LLS survey using HST/ACS by O'Meara et al. 2013, ApJ, 765, 137

        Parameters
        ----------

        Returns
        -------
        lls_survey : IGMSurvey
        """
        # LLS File
        lls_fil = pyigm_path+'/data/LLS/HST/lls_acs_stat_LLS.fits.gz'
        lls = QTable.read(lls_fil)

        # Rename some columns?
        lls.rename_column('QSO_RA', 'RA')
        lls.rename_column('QSO_DEC', 'DEC')

        # Read
        lls_survey = cls.from_sfits(lls)
        lls_survey.ref = 'HST-ACS'

        # QSOs file
        qsos_fil = pyigm_path+'/data/LLS/HST/lls_acs_qsos_sn1020.fits.gz'
        qsos = QTable.read(qsos_fil)
        lls_survey.sightlines = qsos

        # Return
        print('HST-ACS: Loaded')
        return lls_survey

    @classmethod
    def load_HST_WFC3(cls):
        """ Load the LLS survey using HST/WFC3

        by O'Meara et al. 2013, ApJ, 765, 137

        Parameters
        ----------

        Returns
        -------
        lls_survey : IGMSurvey
        """
        # LLS File
        lls_fil = pyigm_path+'/data/LLS/HST/lls_wfc3_stat_LLS.fits.gz'
        lls = QTable.read(lls_fil)

        # Rename some columns?
        lls.rename_column('QSO_RA', 'RA')
        lls.rename_column('QSO_DEC', 'DEC')

        # Read
        lls_survey = cls.from_sfits(lls)
        lls_survey.ref = 'HST-WFC3'

        # QSOs file
        qsos_fil = pyigm_path+'/data/LLS/HST/lls_wfc3_qsos_sn1020.fits.gz'
        qsos = QTable.read(qsos_fil)
        lls_survey.sightlines = qsos

        # Return
        print('HST-WFC3: Loaded')
        return lls_survey

    @classmethod
    def load_HDLLS(cls, grab_spectra=False):
        """ Default sample of LLS (HD-LLS, DR1)

        Parameters
        ----------
        grab_spectra : bool, optional
          Grab 1D spectra?  (155Mb)

        Return
        ------
        lls_survey
        """
        # Pull from Internet (as necessary)
        summ_fil = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_DR1.fits"
        print('HD-LLS: Loading summary file {:s}'.format(summ_fil))

        # Ions
        ions_fil = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_ions.json"
        print('HD-LLS: Loading ions file {:s}'.format(ions_fil))

        # Transitions
        clm_fil = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_clms.json.gz"
        print('HD-LLS: Loading transitions file {:s}'.format(clm_fil))

        # Metallicity
        ZH_fil = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_DR1_dustnhi.hdf5"
        print('HD-LLS: Loading metallicity file {:s}'.format(ZH_fil))

        # Read
        lls_survey = cls.from_sfits(summ_fil)
        names = lls_survey.name
        # Load transitions
        clm_dict = ltu.loadjson(clm_fil)
        for key in clm_dict.keys():
            idx = np.where(key == names)[0]
            if len(idx) != 1:
                raise ValueError("Cannot match this LLS: {:s}".format(key))
            pdb.set_trace()
        # Load ions
        lls_survey.fill_ions(jfile=ions_fil)
        # Load metallicity
        fh5=h5py.File(ZH_fil, 'r')
        # Get coords
        ras = []
        decs = []
        zval = []
        mkeys = fh5['met'].keys()
        mkeys.remove('left_edge_bins')
        for key in mkeys:
            radec, z = key.split('z')
            coord = ltu.radec_to_coord(radec)
            # Save
            zval.append(float(z))
            ras.append(coord.ra.value)
            decs.append(coord.dec.value)
        mcoords = SkyCoord(ra=ras*u.deg, dec=decs*u.deg)

        # Set data path and metallicity
        spath = pyigm_path+"/data/LLS/HD-LLS/Spectra/"
        for lls in lls_survey._abs_sys:
            lls.spec_path = spath
            # Match
            sep = lls.coord.separation(mcoords)
            mt = np.where((sep < 15*u.arcsec) & (np.abs(zval-lls.zabs) < 2e-3))[0]
            if len(mt) == 0:
                pdb.set_trace()
                raise ValueError("Bad match")
            elif len(mt) > 1:  # Take closest
                mt = np.argmin(sep)
            # Save
            lls.metallicity = MetallicityPDF(fh5['met']['left_edge_bins']+
                                             fh5['met']['left_edge_bins'].attrs['BINSIZE']/2.,
                                             fh5['met'][mkeys[mt]])

        # Spectra?
        if grab_spectra:
            specfils = glob.glob(spath+'HD-LLS_J*.fits')
            if len(specfils) < 100:
                import tarfile
                print('HD-LLS: Downloading a 155Mb file.  Be patient..')
                url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_spectra.tar.gz'
                spectra_fil = pyigm_path+'/data/LLS/HD-LLS/HD-LLS_spectra.tar.gz'
                f = urllib2.urlopen(url)
                with open(spectra_fil, "wb") as code:
                    code.write(f.read())
                # Unpack
                print('HD-LLS: Unpacking..')
                outdir = pyigm_path+"/data/LLS/HD-LLS"
                t = tarfile.open(spectra_fil, 'r:gz')
                t.extractall(outdir)
                # Done
                print('HD-LLS: All done')
            else:
                print('HD-LLS: Using files in {:s}'.format(spath))

        return lls_survey

    @classmethod
    def load_SDSS_DR7(cls, sample='stat'):
        """ Load the LLS from the SDSS-DR7 survey (Prochaska+10, ApJ, 718, 391)

        Parameters
        ----------
        sample : str, optional
          LLS sample
            stat : Statistical sample
            all : All LLS
            nonstat : Non-statistical sample


        Returns
        -------
        lls_survey : IGMSurvey

        """
        # LLS File
        lls_fil = pyigm_path+'/data/LLS/SDSS/lls_dr7_stat_LLS.fits.gz'
        print('SDSS-DR7: Loading LLS file {:s}'.format(lls_fil))
        lls = QTable.read(lls_fil)

        # Rename some columns?
        lls.rename_column('QSO_RA', 'RA')
        lls.rename_column('QSO_DEC', 'DEC')

        # Read
        lls_survey = cls.from_sfits(lls)
        lls_survey.ref = 'SDSS-DR7'

        # QSOs file
        qsos_fil = pyigm_path+'/data/LLS/SDSS/lls_dr7_qsos_sn2050.fits.gz'
        print('SDSS-DR7: Loading QSOs file {:s}'.format(qsos_fil))
        qsos = QTable.read(qsos_fil)
        lls_survey.sightlines = qsos

        # All?
        if sample == 'all':
            return lls_survey


        # Stat
        # z_em cut
        zem_min = 3.6
        lowz_q = qsos['ZEM'] < zem_min
        qsos['ZT2'][lowz_q] = 99.99

        # Generate mask
        print('SDSS-DR7: Performing stats (~60s)')
        mask = lls_stat(lls_survey, qsos)
        if sample == 'stat':
            lls_survey.mask = mask
        else:
            lls_survey.mask = ~mask
        # Return
        print('SDSS-DR7: Loaded')
        return lls_survey

    @classmethod
    def load_mage_z3(cls, sample='all'):
        """ Load the LLS table from the z~3 MagE survey

        (Fumagalli et al. 2013, ApJ, 775, 78)

        Parameters
        ----------
        sample : str
          Survey sample
            * all -- All
            * non-color -- Restricts to quasars that were *not* color-selected
            * color -- Restricts to quasars that were color-selected

        Returns
        -------
        lls_survey : IGMSurvey
          Includes all quasars observed in the survey
          And all the LLS

        """
        # LLS File
        survey_fil = pyigm_path+'/data/LLS/HD-LLS/fumagalli13_apj775_78_tab1+2.fits'
        tab = Table.read(survey_fil)

        # Rename some columns
        tab.rename_column('RAJ2000', 'RA')
        tab['RA'].unit = u.deg
        tab.rename_column('DEJ2000', 'DEC')
        tab['DEC'].unit = u.deg
        tab.rename_column('zqso', 'Z_QSO')
        tab.rename_column('zlls', 'Z_LLS')
        tab.rename_column('zend', 'Z_START')  # F13 was opposite of POW10
        tab.rename_column('zstart', 'Z_END')  # F13 was opposite of POW10

        # Cut table
        if sample == 'all':
            pass
        elif sample == 'non-color':
            NC = np.array([True if row['n_Name'][0] == 'N' else False for row in tab])
            tab = tab[NC]
        elif sample == 'color':
            Clr = [True if row['n_Name'][0] == 'C' else False for row in tab]
            tab = tab[Clr]

        # Good LLS
        lls = tab['Z_LLS'] >= tab['Z_START']
        lls_tab = QTable(tab[lls])
        nlls = np.sum(lls)
        # Set NHI to 17.8 (tau>=2)
        lls_tab.add_column(Column([17.8]*nlls, name='NHI'))
        lls_tab.add_column(Column([99.9]*nlls, name='SIGNHI'))

        # Generate survey
        lls_survey = cls.from_sfits(lls_tab)
        lls_survey.ref = 'z3_MagE'
        lls_survey.sightlines = tab

        return lls_survey


    def __init__(self, **kwargs):
        IGMSurvey.__init__(self, 'LLS', **kwargs)

    def cut_nhi_quality(self, sig_cut=0.4):
        """ Cut the LLS on NHI quality.

        Could put this in Absline_Survey

        Parameters
        ----------
        sig_cut : float, optional
            Limit to include as quality

        Returns
        -------
        gdNHI : ndarray
        bdNHI : ndarray
            Indices for those LLS that are good/bad
            gdNHI is a numpy array of indices
            bdNHI is a boolean array
        """
        # Cut
        gdNHI = np.where( (self.sigNHI[:, 0] < sig_cut)
                        & (self.sigNHI[:, 1] < sig_cut))[0]
        # Mask
        bdNHI = (self.NHI == self.NHI)
        bdNHI[gdNHI] = False

        # Return
        return gdNHI, bdNHI


def lls_stat(LLSs, qsos, vprox=3000.*u.km/u.s, maxdz=99.99,
             zem_min=0., NHI_cut=17.5, flg_zsrch=0, dz_toler=0.04,
             LLS_CUT=None, partial=False, prox=False):
    """ Identify the statistical LLS in a survey

    Parameters
    ----------
    vprox
    maxdz
    zem_min
    NHI_cut
    flg_zsrch
    dz_toler
    partial : bool, optional
      Analyze partial LLS? [pLLS]
    prox : bool, optional
      Proximate LLS? [PLLS]

    Returns
    -------
    msk_smpl : bool array
      True = statistical
    """
    from linetools.utils import z_from_v
    # Search redshift
    if flg_zsrch == 0:
        zsrch = qsos['ZT2']
    elif flg_zsrch == 1:
        zsrch = qsos['ZT1']
    elif flg_zsrch == 2:
        zsrch = qsos['ZT0']
    # Modify by LLS along QSO sightline as required
    if LLS_CUT is not None:
        pdb.set_trace()
        #zsrch = zsrch > qsos.zlls[LLS_CUT]

    # LLS
    msk_smpl = LLSs.zem != LLSs.zem
    zmax = z_from_v(qsos['ZEM'], vprox)

    # Make some lists
    lls_coord = LLSs.coord
    lls_zem = LLSs.zem
    lls_zabs = LLSs.zabs
    qsos_coord = SkyCoord(ra=qsos['RA']*u.deg, dec=qsos['DEC']*u.deg)

    for qq, ills in enumerate(LLSs._abs_sys):
        # Two LLS on one sightline?
        small_sep = ills.coord.separation(lls_coord) < 3.6*u.arcsec
        close_zem = np.abs(ills.zem-lls_zem) < 0.03
        close_zabs = np.abs(ills.zabs-lls_zabs) < dz_toler
        if np.sum(small_sep & close_zem & close_zabs) != 1:
            raise ValueError("LLS are probably too close in z")

        # Cut on NHI
        if partial & (ills.NHI > NHI_cut):
            continue
        if ~partial & (ills.NHI <= NHI_cut):
            continue

        # Match to QSO RA, DEC
        idx = np.where( (ills.coord.separation(qsos_coord) < 3.6*u.arcsec) &
                        (np.abs(qsos['ZEM']-ills.zem) < 0.03))[0]
        if len(idx) != 1:
            raise ValueError("Problem with matches")

        # Query redshift
        if ((zsrch[idx] > 0.) &
                (ills.zabs > max(zsrch[idx], qsos['ZEM'][idx] - maxdz) - 1e-4) &
                (qsos['ZEM'][idx] > zem_min)):
            if (~prox) & (ills.zabs < zmax[idx]):  #  Intervening
                msk_smpl[qq] = True
            if prox & (ills.zabs >= zmax[idx]):  # Proximate
                msk_smpl[qq] = True

    # Return
    return msk_smpl
