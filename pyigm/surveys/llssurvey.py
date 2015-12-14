""" Class for LLS Surveys
"""
import numpy as np
import os, imp, glob
import pdb
import urllib2

from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import SkyCoord

from pyigm.surveys.igmsurvey import IGMSurvey

pyigm_path = imp.find_module('pyigm')[1]

class LLSSurvey(IGMSurvey):
    """
    An LLS Survey class
    """

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
        if len(glob.glob(summ_fil)) == 0:
            url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_DR1.fits'
            print('HD-LLS: Grabbing summary file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            with open(summ_fil, "wb") as code:
                code.write(f.read())
            print('HD-LLS: Written to {:s}'.format(summ_fil))
        else:
            print('HD-LLS: Loading summary file {:s}'.format(summ_fil))

        # Ions
        ions_fil = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_ions.json"
        if len(glob.glob(ions_fil)) == 0:
            url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_ions.json'
            print('HD-LLS: Grabbing JSON ion file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            with open(ions_fil, "wb") as code:
                code.write(f.read())
            print('HD-LLS: Written to {:s}'.format(ions_fil))
        else:
            print('HD-LLS: Loading ions file {:s}'.format(summ_fil))

        # Read
        lls_survey = cls.from_sfits(summ_fil)
        # Load ions
        lls_survey.fill_ions(jfile=ions_fil)

        # Set data path (may be None)
        spath = pyigm_path+"/data/LLS/HD-LLS/Spectra/"
        for lls in lls_survey._abs_sys:
            lls.spec_path = spath

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

        Quick and dirty port of XIDL codes..

        Parameters
        ----------
        sample : str, optional
          LLS sample
            stat : Statistical sample
            all : All LLS
            nonstat : Non-statistical sample


        Returns
        -------
        lls : Table
          Table of SDSS LLS

        """
        # LLS File
        lls_fil = pyigm_path+'/data/LLS/SDSS/lls_dr7_stat_LLS.fits.gz'
        if len(glob.glob(lls_fil)) == 0:
            url = 'https://dl.dropboxusercontent.com/u/6285549/LLS/SDSS/lls_dr7_stat_LLS.fits'
            print('SDSS-DR7: Grabbing LLS file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            with open(lls_fil, "wb") as code:
                code.write(f.read())
            print('SDSS-DR7: Written to {:s}'.format(lls_fil))
        else:
            print('SDSS-DR7: Loading LLS file {:s}'.format(lls_fil))
        lls = QTable.read(lls_fil)

        # Rename some columns?
        lls.rename_column('QSO_RA', 'RA')
        lls.rename_column('QSO_DEC', 'DEC')

        # Read
        lls_survey = cls.from_sfits(lls)
        lls_survey.ref = 'SDSS-DR7'

        # All?
        if sample == 'all':
            return lls_survey

        # QSOs file
        qsos_fil = pyigm_path+'/data/LLS/SDSS/lls_dr7_qsos_sn2050.fits.gz'
        if len(glob.glob(qsos_fil)) == 0:
            url = 'https://dl.dropboxusercontent.com/u/6285549/LLS/SDSS/lls_dr7_qsos_sn2050.fits'
            print('SDSS-DR7: Grabbing QSOs file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            with open(qsos_fil, "wb") as code:
                code.write(f.read())
            print('SDSS-DR7: Written to {:s}'.format(qsos_fil))
        else:
            print('SDSS-DR7: Loading QSOs file {:s}'.format(qsos_fil))
        qsos = QTable.read(qsos_fil)

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
