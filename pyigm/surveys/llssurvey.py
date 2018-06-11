""" Class for LLS Surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import imp, glob
import pdb
import h5py
import json
try:
    from urllib2 import urlopen # Python 2.7
except ImportError:
    from urllib.request import urlopen


from astropy.table import Column, Table
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky

from linetools import utils as ltu

from pyigm.surveys.igmsurvey import IGMSurvey
from pyigm.metallicity.pdf import MetallicityPDF
from pyigm.surveys import utils as pyisu


pyigm_path = imp.find_module('pyigm')[1]

class LLSSurvey(IGMSurvey):
    """
    An LLS Survey class
    """

    @classmethod
    def load_HST_ACS(cls, tau_LL=2, sample='stat'):
        """ Load the LLS survey using HST/ACS by O'Meara et al. 2013, ApJ, 765, 137

        Parameters
        ----------
        tau_LL : int, optional
          Sets sample
        sample : str, optional

        Returns
        -------
        lls_survey : IGMSurvey
        """
        # LLS File
        lls_fil = pyigm_path+'/data/LLS/HST/lls_acs_stat_LLS.fits.gz'
        lls = Table.read(lls_fil)

        # Rename some columns?
        lls.rename_column('QSO_RA', 'RA')
        lls.rename_column('QSO_DEC', 'DEC')

        # Generate coords
        scoords = [lls['RA'][ii]+' '+lls['DEC'][ii] for ii in range(len(lls))]
        coords = SkyCoord(scoords, unit=(u.hourangle, u.deg))

        # Read
        lls_survey = cls.from_sfits(lls, coords=coords)
        lls_survey.ref = 'HST-ACS'

        # QSOs file
        qsos_fil = pyigm_path+'/data/LLS/HST/lls_acs_qsos_sn1020.fits.gz'
        qsos = Table.read(qsos_fil)
        lls_survey.sightlines = qsos

        # Z_START, Z_END
        #   Using tau=2
        if tau_LL == 2:
            lls_survey.sightlines['Z_START'] = np.maximum(qsos['ZT2'], qsos['ZLLS'][:,0])
        else:
            pdb.set_trace()
        # zend
        zend = ltu.z_from_dv(-3000*u.km/u.s*np.ones(len(qsos)),
                             lls_survey.sightlines['ZEM'])
        lls_survey.sightlines['Z_END'] = zend

        # Stat me
        mask = lls_stat(lls_survey)
        if sample == 'stat':
            lls_survey.mask = mask
        else:
            lls_survey.mask = ~mask

        # Return
        print('HST-ACS: Loaded')
        return lls_survey

    @classmethod
    def load_HST_WFC3(cls, tau_LL=2, sample='stat'):
        """ Load the LLS survey using HST/WFC3

        by O'Meara et al. 2013, ApJ, 765, 137

        Parameters
        ----------
        tau_LL : int, optional
          Sets sample
        sample : str, optional

        Returns
        -------
        lls_survey : IGMSurvey
        """
        # LLS File
        lls_fil = pyigm_path+'/data/LLS/HST/lls_wfc3_stat_LLS.fits.gz'
        lls = Table.read(lls_fil)

        # Rename some columns?
        lls.rename_column('QSO_RA', 'RA')
        lls.rename_column('QSO_DEC', 'DEC')

        # Generate coords
        scoords = [lls['RA'][ii]+' '+lls['DEC'][ii] for ii in range(len(lls))]
        coords = SkyCoord(scoords, unit=(u.hourangle, u.deg))

        # Read
        lls_survey = cls.from_sfits(lls, coords=coords)
        lls_survey.ref = 'HST-WFC3'

        # QSOs file
        qsos_fil = pyigm_path+'/data/LLS/HST/lls_wfc3_qsos_sn1020.fits.gz'
        qsos = Table.read(qsos_fil)
        lls_survey.sightlines = qsos

        # Z_START, Z_END
        #   Using tau=2
        if tau_LL == 2:
            lls_survey.sightlines['Z_START'] = np.maximum(qsos['ZT2'], qsos['ZLLS'][:,0])
        else:
            pdb.set_trace()
        # zend
        zend = ltu.z_from_dv(-3000*u.km/u.s*np.ones(len(qsos)),
                             lls_survey.sightlines['ZEM'])
        lls_survey.sightlines['Z_END'] = zend

        # Stat me
        mask = lls_stat(lls_survey)
        if sample == 'stat':
            lls_survey.mask = mask
        else:
            lls_survey.mask = ~mask

        # Return
        print('HST-WFC3: Loaded')
        return lls_survey

    @classmethod
    def load_lowz(cls, grab_spectra=False, isys_path=None):
        """ LLS from Wotta+16 (includes Lehner+13)
        Updated to include systems excluded by Wotta+16 (high NHI)

        Parameters
        ----------
        grab_spectra : bool, optional
          Not implemented

        Return
        ------
        lls_survey
        """

        # System files
        l13_files = pyigm_path+'/data/LLS/Literature/lehner13.tar.gz'
        w16_files = pyigm_path+'/data/LLS/Literature/wotta16.tar.gz'


        # Load systems via the sys tarball.  Includes transitions
        L13_survey = pyisu.load_sys_files(l13_files, 'LLS', ref='Lehner+13')
        W16_survey = pyisu.load_sys_files(w16_files, 'LLS', ref='Wotta+16')
        lowz_LLS = L13_survey+W16_survey

        # Spectra?
        if grab_spectra:
            raise NotImplementedError("NOPE")
            specfils = glob.glob(spath+'HD-LLS_J*.fits')
            if len(specfils) < 100:
                import tarfile
                print('HD-LLS: Downloading a 155Mb file.  Be patient..')
                url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_spectra.tar.gz'
                spectra_fil = pyigm_path+'/data/LLS/HD-LLS/HD-LLS_spectra.tar.gz'
                f = urlopen(url)
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

        return lowz_LLS

    @classmethod
    def load_HDLLS(cls, load_sys=True, grab_spectra=False, isys_path=None):
        """ Default sample of LLS (HD-LLS, DR1)

        Parameters
        ----------
        grab_spectra : bool, optional
          Grab 1D spectra?  (155Mb)
        load_sys : bool, optional
          Load systems using the sys tarball
        isys_path : str, optional
          Read system files from this path

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

        # System files
        sys_files = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_sys.tar.gz"

        # Transitions
        #clm_files = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_clms.tar.gz"

        # Metallicity
        ZH_fil = pyigm_path+"/data/LLS/HD-LLS/HD-LLS_DR1_dustnhi.hdf5"
        print('HD-LLS: Loading metallicity file {:s}'.format(ZH_fil))

        # Load systems via the sys tarball.  Includes transitions
        if load_sys:  # This approach takes ~120s
            if isys_path is not None:
                lls_survey = pyisu.load_sys_files(isys_path, 'LLS',sys_path=True, ref='HD-LLS')
            else:
                lls_survey = pyisu.load_sys_files(sys_files, 'LLS', ref='HD-LLS', use_coord=True)
            lls_survey.build_all_abs_sys()
            lls_survey.fill_ions(use_components=True)
        else:
            # Read
            lls_survey = cls.from_sfits(summ_fil)
            # Load ions
            lls_survey.build_all_abs_sys()
            lls_survey.fill_ions(jfile=ions_fil)
        lls_survey.ref = 'HD-LLS'

        """
        # Load transitions
        if not skip_trans:
            print('HD-LLS: Loading transitions from {:s}'.format(clm_files))
            tar = tarfile.open(clm_files)
            for member in tar.getmembers():
                if '.' not in member.name:
                    print('Skipping a likely folder: {:s}'.format(member.name))
                    continue
                # Extract
                f = tar.extractfile(member)
                tdict = json.load(f)
                # Find system
                i0 = member.name.rfind('/')
                i1 = member.name.rfind('_clm')
                try:
                    idx = names.index(member.name[i0+1:i1])
                except ValueError:
                    print('Skipping {:s}, not statistical in DR1'.format(member.name[i0+1:i1]))
                    continue
                # Fill up
                lls_survey._abs_sys[idx].load_components(tdict)
                lls_survey._abs_sys[idx]._components = lls_survey._abs_sys[idx].subsys['A']._components
        """


        # Load metallicity
        fh5=h5py.File(ZH_fil, 'r')
        ras = []
        decs = []
        zval = []
        mkeys = list(fh5['met'].keys())  # Python 3
        mkeys.remove('left_edge_bins')
        for key in mkeys:
            radec, z = key.split('z')
            coord = ltu.radec_to_coord(radec)
            # Save
            zval.append(float(z))
            ras.append(coord.ra.value)
            decs.append(coord.dec.value)
        mcoords = SkyCoord(ras, decs, unit='deg')
        zval = np.array(zval)

        # Set data path and metallicity
        spath = pyigm_path+"/data/LLS/HD-LLS/Spectra/"
        for kk in range(lls_survey.nsys):
            lls = lls_survey.abs_sys(kk)
            lls.spec_path = spath
            # Match
            sep = lls.coord.separation(mcoords)
            mt = np.where((sep < 15*u.arcsec) & (np.abs(zval-lls.zabs) < 2e-3))[0]
            if len(mt) == 0:
                pdb.set_trace()
                raise ValueError("Bad match")
            elif len(mt) > 1:  # Take closest
                mt = np.argmin(sep)
            else:
                mt = mt[0]
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
                f = urlopen(url)
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
    def load_ribaudo(cls, sample='stat'):
        pass

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
        lls = Table.read(lls_fil)

        # Rename some columns?
        lls.rename_column('QSO_RA', 'RA')
        lls.rename_column('QSO_DEC', 'DEC')

        # Instantiate
        scoords = [lls['RA'][ii]+' '+lls['DEC'][ii] for ii in range(len(lls))]
        coords = SkyCoord(scoords, unit=(u.hourangle, u.deg))
        lls_survey = cls.from_sfits(lls, coords=coords)
        lls_survey.ref = 'SDSS-DR7'

        # QSOs file
        qsos_fil = pyigm_path+'/data/LLS/SDSS/lls_dr7_qsos_sn2050.fits.gz'
        print('SDSS-DR7: Loading QSOs file {:s}'.format(qsos_fil))
        qsos = Table.read(qsos_fil)
        lls_survey.sightlines = qsos

        # All?
        if sample == 'all':
            return lls_survey


        # Stat
        # z_em cut
        zem_min = 3.6
        lowz_q = qsos['ZEM'] < zem_min
        qsos['ZT2'][lowz_q] = 99.99

        # Survey path
        #   Using tau=2
        lls_survey.sightlines['Z_START'] = np.maximum(qsos['ZT2'], qsos['ZLLS'])
        # Trim bad ones
        bad_s = np.any([lls_survey.sightlines['Z_START'] <= 0.,
                        lls_survey.sightlines['FLG_QSO'] != 0],axis=0)
        lls_survey.sightlines['Z_START'][bad_s] = 99.99
        # zend
        zend = ltu.z_from_dv(-3000*u.km/u.s*np.ones(len(qsos)),
                             lls_survey.sightlines['ZEM'])
        lls_survey.sightlines['Z_END'] = zend

        # Generate mask
        print('SDSS-DR7: Performing stats')
        mask = lls_stat(lls_survey)
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
        tab.rename_column('zlls', 'ZLLS')
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
        gdlls = tab['ZLLS'] >= tab['Z_START']
        lls_tab = Table(tab[gdlls])
        nlls = np.sum(gdlls)
        # Set NHI to 17.5 (tau>=2)
        lls_tab.add_column(Column([17.5]*nlls, name='NHI'))
        lls_tab.add_column(Column([99.9]*nlls, name='SIGNHI'))
        lls_tab.rename_column('ZLLS', 'Z_LLS')

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


def lls_stat(LLSs, maxdz=99.99, zem_min=0., NHI_cut=17.5,
             partial=False, prox=False):
    """ Identify the statistical LLS in a survey

    Parameters
    ----------
    LLSs
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
    qsos = LLSs.sightlines

    # Search redshift
    # Turn of LLS mask
    LLSs.mask = None

    # LLS
    msk_smpl = LLSs.zem != LLSs.zem

    # Make some lists
    lls_coord = LLSs.coords
    qsos_coord = SkyCoord(ra=qsos['RA'], dec=qsos['DEC'], unit='deg')

    # Match em
    qidx, d2d, _ = match_coordinates_sky(lls_coord, qsos_coord, nthneighbor=1)
    try:
        assert np.min(d2d) < 3.6*u.arcsec  # Not sure why this value..
    except:
        pdb.set_trace()

    for qq, zabs, zem, NHI in zip(range(LLSs.nsys), LLSs.zabs, LLSs.zem, LLSs.NHI):
        idx = qidx[qq]

        '''
        small_sep = coord.separation(lls_coord) < 3.6*u.arcsec
        close_zem = np.abs(zem-lls_zem) < 0.03
        close_zabs = np.abs(zabs-lls_zabs) < dz_toler
        if np.sum(small_sep & close_zem & close_zabs) != 1:
            raise ValueError("LLS are probably too close in z")
        '''

        # Cut on NHI
        if partial & (NHI > NHI_cut):
            continue
        if (~partial) & (NHI <= NHI_cut):
            continue

        # Match to QSO RA, DEC
        assert np.abs(qsos['ZEM'][idx]-zem) < 0.03

        # Query redshift
        if (zabs > max(qsos['Z_START'][idx], qsos['Z_END'][idx] - maxdz) - 1e-4) & (
                qsos['ZEM'][idx] > zem_min):
            if (~prox) & (zabs < qsos['Z_END'][idx]):  #  Intervening
                msk_smpl[qq] = True
            if prox & (zabs >= qsos['Z_END'][idx]):  # Proximate
                msk_smpl[qq] = True

    # Return
    return msk_smpl

