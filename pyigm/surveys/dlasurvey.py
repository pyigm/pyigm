""" Class for DLA Surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import imp, glob
import pdb
import warnings

try:
    from urllib2 import urlopen # Python 2.7
except ImportError:
    from urllib.request import urlopen

from pkg_resources import resource_filename

from astropy.table import QTable, Column, Table, vstack
from astropy import units as u
from astropy.stats import poisson_conf_interval as aspci
from astropy import constants as const
from astropy.cosmology import core as acc

from linetools import utils as ltu

from pyigm.surveys.igmsurvey import IGMSurvey
from pyigm.surveys import utils as pyisu
from pyigm import utils as pyigmu

pyigm_path = imp.find_module('pyigm')[1]

# Class for DLA Survey
class DLASurvey(IGMSurvey):
    """An DLA Survey class

    Attributes:

    """
    @classmethod
    def load_HST16(cls, sample='stat'):
        """ HST Survey by Neeleman+16

        Parameters
        ----------
        sample : str, optional

        Returns
        -------
        dla_survey

        """
        # Read DLAs
        dat_file = resource_filename('pyigm', '/data/DLA/HST/HSTDLA.dat')
        dlas = Table.read(dat_file, format='ascii')

        # Read Quasars
        #qsos = Table.read(pyigm_path + '/all_qso_table.txt', format='ascii')

        # Read Sightlines
        srvy_file = resource_filename('pyigm', '/data/DLA/HST/hstpath.dat')
        survey = Table.read(srvy_file, format='ascii')

        # Add info to DLA table
        ras, decs, zems = [], [], []
        for dla in dlas:
            mt = np.where(survey['QSO'] == dla['NAME'])[0]
            if len(mt) == 0:
                pdb.set_trace()
                raise ValueError("Uh oh")
            else:
                mt = mt[0]
            # Generate RA/DEC
            row = survey[mt]
            coord = ltu.radec_to_coord('J{:02d}{:02d}{:f}{:s}{:02d}{:02d}{:f}'.format(
                row['RAh'], row['RAm'], row['RAs'], row['DE-'], row['DEd'], row['DEm'],
                row['DEs']))
            ras.append(coord.ra.value)
            decs.append(coord.dec.value)
            # zem
            zems.append(row['ZEM'])
        dlas['RA'] = ras
        dlas['DEC'] = decs
        dlas['QSO_ZEM'] = zems

        # Instantiate
        dla_survey = cls.from_sfits(dlas)
        dla_survey.ref = 'Neeleman+16'

        # Fiddle a bit
        survey.rename_column('STTMIN', 'Z_START')
        survey.rename_column('STTMAX', 'Z_END')
        stat = survey['Z_END'] > 0
        stat_survey = survey[stat]
        # Restrict to statistical sightlines
        if sample == 'stat':
            stat_survey = stat_survey[stat_survey['F_STT'] == 1]
        ras, decs, zems = [], [], []
        for row in stat_survey:
            coord = ltu.radec_to_coord('J{:02d}{:02d}{:f}{:s}{:02d}{:02d}{:f}'.format(
                    row['RAh'], row['RAm'], row['RAs'], row['DE-'], row['DEd'], row['DEm'],
                    row['DEs']))
            ras.append(coord.ra.value)
            decs.append(coord.dec.value)
        stat_survey['RA'] = ras
        stat_survey['DEC'] = decs
        stat_survey['FLG_BAL'] = 0

        # Sightlines
        dla_survey.sightlines = stat_survey

        # Stat?
        if sample in ['all', 'all_sys']:
            return dla_survey
        mask = dla_stat(dla_survey, stat_survey)
        if sample == 'stat':
            dla_survey.mask = mask & (dlas['STAT_FLG'] == 1)
        else:
            dla_survey.mask = ~mask
        # Return
        return dla_survey

    @classmethod
    def load_H100(cls, grab_spectra=False, load_sys=True, isys_path=None):
        """ Sample of unbiased HIRES DLAs compiled and analyzed by Neeleman+13

        Neeleman, M. et al. 2013, ApJ, 769, 54

        Parameters
        ----------
        grab_spectra : bool, optional
          Grab 1D spectra?  (141Mb)
        load_sys : bool, optional
          Load systems? (takes ~60s for 100 systems)
        skip_trans : bool, optional
          Skip loading transitions (takes ~60s)?
        isys_path : str, optional
          Read system files from this path

        Return
        ------
        dla_survey
        """
        # Pull from Internet (as necessary)
        summ_fil = resource_filename('pyigm', "/data/DLA/H100/H100_DLA.fits")
        print('H100: Loading summary file {:s}'.format(summ_fil))

        # Ions
        ions_fil = resource_filename('pyigm', "/data/DLA/H100/H100_DLA_ions.json")
        print('H100: Loading ions file {:s}'.format(ions_fil))

        # Transitions
        trans_fil = resource_filename('pyigm', "/data/DLA/H100/H100_DLA_clms.tar.gz")

        # System files
        sys_files = resource_filename('pyigm', "/data/DLA/H100/H100_DLA_sys.tar.gz")

        if load_sys:  # This approach takes ~90s
            print('H100: Loading systems.  This takes ~90s')
            if isys_path is not None:
                dla_survey = pyisu.load_sys_files(isys_path, 'DLA', sys_path=True, use_coord=True)
            else:
                dla_survey = pyisu.load_sys_files(sys_files, 'DLA', use_coord=True)
            dla_survey.fill_ions(use_components=True)
        else:
            # Read
            dla_survey = cls.from_sfits(summ_fil)
            # Load ions
            dla_survey.fill_ions(jfile=ions_fil)

        # Metallicities
        tbl2_file = resource_filename('pyigm', "/data/DLA/H100/H100_table2.dat")
        tbl2 = Table.read(tbl2_file, format='cds')
        # Parse for matching
        if load_sys:
            names = dla_survey.name
            qsonames = []
            zabs = []
            for name in names:
                prs = name.split('_')
                qsonames.append(prs[0])
                try:
                    zabs.append(float(prs[-1][1:]))
                except ValueError:
                    pdb.set_trace()
            qsonames = np.array(qsonames)
            zabs = np.array(zabs)
        else:
            tbl3 = Table.read(summ_fil)
            qsonames = tbl3['QSO']
            zabs = dla_survey.zabs
        # Match
        for ii, iqso, izabs in zip(range(len(tbl2)), tbl2['QSO'], tbl2['zabs']):
            mt = np.where((qsonames == iqso) & (np.abs(izabs-zabs) < 1e-3))[0]
            if len(mt) == 0:
                pdb.set_trace()
            elif len(mt) != 1:
                pdb.set_trace()
            # Metallicity
            dla_survey._abs_sys[mt[0]].ZH = tbl2['[M/H]'][ii]
            dla_survey._abs_sys[mt[0]].sig_ZH = tbl2['e_[M/H]'][ii]
            if tbl2['M'][ii] in ['S','Si','O']:
                dla_survey._abs_sys[mt[0]].flag_ZH = 1  # Alpha
            elif tbl2['M'][ii] in ['Zn']:
                dla_survey._abs_sys[mt[0]].flag_ZH = 2  # Zn
            elif tbl2['M'][ii] in ['Fe']:
                dla_survey._abs_sys[mt[0]].flag_ZH = 4  # Fe
            else:
                raise ValueError("Bad metal")
            dla_survey._abs_sys[mt[0]].elm_Z = tbl2['M'][ii]
            # Kin
            dla_survey._abs_sys[mt[0]].kin['dv'] = tbl2['dv'][ii]
            dla_survey._abs_sys[mt[0]].kin['trans'] = tbl2['trans'][ii]
            dla_survey._abs_sys[mt[0]].selection = tbl2['Select'][ii]

        dla_survey.ref = 'Neeleman+13'

        spath = pyigm_path+"/data/DLA/H100/Spectra/"
        for dla in dla_survey._abs_sys:
            dla.spec_path = spath

        # Spectra?
        if grab_spectra:
            warnings.warn("All of these spectra are in igmspec at https://github.com/specdb/specdb")
            specfils = glob.glob(spath+'H100_J*.fits')
            if len(specfils) < 100:
                import tarfile
                print('H100: Downloading a 141Mb file.  Be patient..')
                url = 'https://dl.dropboxusercontent.com/u/6285549/DLA/H100/H100_DLA_spectra.tar.gz'
                spectra_fil = pyigm_path+'/data/DLA/H100/H100_spectra.tar.gz'
                f = urlopen(url)
                with open(spectra_fil, "wb") as code:
                    code.write(f.read())
                # Unpack
                print('H100: Unpacking..')
                outdir = pyigm_path+"/data/DLA/H100"
                t = tarfile.open(spectra_fil, 'r:gz')
                t.extractall(outdir)
                # Remove
                print('H100: Removing tar file')
                os.remove(spectra_fil)
                # Done
                print('H100: All done')
            else:
                print('H100: Using files in {:s}'.format(spath))

        return dla_survey

    @classmethod
    def load_SDSS_DR5(cls, sample='stat'):
        """ Load the DLA from the SDSS-DR5 survey

        (Prochaska & Wolfe 2009, ApJ, 696, 1543)

        Parameters
        ----------
        sample : str, optional
          DLA sample
            stat : Statistical sample
            all : All DLA (NHI >= 20.3)
            all_sys : All systems identified -- Returns an LLSSurvey instead
            nonstat : Non-statistical sample


        Returns
        -------
        dla_survey : DLASurvey

        """
        from .llssurvey import LLSSurvey
        import warnings

        # LLS File
        dla_fil = pyigm_path+'/data/DLA/SDSS_DR5/dr5_alldla.fits.gz'
        print('SDSS-DR5: Loading DLA file {:s}'.format(dla_fil))
        dlas = QTable.read(dla_fil)

        # Rename some columns?
        dlas.rename_column('QSO_RA', 'RA')
        dlas.rename_column('QSO_DEC', 'DEC')

        # Cut on NHI
        if sample != 'all_sys':
            gd_dla = dlas['NHI'] >= 20.3
            dla_survey = cls.from_sfits(dlas[gd_dla])
        else:
            warnings.warn("Loading an LLSSurvey not a DLASurvey")
            dla_survey = LLSSurvey.from_sfits(dlas)

        # Read
        dla_survey.ref = 'SDSS-DR5 (PW09)'

        # g(z) file
        qsos_fil = pyigm_path+'/data/DLA/SDSS_DR5/dr5_dlagz_s2n4.fits'
        print('SDSS-DR5: Loading QSOs file {:s}'.format(qsos_fil))
        qsos = QTable.read(qsos_fil)
        qsos.rename_column('Z1', 'Z_START')
        qsos.rename_column('Z2', 'Z_END')
        # Reformat
        new_cols = []
        for key in qsos.keys():
            if key in ['GZZ', 'GZV']:
                continue
            # New one
            new_cols.append(Column(qsos[key].flatten(), name=key))
        newqsos = QTable(new_cols)
        newqsos['RA'].unit = u.deg
        newqsos['DEC'].unit = u.deg
        dla_survey.sightlines = newqsos

        # All?
        if sample in ['all', 'all_sys']:
            return dla_survey


        # Stat
        # Generate mask
        print('SDSS-DR5: Performing stats')
        mask = dla_stat(dla_survey, newqsos)
        if sample == 'stat':
            dla_survey.mask = mask
        else:
            dla_survey.mask = ~mask
        # Return
        print('SDSS-DR5: Loaded')
        return dla_survey

    @classmethod
    def load_lit(cls, dla_fil, qsos_fil, ref, sample='stat', fmt=None,
        Pdla_fil=None, **kwargs):
        """ Load the DLA from a literature sample using the files
        provided by Ruben (see Sanchez-Ramirez et al. 2016, MNRAS, 456, 4488)

        Parameters
        ----------
        dla_fil : str or Table
          Name of file containting a Table (or the Table itself) on DLAs
        qsos_fil : str or Table
          Name of file containting a Table (or the Table itself) on QSO sightlines
        fmt : str, optional
          Format for Table.read()
        sample : str, optional
          DLA sample
            stat : Statistical sample
            all : All LLS
            nonstat : Non-statistical sample
        **kwargs : optional
          Passed to dla_stat()

        Returns
        -------
        dla_survey : DLASurvey

        """
        # DLA files
        stat_dlas = Table.read(dla_fil, format=fmt)
        if Pdla_fil is not None:
            Pdlas = Table.read(Pdla_fil)
            dlas = vstack([stat_dlas,Pdlas])
        else:
            dlas = stat_dlas

        # Rename some columns?
        try:
            dlas.rename_column('Dec', 'DEC')
        except KeyError:
            pass
        else:
            dlas.rename_column('logN', 'NHI')

        # Cut on NHI
        gd_dla = dlas['NHI'] >= 20.3

        # Read
        dla_survey = cls.from_sfits(dlas[gd_dla])
        dla_survey.ref = ref

        # g(z) file
        print('Loading QSOs file {:s}'.format(qsos_fil))
        qsos = Table.read(qsos_fil, format=fmt)
        try:
            qsos.rename_column('zmin', 'Z_START')
        except KeyError:
            pass
        else:
            qsos.rename_column('zmax', 'Z_END')
            qsos.rename_column('Dec', 'DEC')
            qsos.rename_column('zem', 'ZEM')
        dla_survey.sightlines = qsos

        # BAL?
        if 'FLG_BAL' not in qsos.keys():
            qsos['FLG_BAL'] = 0

        # All?
        if sample == 'all':
            return dla_survey

        # Stat
        # Generate mask  (True = good)
        mask = dla_stat(dla_survey, qsos, **kwargs)
        if sample == 'stat':
            dla_survey.mask = mask
        else:
            dla_survey.mask = ~mask
        # Return
        print('Loaded survey')
        return dla_survey

    @classmethod
    def load_P03(cls, sample='stat'):
        """ Load the DLA from the Peroux+03 survey

        (Peroux et al. 2003, MNRAS, 346, 1103)
        IUE dataset has been removed (see Sanchez-Ramirez)
        Errors and duplicates cleaned by Sanchez-Ramirez
        Adopts a 5000km/s cutoff

        Parameters
        ----------
        sample : str, optional
          DLA sample
            stat : Statistical sample
            all : All LLS
            nonstat : Non-statistical sample

        Returns
        -------
        dla_survey : DLASurvey

        """
        # DLA files
        dla_fil = pyigm_path+'/data/DLA/P03/P03_abs.fit'
        ref = 'P03'
        qsos_fil = pyigm_path+'/data/DLA/P03/P03_zpath.fit'
        #
        dla_survey = cls.load_lit(dla_fil, qsos_fil, ref, sample=sample, skip_zem=True)
        return dla_survey


    @classmethod
    def load_G09(cls, sample='stat'):
        """ Load the DLA from the Guimaraes+09 survey

        (Guimaraes et al. 2009, A&A, 508, 133)
        Adopts a 5000km/s cutoff

        Parameters
        ----------
        sample : str, optional
          DLA sample
            stat : Statistical sample
            all : All LLS
            nonstat : Non-statistical sample

        Returns
        -------
        dla_survey : DLASurvey

        """
        # DLA files
        dla_fil = pyigm_path+'/data/DLA/G09/G09_abs.fit'
        Pdla_fil = pyigm_path+'/data/DLA/G09/G09_pabs.fit'
        ref = 'G09'
        qsos_fil = pyigm_path+'/data/DLA/G09/G09_zpath.fit'
        #
        dla_survey = cls.load_lit(dla_fil, qsos_fil, ref,
                                  Pdla_fil=Pdla_fil, sample=sample, skip_zem=True)
        return dla_survey

    @classmethod
    def load_GGG(cls, sample='stat'):
        """ Load the DLA from GGG

        (Crighton et al. 2015, MNRAS, 452, 217
        http://adsabs.harvard.edu/abs/2015MNRAS.452..217C)

        Parameters
        ----------
        sample : str, optional

        Returns
        -------
        dla_survey : DLASurvey
        """
        # DLA files
        dla_fil = pyigm_path+'/data/DLA/GGG/GGG_DLA.dat'
        ref = 'GGG'
        qsos_fil = pyigm_path+'/data/DLA/GGG/GGG_QSO.dat'
        #
        dla_survey = cls.load_lit(dla_fil, qsos_fil, ref, sample=sample, fmt='ascii')
        return dla_survey

    @classmethod
    def load_XQ100(cls, sample='stat'):
        """ Load the DLA from XQ-100

        (Sanchez-Ramirez et al. 2016, MNRAS, 456, 4488)
        http://adsabs.harvard.edu/abs/2016MNRAS.456.4488S

        Parameters
        ----------
        sample : str, optional
          DLA sample
            stat : Statistical sample
            all : All DLA
            nonstat : Non-statistical sample

        Returns
        -------
        dla_survey : DLASurvey
        """
        # DLA files
        dla_fil = pyigm_path+'/data/DLA/XQ-100/XQ100_abs.fit'
        Pdla_fil = pyigm_path+'/data/DLA/XQ-100/XQ100_pabs.fit'
        ref = 'XQ-100'
        qsos_fil = pyigm_path+'/data/DLA/XQ-100/XQ100_zpath.fit'
        #
        dla_survey = cls.load_lit(dla_fil, qsos_fil, ref,Pdla_fil=Pdla_fil,
                                  sample=sample, skip_zem=True)
        return dla_survey

    @classmethod
    def neeleman13_tree(cls):
        """ Read Neeleman+13 data from the DLA tree (deprecated)
        Returns
        -------
        dlasurvey : IGMSurvey
        """
        # Default sample of DLA:  Neeleman
        if os.getenv('DLA') is None:
            print('Need to grab the DLA tree from JXP')
            return None
        dlasurvey = cls.from_flist('Lists/Neeleman13.lst', tree=os.environ.get('DLA'))
        dlasurvey.ref = 'Neeleman+13'

        # Return
        return dlasurvey

    def __init__(self, **kwargs):
        IGMSurvey.__init__(self, 'DLA', **kwargs)

        # define the cosmology (for H0)
        try:
            _ = self.cosmo
        except ValueError:
            self.cosmo = acc.FlatLambdaCDM(70., 0.3)

        # Load fitted quantities
        self.fitted = self.load_fitted()

    def binned_loz(self, zbins, NHI_mnx=(20.3, 23.00)):
        """ Calculate l(z) empirically in zbins for an interval in NHI
        Wrapper on lox

        Parameters
        ----------
        zbins : list
          Defines redshift intervals
          e.g.  [2., 2.5, 3., 4.]
        NHI_mnx : tuple, optional
          min/max of NHI for evaluation

        Returns
        -------
        lz, sig_lz_lower, sig_lz_upper : ndarray
        """
        return self.calculate_lox(zbins, NHI_mnx=NHI_mnx, use_Dz=True)

    def binned_lox(self, zbins, NHI_mnx=(20.3, 23.00), use_Dz=False):
        """ Calculate l(X) in zbins for an interval in NHI
        Parameters
        ----------
        zbins : list
          Defines redshift intervals
          e.g.  [2., 2.5, 3., 4.]
        NHI_mnx : tuple, optional
          min/max of NHI for evaluation
        use_gz : bool, optional
          Use gz instead of gX.
          This forces the calculation of l(z) instead of l(X)

        Returns
        -------
        lX, sig_lX_lower, sig_lX_upper : ndarray
        """
        # assign the nhbins
        nhbins = np.array(NHI_mnx)

        # generate the fN components
        fncomp = self.__generate_fncomp__(nhbins, zbins)

        # get the absorption path length
        dXtot = self.__find_dXtot__(zbins, calc_Dz=use_Dz)

        # total number of absorbers + poisson uncertainty
        Ntot = fncomp.sum(axis=0)
        Nunc = aspci(Ntot, interval='frequentist-confidence')

        # l(X)
        lX = Ntot / dXtot
        lX_lo = Nunc[0, :] / dXtot
        lX_hi = Nunc[1, :] / dXtot

        return lX, lX - lX_lo, lX_hi - lX

    def binned_rhoHI(self, zbins, nhbins=(20.3, 23.), nboot=1000):
        """ Calculate the mass density in HI

        Parameters
        ----------
        zbins : list
        nhbins : list

        Returns
        -------
        rhoHI : ndarray
          Evaluation of HI mass density, with units
        rhoHI_lo : ndarray
          Error estimate (low side)
        rhoHI_hi : ndarray
          Error estimate (high side)
        """
        # generate the fN components
        fncomp = self.__generate_fncomp__(nhbins, zbins)

        # get the absorption path length
        dXtot = self.__find_dXtot__(zbins)

        # get the total column density per zbin
        NHtot = self.__find_NHtot__(zbins, NH_mnx=(np.min(nhbins), np.max(nhbins)))

        # bootstrap NH_average uncertainty
        #NHunc = self.__bootstrap_rhohi__(fncomp, nhbins, zbins, nboot=nboot)
        NHunc = 1e20

        # total number of absorbers + poisson uncertainty
        Ntot = fncomp.sum(axis=0)
        Nunc = aspci(Ntot, interval='frequentist-confidence')

        frac_unc = np.sqrt(np.power(abs(Nunc - Ntot) / Ntot, 2) +
                           np.power(np.array([NHunc / (NHtot / Ntot), ] * 2), 2))
        # rho_HI
        rhoHI = NHtot / dXtot
        rhoHI_lo = rhoHI * frac_unc[0, :]
        rhoHI_hi = rhoHI * frac_unc[1, :]

        # Constants
        rhoHI = rhoHI * (const.m_p.cgs * self.cosmo.H0 /
                         const.c.cgs / (u.cm ** 2)).to(u.Msun / u.Mpc ** 3)
        rhoHI_lo = rhoHI_lo * (const.m_p.cgs * self.cosmo.H0 /
                               const.c.cgs / (u.cm ** 2)).to(u.Msun / u.Mpc ** 3)
        rhoHI_hi = rhoHI_hi * (const.m_p.cgs * self.cosmo.H0 /
                               const.c.cgs / (u.cm ** 2)).to(u.Msun / u.Mpc ** 3)

        return rhoHI, rhoHI_lo, rhoHI_hi

    def binned_fn(self, nhbins, zbins, log=False):
        """ Calculate f(N,X) empirically in bins of NHI and z

        Parameters
        ----------
        nhbins : list
        zbins : list
        log : bool, optional
          Report log10 values?

        Returns
        -------
        fn : ndarray
          log10 f(N,X)
        fn_lo : ndarray
          error in fn (low side)
        fn_hi : ndarray
          error in fn (high side)

        """

        # generate the fN components
        fncomp = self.__generate_fncomp__(nhbins, zbins)

        # calculate the uncertainty on the bins
        fnunc = aspci(fncomp, interval='frequentist-confidence')

        # get the absorption path length
        dXtot = self.__find_dXtot__(zbins)

        # find the nhi bin size
        dNHI = np.power(10, nhbins[1:]) - np.power(10, nhbins[:-1])

        # calculate the fN values
        fn = np.transpose(np.transpose(fncomp / dXtot) / dNHI)
        fn_lo = np.transpose(np.transpose(fnunc[0] / dXtot) / dNHI)
        fn_hi = np.transpose(np.transpose(fnunc[1] / dXtot) / dNHI)

        if log:
            return np.log10(fn), np.log10(fn) - np.log10(fn_lo), np.log10(fn_hi) - np.log10(fn)
        else:
            return fn, fn - fn_lo, fn_hi - fn

    def fitted_lz(self, z):
        """
        Parameters
        ----------
        z

        Returns
        -------
        loz : float or ndarray

        """

    def load_fitted(self):
        """ Load the fit info
        Returns
        -------
        dla_fits : dict

        """
        # File
        fits_file = resource_filename('pyigm','data/DLA/dla_fits.json')
        dla_fits = ltu.loadjson(fits_file)
        # Return
        return dla_fits


    def __generate_fncomp__(self, nhbins, zbins):
        """  Generate binned evaluation of f(NHI,X)

        Parameters
        ----------
        nhbins : list
          Defines NHI bins for f(NHI,X) evaluation, e.g.
            [20.3, 20.6, 21.0, 21.5, 23.]
        zbins : list

        Returns
        -------
        fncomp : ndarray
          f(NHI,X)

        """
        # calculate the total absorption path length g(X) from g(z)
        #z, gX = self.__calculate_gX__()

        # create the fn array
        zabs = self.__getattr__('zabs')
        nhi = self.__getattr__('NHI')
        fncomp = np.histogram2d(nhi, zabs, bins=[nhbins, zbins])[0]

        return fncomp

    def __find_dXtot__(self, zbins, calc_Dz=False):
        """ Calculate DX in zbins
        Parameters
        ----------
        zbins : list
        calc_Dz : bool, optional
          Return Dztot instead of DXtot

        Returns
        -------
        dXtot : ndarray
          dX for the full survey

        """
        # get z, g(z)
        z, gz = self.calculate_gz()
        dz = z[1] - z[0]
        #
        if not calc_Dz:
            dXdz = pyigmu.cosm_xz(z, cosmo=self.cosmo, flg_return=1)
        else:
            dXdz = 1.

        dXtot = np.zeros(len(zbins) - 1)
        for kk in range(len(zbins) - 1):
            # the indices of values within the redshift range
            idx = np.where((z >= zbins[kk]) & (z < zbins[kk + 1]))
            dXtot[kk] = np.sum((gz*dz*dXdz)[idx])

        return dXtot

    def __find_NHtot__(self, zbins, NH_mnx):
        """ Calculate the summed NHI

        Parameters
        ----------
        zbins : list
        NH_mnx : list or tuple
          Min/max for the sum

        Returns
        -------
        NHtot : ndarray
          Sum in the zbin intervals
        """

        # import the values from the survey
        zabs = self.__getattr__('zabs')
        nhi = self.__getattr__('NHI')

        # trim the data to only the NHI values in the range
        idx = np.where((nhi >= NH_mnx[0]) & (nhi < NH_mnx[1]))
        zabs = zabs[idx]
        nhi = nhi[idx]

        # find the total column density
        NHtot = np.zeros(len(zbins) - 1)
        for kk in range(len(zbins) - 1):
            idx2 = np.where((zabs >= zbins[kk]) & (zabs < zbins[kk + 1]))
            # total column density in the zbin
            NHtot[kk] = np.sum(np.power(10., nhi[idx2]))

        return NHtot

    '''
    def __calculate_gX__(self):
        """ Calculate g(X) from g(z)
        Returns
        -------
        z : ndarray
          Redshifts where g(z) is evaluated
        gX : ndarray
          g(X)
        """

        # calculate the total absorption path length g(X) from g(z)
        z, gz = self.calculate_gz()
        dXdz = pyigmu.cosm_xz(z, cosmo=self.cosmo, flg_return=1)
        gX = gz / dXdz #* dz  # THIS LOOKS WRONG TO ME!

        return z, gX
    '''

    def __bootstrap_rhohi__(self, fncomp, nhbins, zbins, nboot=1000):
        """ Calculate standard deviation in <NHI> with a bootstrap technique

        Parameters
        ----------
        fncomp
        nhbins
        zbins
        nboot

        Returns
        -------
        NHunc : float
          Standard deviation in <NHI> boostratp realizations

        """

        # calculate the uncertainty on fncomp
        fnunc = aspci(fncomp, interval='frequentist-confidence')

        # create a nboot, len(nhbins)-1, len(zbins)-1, array
        tfncomp = np.array([fncomp, ] * nboot)

        # populate the array with different realizations of the distribution
        randarr = np.random.randn(nboot, len(nhbins) - 1, len(zbins) - 1)
        tfncomp = (tfncomp + randarr.clip(0, 99) * (fnunc[1, :, :] - fncomp) +
                   randarr.clip(-99, 0) * (fncomp - fnunc[0, :, :]))

        # calculate the average column density
        Ntot = np.sum(tfncomp, axis=1)
        NH_avg = np.power(10., (nhbins[:-1] + 0.2 * (nhbins[1:] - nhbins[:-1])))  ##NOTE THE KLUDGE!
        NH_tavg = np.dot(NH_avg, tfncomp) / Ntot

        # trim the 1 sigma edges off and find the min/max
        # trimmed=scipy.stats.trimboth(NH_tot,scipy.special.erfc(1),axis=0)
        # NHunc=np.array([trimmed.min(axis=0),trimmed.max(axis=0)])

        # Should we do this in log space??  Does it matter??

        # just calculate the standard deviation of the distribution
        NHunc = np.std(NH_tavg, axis=0)

        return NHunc


def dla_stat(DLAs, qsos, vprox=None, buff=3000.*u.km/u.s,
             zem_min=0., flg_zsrch=0, vmin=0.*u.km/u.s,
             LLS_CUT=None, partial=False, prox=False,
             zem_tol=0.03, skip_zem=False):
    """ Identify the statistical DLA in a survey
    Note that this algorithm ignores any existing mask

    Parameters
    ----------
    DLAs : DLASurvey
    qsos : Table
      keys must include RA, DEC, ZEM, Z_START
    vmin : Quantity
    vprox
    maxdz
    zem_min
    buff : Quantity
      Buffer velocity in Proximate analysis [not ready for this]
    NHI_cut
    flg_zsrch
    dz_toler
    partial : bool, optional
      Analyze partial LLS? [pLLS]
    prox : bool, optional
      Proximate LLS? [PLLS]
    zem_tol : float, optional
      Tolerance in zem
    skip_zem : bool, optional
      Skip check on zem?? -- For proximates
      XQ-100 needs to be skipped

    Returns
    -------
    msk_smpl : bool array
      True = statistical
    """
    import warnings
    from astropy.coordinates import SkyCoord, match_coordinates_sky
    # Check for mask
    if DLAs.mask is not None:
        warnings.warn("Resetting mask to None.  Be careful here")
        DLAs.mask = None
    # DLA
    msk_smpl = DLAs.zem != DLAs.zem
    #zmax = ltu.z_from_dv(vprox, qsos['ZEM'])
    zmin = ltu.z_from_dv(vmin*np.ones(len(qsos)), qsos['Z_START'].data) # vmin must be array-like to be applied to each individual qsos['Z_START']

    # Make some lists
    qsos_coord = SkyCoord(ra=qsos['RA'], dec=qsos['DEC'], unit='deg')
    dla_coord = DLAs.coord

    idx, d2d, d3d = match_coordinates_sky(dla_coord, qsos_coord, nthneighbor=1)
    close = d2d < 1.*u.arcsec

    for qq, idla in enumerate(DLAs._abs_sys):
        # In stat?
        if close[qq]:
            if (np.abs(idla.zem-qsos['ZEM'][idx[qq]]) < zem_tol) or (skip_zem):
                if ((idla.zabs >= zmin[idx[qq]]) &
                        (idla.zabs <= qsos['Z_END'][idx[qq]]) & (qsos[idx[qq]]['FLG_BAL'] != 2)):
                        msk_smpl[qq] = True
    # Return
    return msk_smpl


def load_dla_surveys():
    """ Load up a select set of the DLA surveys
    for statistical analysis

    Returns
    -------
    surveys : list of DLASurvey objects

    """
    # Load surveys
    print("Loading DLA surveys...")
    print('Loading HST')
    hst = DLASurvey.load_HST16()
    print('Loading SDSS-DR5')
    sdss = DLASurvey.load_SDSS_DR5()
    print('Loading XQ100')
    xq100 = DLASurvey.load_XQ100()
    print('Loading GGG')
    ggg = DLASurvey.load_GGG()

    # Return
    surveys = (hst, sdss, xq100, ggg)
    return surveys


def load_dla_fits(fit_file=None):
    if fit_file is None:
        fit_file = resource_filename('pyigm', 'data/DLA/dla_fits.json')
    if os.path.isfile(fit_file):
        dla_fits = ltu.loadjson(fit_file)
    else:
        dla_fits = {}
    # Return
    return dla_fits, fit_file


def update_dla_fits(new_fits):
    import datetime
    import getpass
    # Load existing
    dla_fits, fit_file = load_dla_fits()

    # Write fit
    date = str(datetime.date.today().strftime('%Y-%b-%d'))
    user = getpass.getuser()
    #
    for key in new_fits:
        new_fits[key]['CreationDate'] = date
        new_fits[key]['User'] = user
        # Add
        if key not in dla_fits.keys():
            dla_fits[key] = {}
        for subkey in new_fits[key]:
            dla_fits[key][subkey] = new_fits[key][subkey]
    # Write
    pdb.set_trace()
    jdfits = ltu.jsonify(dla_fits)
    ltu.savejson(fit_file, jdfits, easy_to_read=True, overwrite=True)
    print("Wrote: {:s}".format(fit_file))
