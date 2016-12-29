""" Class for DLA Surveys
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import imp, glob
import pdb
import urllib2
import json

from astropy.table import QTable, Column, Table, vstack
from astropy import units as u
#from astropy.coordinates import SkyCoord
from linetools import utils as ltu
from pyigm.surveys.igmsurvey import IGMSurvey
from pyigm.surveys import utils as pyisu

pyigm_path = imp.find_module('pyigm')[1]


# Class for DLA Survey
class DLASurvey(IGMSurvey):
    """An DLA Survey class

    Attributes:

    """
    @classmethod
    def load_H100(cls, grab_spectra=False, load_sys=True, isys_path=None):  #skip_trans=True):
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
        summ_fil = pyigm_path+"/data/DLA/H100/H100_DLA.fits"
        print('H100: Loading summary file {:s}'.format(summ_fil))

        # Ions
        ions_fil = pyigm_path+"/data/DLA/H100/H100_DLA_ions.json"
        print('H100: Loading ions file {:s}'.format(ions_fil))

        # Transitions
        trans_fil = pyigm_path+"/data/DLA/H100/H100_DLA_clms.tar.gz"

        # System files
        sys_files = pyigm_path+"/data/DLA/H100/H100_DLA_sys.tar.gz"

        if load_sys:  # This approach takes ~120s
            print('H100: Loading systems.  This takes ~120s')
            if isys_path is not None:
                dla_survey = pyisu.load_sys_files(isys_path, 'DLA', sys_path=True)
            else:
                dla_survey = pyisu.load_sys_files(sys_files, 'DLA')
            dla_survey.fill_ions(use_components=True)
        else:
            # Read
            dla_survey = cls.from_sfits(summ_fil)
            # Load ions
            dla_survey.fill_ions(jfile=ions_fil)

        dla_survey.ref = 'Neeleman+13'

        """
        # Load transitions
        names = list(dla_survey.name)
        if not skip_trans:
            print('H100: Loading transitions file {:s}'.format(trans_fil))
            tar = tarfile.open(trans_fil)
            for member in tar.getmembers():
                if '.' not in member.name:
                    print('Skipping a likely folder: {:s}'.format(member.name))
                    continue
                # Extract
                f = tar.extractfile(member)
                # Need to fix for 3.4
                tdict = json.load(f)
                # Find system
                i0 = member.name.rfind('/')
                i1 = member.name.rfind('_clm')
                try:
                    idx = names.index(member.name[i0+1:i1])
                except ValueError:
                    pdb.set_trace()
                # Fill up
                dla_survey._abs_sys[idx].load_components(tdict)
        """


        spath = pyigm_path+"/data/DLA/H100/Spectra/"
        for dla in dla_survey._abs_sys:
            dla.spec_path = spath

        # Spectra?
        if grab_spectra:
            pdb.set_trace()  # USE IGMSPEC!!
            specfils = glob.glob(spath+'H100_J*.fits')
            if len(specfils) < 100:
                import tarfile
                print('H100: Downloading a 141Mb file.  Be patient..')
                url = 'https://dl.dropboxusercontent.com/u/6285549/DLA/H100/H100_DLA_spectra.tar.gz'
                spectra_fil = pyigm_path+'/data/DLA/H100/H100_spectra.tar.gz'
                f = urllib2.urlopen(url)
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
    def load_lit(cls, dla_fil, qsos_fil, ref, sample='stat', Pdla_fil=None):
        """ Load the DLA from a literature sample using the files
        provided by Ruben (see Sanchez-Ramirez et al. 2016, MNRAS, 456, 4488)

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
        stat_dlas = Table.read(dla_fil)
        if Pdla_fil is not None:
            Pdlas = Table.read(Pdla_fil)
            dlas = vstack([stat_dlas,Pdlas])
        else:
            dlas = stat_dlas

        # Rename some columns?
        dlas.rename_column('Dec', 'DEC')
        dlas.rename_column('logN', 'NHI')

        # Cut on NHI
        gd_dla = dlas['NHI'] >= 20.3

        # Read
        dla_survey = cls.from_sfits(dlas[gd_dla])
        dla_survey.ref = ref

        # g(z) file
        print('Loading QSOs file {:s}'.format(qsos_fil))
        qsos = Table.read(qsos_fil)
        qsos.rename_column('zmin', 'Z_START')
        qsos.rename_column('zmax', 'Z_END')
        qsos.rename_column('Dec', 'DEC')
        qsos.rename_column('zem', 'ZEM')
        dla_survey.sightlines = qsos

        # All?
        if sample == 'all':
            return dla_survey

        # Stat
        # Generate mask
        mask = dla_stat(dla_survey, qsos)
        if sample == 'stat':
            dla_survey.mask = mask
        else:
            dla_survey.mask = ~mask
        # Return
        print('Loaded')
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
        dla_survey = cls.load_lit(dla_fil, qsos_fil, ref, sample=sample)
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
                                  Pdla_fil=Pdla_fil, sample=sample)
        return dla_survey

    @classmethod
    def load_XQ100(cls, sample='stat'):
        """ Load the DLA from XQ-100

        (Sanchez-Ramirez et al. 2016, MNRAS, 456, 4488)

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
                                  sample=sample)
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


def dla_stat(DLAs, qsos, vprox=None, buff=3000.*u.km/u.s,
             zem_min=0., flg_zsrch=0, vmin=0.*u.km/u.s,
             LLS_CUT=None, partial=False, prox=False,
             zem_tol=0.03):
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
    qsos_coord = SkyCoord(ra=qsos['RA'], dec=qsos['DEC'])
    dla_coord = DLAs.coord

    idx, d2d, d3d = match_coordinates_sky(dla_coord, qsos_coord, nthneighbor=1)
    close = d2d < 1.*u.arcsec

    for qq, idla in enumerate(DLAs._abs_sys):
        # In stat?
        if close[qq]:
            if np.abs(idla.zem-qsos['ZEM'][idx[qq]]) < zem_tol:
                if ((idla.zabs >= zmin[idx[qq]]) &
                        (idla.zabs <= qsos['Z_END'][idx[qq]]) & (qsos[idx[qq]]['FLG_BAL'] != 2)):
                        msk_smpl[qq] = True
    # Return
    return msk_smpl
