"""  Module for the COS-Halos survey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob
import pdb
import warnings
import h5py
import json, yaml

from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table import Table, Column

from linetools.spectra import io as lsio
from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.analysis import abskin as laak
from linetools.isgm.abscomponent import AbsComponent

from pyigm.metallicity.pdf import MetallicityPDF, DensityPDF, GenericPDF
from pyigm.cgm.cgmsurvey import CGMAbsSurvey
from pyigm.field.galaxy import Galaxy
from .cgm import CGMAbsSys
from pyigm.abssys.igmsys import IGMSystem
import pyigm

class COSHalos(CGMAbsSurvey):
    """Inherits CGM Abs Survey

    Parameters:
    -----------
    fits_path : str, optional
      Path to the FITS data files for COS-Halos
    cdir : str, optional
      Path to the COS-Halos Dropbox
    kin_init_file : str, optional
      Path to kinematics file
    """
    def __init__(self, cdir=None, fits_path=None, load=True, **kwargs):
        CGMAbsSurvey.__init__(self)
        self.survey = 'COS-Halos'
        self.ref = 'Tumlinson+11; Werk+12; Tumlinson+13; Werk+13; Werk+14'
        #
        if cdir is None:
            self.cdir = pyigm.__path__[0]+'/data/CGM/COS_Halos/'
        else:
            self.cdir = cdir
        # Summary Tables
        if fits_path is None:
            self.fits_path = self.cdir+'/Summary/'
        else:
            self.fits_path = fits_path
        try:
            self.werk14_cldy = Table.read(self.fits_path+'coshaloscloudysol_newphi.fits')
        except IOError:
            self.werk14_cldy = None
        # Kinematics
        self.kin_init_file = self.cdir+'/Kin/coshalo_kin_driver.dat'
        # Init?
        if load:
            self.load_sys(**kwargs)

    def load_single_fits(self, inp, skip_ions=False, verbose=True, **kwargs):
        """ Load a single COS-Halos sightline
        Appends to cgm_abs list

        Parameters
        ----------
        inp : tuple or str
          if tuple -- (field,gal_id)
            field: str 
              Name of field (e.g. 'J0226+0015')
            gal_id: str 
              Name of galaxy (e.g. '268_22')
        skip_ions : bool, optional
          Avoid loading the ions (not recommended)
        verbose : bool, optional
        """
        # Parse input
        if isinstance(inp, basestring):
            fil = inp
        elif isinstance(inp, tuple):
            field, gal_id = inp
            tmp = self.fits_path+'/'+field+'.'+gal_id+'.fits.gz'
            fils = glob.glob(tmp)
            if len(fils) != 1:
                raise IOError('Bad field, gal_id: {:s}'.format(tmp))
            fil = fils[0]
        else:
            raise IOError('Bad input to load_single')

        # Read COS-Halos file
        if verbose:
            print('cos_halos: Reading {:s}'.format(fil))
        hdu = fits.open(fil)
        summ = Table(hdu[1].data)
        galx = Table(hdu[2].data)
        # Instantiate the galaxy
        gal = Galaxy((galx['RA'][0], galx['DEC'][0]), z=summ['ZFINAL'][0])
        gal.field = galx['FIELD'][0]
        gal.gal_id = galx['GALID'][0]
        # Galaxy properties
        gal.halo_mass = summ['LOGMHALO'][0]
        gal.stellar_mass = summ['LOGMFINAL'][0]
        gal.rvir = galx['RVIR'][0]
        gal.MH = galx['ABUN'][0]
        gal.flag_MH = galx['ABUN_FLAG'][0]
        gal.sdss_phot = [galx[key][0] for key in ['SDSSU','SDSSG','SDSSR','SDSSI','SDSSZ']]
        gal.sdss_phot_sig = [galx[key][0] for key in ['SDSSU_ERR','SDSSG_ERR','SDSSR_ERR','SDSSI_ERR','SDSSZ_ERR']]
        gal.sfr = (galx['SFR_UPLIM'][0], galx['SFR'][0],
                                       galx['SFR_FLAG'][0]) # FLAG actually gives method used
        gal.ssfr = galx['SSFR'][0]
        # Instantiate the IGM System
        igm_sys = IGMSystem((galx['QSORA'][0], galx['QSODEC'][0]),
                            summ['ZFINAL'][0], [-600, 600.]*u.km/u.s,
                            abs_type='CGM')
        igm_sys.zqso = galx['ZQSO'][0]
        # Instantiate
        cgabs = CGMAbsSys(gal, igm_sys, name=gal.field+'_'+gal.gal_id, **kwargs)
        # EBV
        cgabs.ebv = galx['EBV'][0]
        # Ions
        if skip_ions is True:
            # NHI
            dat_tab = Table(hdu[3].data)
            #if dat_tab['Z'] != 1:
            #    raise ValueError("Uh oh")
            cgabs.igm_sys.NHI = dat_tab['CLM'][0]
            cgabs.igm_sys.sig_NHI = dat_tab['SIG_CLM'][0]
            cgabs.igm_sys.flag_NHI = dat_tab['FLG_CLM'][0]
            self.cgm_abs.append(cgabs)
            return
        all_Z = []
        all_ion = []
        for jj in range(summ['NION'][0]):
            iont = hdu[3+jj].data
            if jj == 0: # Generate new Table
                dat_tab = Table(iont)
            else:
                try:
                    dat_tab.add_row(Table(iont)[0])
                except:
                    pdb.set_trace()
            all_Z.append(iont['ZION'][0][0])
            all_ion.append(iont['ZION'][0][1])
            # AbsLines
            abslines = []
            ntrans = len(np.where(iont['LAMBDA'][0] > 1.)[0])
            for kk in range(ntrans):
                flg = iont['FLG'][0][kk]
                # Fill in
                aline = AbsLine(iont['LAMBDA'][0][kk]*u.AA, closest=True)
                aline.attrib['flag_origCH'] = int(flg)
                aline.attrib['EW'] = iont['WOBS'][0][kk]*u.AA/1e3  # Observed
                aline.attrib['sig_EW'] = iont['SIGWOBS'][0][kk]*u.AA/1e3
                if aline.attrib['EW'] > 3.*aline.attrib['sig_EW']:
                    aline.attrib['flag_EW'] = 1
                else:
                    aline.attrib['flag_EW'] = 3
                # Force an upper limit (i.e. from a blend)
                if (flg == 2) or (flg == 4) or (flg == 6):
                    aline.attrib['flag_EW'] = 3
                #
                aline.analy['vlim'] = [iont['VMIN'][0][kk],iont['VMAX'][0][kk]]*u.km/u.s
                aline.attrib['z'] = igm_sys.zabs
                aline.attrib['coord'] = igm_sys.coord
                # Check f
                if (np.abs(aline.data['f']-iont['FVAL'][0][kk])/aline.data['f']) > 0.001:
                    Nscl = iont['FVAL'][0][kk] / aline.data['f']
                    flag_f = True
                else:
                    Nscl = 1.
                    flag_f = False
                # Colm
                if ((flg % 2) == 0) or (flg == 15) or (flg == 13):
                    flgN = 0
                    print('Skipping column contribution from {:g} as NG for a line; flg={:d}'.format(iont['LAMBDA'][0][kk],flg))
                elif (flg == 1) or (flg == 3):
                    flgN = 1
                elif (flg == 5) or (flg == 7):
                    flgN = 3
                elif (flg == 9) or (flg == 11):
                    flgN = 2
                else:
                    pdb.set_trace()
                    raise ValueError("Bad flag!")
                if flgN == 3:
                    aline.attrib['logN'] = iont['LOGN2SIG'][0][kk] + np.log10(Nscl)
                    aline.attrib['sig_logN'] = 9.
                elif flgN == 0:  # Not for N measurement
                    pass
                else:
                    aline.attrib['logN'] = iont['LOGN'][0][kk] + np.log10(Nscl)
                    aline.attrib['sig_logN'] = iont['SIGLOGN'][0][kk]
                aline.attrib['flag_N'] = int(flgN)
                #pdb.set_trace()
                if flgN != 0:
                    _,_ = ltaa.linear_clm(aline.attrib)
                # Append
                abslines.append(aline)
            # Component
            if len(abslines) == 0:
                comp = AbsComponent(cgabs.igm_sys.coord,
                                    (iont['ZION'][0][0],iont['ZION'][0][1]),
                                    igm_sys.zabs, igm_sys.vlim)

            else:
                comp = AbsComponent.from_abslines(abslines, chk_vel=False)
                if comp.Zion != (1,1):
                    comp.synthesize_colm()  # Combine the abs lines
                    if np.abs(comp.logN - float(iont['CLM'][0])) > 0.15:
                        print("New colm for ({:d},{:d}) and sys {:s} is {:g} different from old".format(
                            comp.Zion[0], comp.Zion[1], cgabs.name, comp.logN - float(iont['CLM'][0])))
                    if comp.flag_N != iont['FLG_CLM'][0]:
                        if comp.flag_N == 0:
                            pass
                        else:
                            print("New flag for ({:d},{:d}) and sys {:s} is different from old".format(
                                comp.Zion[0], comp.Zion[1], cgabs.name))
                            pdb.set_trace()
            #_,_ = ltaa.linear_clm(comp)
            cgabs.igm_sys.add_component(comp)
        self.cgm_abs.append(cgabs)

        # Add Z,ion
        dat_tab.add_column(Column(all_Z,name='Z'))
        dat_tab.add_column(Column(all_ion,name='ion'))
        # Rename
        dat_tab.rename_column('LOGN','indiv_logN')
        dat_tab.rename_column('SIGLOGN','indiv_sig_logN')
        dat_tab.rename_column('CLM','logN')
        dat_tab.rename_column('SIG_CLM','sig_logN')
        dat_tab.rename_column('FLG_CLM','flag_N')
        # Set
        self.cgm_abs[-1].igm_sys._ionN = dat_tab
        # NHI
        HI = (dat_tab['Z'] == 1) & (dat_tab['ion'] == 1)
        if np.sum(HI) > 0:
            self.cgm_abs[-1].igm_sys.NHI = dat_tab[HI]['logN'][0]
            self.cgm_abs[-1].igm_sys.sig_NHI = dat_tab[HI]['sig_logN'][0]
            self.cgm_abs[-1].igm_sys.flag_NHI = dat_tab[HI]['flag_N'][0]
        else:
            warnings.warn("No HI measurement for {}".format(self.cgm_abs[-1]))
            self.cgm_abs[-1].igm_sys.flag_NHI = 0

        #if self.cgm_abs[-1].name == 'J0950+4831_177_27':
        #    pdb.set_trace()

    def load_mega(self, data_file=None, cosh_dct=None, test=False, **kwargs):
        """ Load the data for COS-Halos from FITS files taken from the mega structure

        Parameters
        ----------
        data_file : string
          Name of data file
        pckl_fil : string
          Name of file for pickling
        """
        warnings.warn("This method will be DEPRECATED")
        # Loop
        if test is True:
            cos_files = glob.glob(self.fits_path+'/J09*.fits.gz')  # For testing
        elif 'Dwarfs' in self.fits_path:  # COS-Dwarfs
            cos_files = glob.glob(self.fits_path+'/*.fits')
        else:  # COS-Halos
            cos_files = glob.glob(self.fits_path+'/J*.fits.gz')
        # Read
        for fil in cos_files:
            self.load_single_fits(fil, **kwargs)
        # Werk+14
        if ('Halos' in self.fits_path) and (self.werk14_cldy is not None):
            self.load_werk14()

    def load_werk14(self):
        """ Load up the Werk+14 results
        """
        cldy_names = np.array([row['FIELD'].strip()+'_'+row['GALID'].strip() for row in self.werk14_cldy])
        for cgm_abs in self.cgm_abs:
            igm_sys = cgm_abs.igm_sys
            # Metallicity
            igm_sys.ZH = -99.
            #mtc = np.where((self.werk14_cldy['GALID'] == gal.gal_id) &
            #               (self.werk14_cldy['FIELD'] == gal.field))[0]
            mtc = np.where(cldy_names == cgm_abs.name)[0]
            if len(mtc) == 1:
                # Metallicity
                igm_sys.werk14_ZH = self.werk14_cldy['ZBEST'][mtc][0]
                igm_sys.werk14_ZHmnx = [self.werk14_cldy['ZMIN'][mtc][0],
                                        self.werk14_cldy['ZMAX'][mtc][0]]
                igm_sys.ZH = igm_sys.werk14_ZH
                # NHI
                igm_sys.werk14_NHI = self.werk14_cldy['NHI_BEST'][mtc][0]
                # NH
                igm_sys.werk14_NH = np.log10(self.werk14_cldy['NH_BEST'][mtc][0])
                igm_sys.werk14_NHmnx = [np.log10(self.werk14_cldy['NH_LOW'][mtc][0]),
                                        np.log10(self.werk14_cldy['NH_HIGH'][mtc][0])]
            else:
                print('No Werk+14 Cloudy solution for {:s}'.format(cgm_abs.name))
            # Hand edit for J0914+2823_41_27
            if cgm_abs.name == 'J0914+2823_41_27':
                igm_sys.werk14_ZH = -0.8

    def load_sys(self, tfile=None, empty=True, debug=False, **kwargs):
        """ Load the COS-Halos survey from JSON files

        Empties the list

        Parameters
        ----------
        tfile : str, optional
        empty : bool, optional
          Empty the list
        debug : bool, optional
          Only load the first 5

        Returns
        -------

        """
        import tarfile
        import json
        from linetools.lists.linelist import LineList
        llist = LineList('ISM')

        # Tar file
        if tfile is None:
            tarfiles = glob.glob(self.cdir+'cos-halos_systems.v*.tar.gz')
            tarfiles.sort()
            tfile = tarfiles[-1]
        print("Be patient, using {:s} to load".format(tfile))
        # Empty
        if empty:
            self.cgm_abs = []
        # Load
        tar = tarfile.open(tfile)
        for kk, member in enumerate(tar.getmembers()):
            if '.' not in member.name:
                print('Skipping a likely folder: {:s}'.format(member.name))
                continue
            # Debug
            if debug and (kk == 5):
                break
            # Extract
            f = tar.extractfile(member)
            tdict = json.load(f)
            # Generate
            cgmsys = CGMAbsSys.from_dict(tdict, chk_vel=False, chk_sep=False, chk_data=False,
                                         use_coord=True, use_angrho=True,
                                         linelist=llist, **kwargs)
            self.cgm_abs.append(cgmsys)
        tar.close()
        # Werk+14
        if ('Halos' in self.fits_path) and (self.werk14_cldy is not None):
                self.load_werk14()

    def load_mtl_pdfs(self, ZH_fil, keep_all=False):
        """ Load the metallicity PDFs from an input file (usually hdf5)

        Parameters
        ----------
        ZH_fil : str
        keep_all : bool, optional
          Save the full metallicity data?

        """
        fh5=h5py.File(ZH_fil, 'r')
        mkeys = fh5['met'].keys()
        mkeys.remove('left_edge_bins')
        mkeys.remove('right_edge_bins')
        mkeys = np.array(mkeys)

        # Loop
        for cgm_abs in self.cgm_abs:
            # Match?
            mt = np.where(mkeys == cgm_abs.name)[0]
            if len(mt) == 0:
                print('No metallicity info for {:s}'.format(cgm_abs.name))
                print('Skipping..')
                continue
            # Z/H
            cgm_abs.igm_sys.metallicity = MetallicityPDF(fh5['met']['left_edge_bins']+
                                         fh5['met']['left_edge_bins'].attrs['BINSIZE']/2.,
                                         fh5['met'][mkeys[mt][0]])
            cgm_abs.igm_sys.metallicity.inputs = {}
            for key in fh5['inputs'][cgm_abs.name]:
                cgm_abs.igm_sys.metallicity.inputs[key] = fh5['inputs'][cgm_abs.name][key].value
            # NHI
            cgm_abs.igm_sys.NHIPDF = GenericPDF(fh5['col']['left_edge_bins']+
                                                 fh5['col']['left_edge_bins'].attrs['BINSIZE']/2.,
                                                 fh5['col'][mkeys[mt][0]])
            # Density
            cgm_abs.igm_sys.density = DensityPDF(fh5['dens']['left_edge_bins']+
                                                         fh5['dens']['left_edge_bins'].attrs['BINSIZE']/2.,
                                                         fh5['dens'][mkeys[mt][0]])

    
    ########################## ##########################
    def load_abskin(self, flg=1, kin_file=None, kin_init_file=None):
        """ Load the absorption-line kinematic data for COS-Halos (or COS-Dwarfs)
        Calculate from scratch if needed

        Parameters
        ----------
        flg: int, optional 
          Flag indicating how to load the data
            0 = Load from file
            1 = Generate
        kin_init_file: str
          Name of kinematics driver file
        kin_file: str
          Name of kinematics output file [First made for John Forbes]
        """
    
        if flg == 0: # Load
            if kin_file is None:
                kin_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/Kin/'+
                                                  'COS-Halos_kin.fits')
            hdu = fits.open(kin_file)
            # Metals
            metals = Table(hdu[1].data)
            for row in metals:
                mt = np.where( (row['field']==self.field) & 
                    (row['gal_id']==self.gal_id))[0]
                pdb.set_trace()

        elif flg == 1: # Generate
            # Read init file
            if kin_init_file is None:
                kin_init_file = self.kin_init_file
            kin_init = ascii.read(kin_init_file,guess=False)
    
            # Loop to my loop
            fgal = zip(self.field, self.gal_id)
            for qq,cgm_abs in enumerate(self.cgm_abs):
                # Match to kin_init
                mt = np.where((kin_init['QSO'] == cgm_abs.field) &
                               (kin_init['Galaxy'] == cgm_abs.gal_id))[0]
                if len(mt) == 0:
                    warnings.warn('load_kin: No kinematics for {:s}, {:s}'.format(cgm_abs.field,
                                                                          cgm_abs.gal_id))
                    continue
                mt = mt[0]

                # Metals
                if kin_init['flgL'][mt] > 0:
                    wrest = kin_init['mtl_wr'][mt]*u.AA 
                    if wrest.value <= 1:
                        pdb.set_trace()
                    spec = self.load_bg_cos_spec(qq, wrest)
                    vmnx = (kin_init['L_vmn'][mt], kin_init['L_vmx'][mt])*u.km/u.s
                    # Process
                    aline = AbsLine(wrest)
                    aline.analy['spec'] = spec
                    aline.analy['vlim'] = vmnx
                    aline.attrib['z'] = cgm_abs.igm_sys.zabs
                    fx, sig, cdict = aline.cut_spec()
                    # Kin
                    stau = laak.generate_stau(cdict['velo'], fx, sig)
                    cgm_abs.igm_sys.kin['Metal'] = laak.pw97_kin(cdict['velo'], stau, per=0.07)
                    cgm_abs.igm_sys.kin['Metal'].update(laak.cgm_kin(cdict['velo'], stau, per=0.07))
                    # Save spec
                    #cgm_abs.igm_sys.kin['Metal']['spec'] = spec
                else:
                    cgm_abs.igm_sys.kin['Metal'] = {}

                # HI
                if kin_init['flgH'][mt] > 0:
                    wrest = kin_init['HI_wrest'][mt]*u.AA 
                    if wrest.value <= 1:
                        pdb.set_trace()
                    spec = self.load_bg_cos_spec( qq, wrest )
                    if spec is None:
                        pdb.set_trace()
                    vmnx = (kin_init['HIvmn'][mt], kin_init['HIvmx'][mt])*u.km/u.s
                    # Process
                    aline = AbsLine(wrest)
                    aline.analy['spec'] = spec
                    aline.analy['vlim'] = vmnx
                    aline.attrib['z'] = cgm_abs.igm_sys.zabs
                    fx, sig, cdict = aline.cut_spec()
                    # Kin
                    stau = laak.generate_stau(cdict['velo'], fx, sig)
                    cgm_abs.igm_sys.kin['HI'] = laak.pw97_kin(cdict['velo'], stau, per=0.07)
                    cgm_abs.igm_sys.kin['HI'].update(laak.cgm_kin(cdict['velo'], stau, per=0.07))
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.igm_sys.kin['HI'] = {}

    def load_gal_spec(self, inp):
        """ Load the galaxy spectrum

        Parameters
        ----------
        inp : int or tuple
          int -- Index of the cgm_abs list
          tuple -- (field,gal_id)

        Returns
        ----------
        spec : XSpectrum1D
          Splices the blue and red side for LRIS
        """
        from linetools.spectra import utils as ltsu
        # Init
        cgm_abs = self[inp]
        # Directories
        galdir = self.cdir+'/Galaxies/'
        #fielddir = 'fields/'+cgm_abs.field+'/'
        #sysdir = cgm_abs.gal_id+'/spec1d/'
        sysname = cgm_abs.galaxy.field+'_'+cgm_abs.galaxy.gal_id

        # Find files
        lris_files = glob.glob(galdir+sysname+'*corr.fits.gz')
        if len(lris_files) == 0:
            raise ValueError('No LRIS files! {:s}'.format(galdir+sysname))
        elif len(lris_files) == 2:
            lris_files.sort()
            specb = lsio.readspec(lris_files[0]) 
            specr = lsio.readspec(lris_files[1]) 
            spec = ltsu.splice_two(specb, specr)
        else:
            raise ValueError('Not sure what happened')

        # Return
        return spec


    def load_bg_cos_spec(self, inp, wrest):
        """ Load the absorption-line kinematic data for COS-Halos
        Calculate from scratch if needed

        Parameters
        ----------
        inp : int or tuple
          int -- Index of the cgm_abs list
          tuple -- (field,gal_id)
        wrest : Quantity
          Rest wavelength for spectrum of interest
    
        JXP on 11 Dec 2014
        """
        cgm_abs = self[inp]
        # Directories
        sysdir = cgm_abs.galaxy.gal_id+'_z{:5.3f}'.format(cgm_abs.galaxy.z)
        sysname = cgm_abs.galaxy.field+'_'+sysdir

        # Transition
        templ_fil = self.cdir+'/Targets/system_template.lst'
        tab = ascii.read(templ_fil)
        mt = np.argmin(np.abs(tab['col1']-wrest.value))
        if np.abs(tab['col1'][mt]-wrest.value) > 1e-2:
            raise ValueError('get_coshalo_spec: wrest={:g} not found!'.format(wrest))
        trans = tab['col2'][mt]+tab['col3'][mt]

        # Read
        slicedir = self.cdir+'/Targets/fitting/'
        slicename = sysname+'_'+trans+'_slice.fits'
        try:
            spec = lsio.readspec(slicedir+slicename,
                                 flux_tag='FNORM', sig_tag='ENORM')
        except IOError:
            warnings.warn("File {:s} not found".format(slicedir+slicename))
            return None
        # Fill velocity
        spec.velo = spec.relative_vel((cgm_abs.galaxy.z+1)*wrest)
    
        #spec.qck_plot()
        return spec

    def stack_plot(self, inp, use_lines=None, ymnx=None, add_lines=None, **kwargs):
        """ Generate a stack plot of the key lines for a given COS-Halos system
        Parameters
        ----------
        inp : int or tuple
          int -- Index of the cgm_abs list
          tuple -- (field,gal_id)
        add_lines : list, optional
          List of additional lines to plot
        """
        # Init
        from linetools.analysis import plots as ltap
        if ymnx is None:
            ymnx=(-0.1,1.2)
        cgm_abs = self[inp]
        abs_lines = []
        # Setup the lines (defaults to a key seto)
        if use_lines is None:
            use_lines = [1215.6700, 1025.7223, 1334.5323, 977.020, 1031.9261, 1037.6167,
                         1260.4221, 1206.500, 1393.7550, 2796.352]*u.AA
            if add_lines is not None:
                use_lines = list(use_lines.value) + add_lines
                use_lines.sort()
                use_lines = use_lines*u.AA
        for iline in use_lines:
            spec = self.load_bg_cos_spec(inp, iline)
            if spec is None:
                print('Skipping {:g}. Assuming no coverage'.format(iline))
            aline = AbsLine(iline, closest=True)
            aline.analy['spec'] = spec
            aline.attrib['z'] = cgm_abs.galaxy.z
            abs_lines.append(aline)
        # Execute
        ltap.stack_plot(abs_lines, vlim=[-400., 400]*u.km/u.s, ymnx=ymnx, **kwargs)

    def write_survey(self, outfil='COS-Halos_sys.tar.gz'):
        """ Write the survey to a tarball of JSON files

        Parameters
        ----------
        outfil : str, optional
        """
        self.to_json_tarball(outfil)

    def __getitem__(self, inp):
        """Grab CgmAbs Class from the list

        Parameters:
        -----------
        ion: tuple
          tuple:  (field,gal_id)
          str: field_gal_id

        Returns:
        ----------
        cgm_abs
        """
        if isinstance(inp,int):
            return self.cgm_abs[inp]
        elif isinstance(inp,tuple):
            if not isinstance(inp[0], basestring):
                raise IOError("Bad input")
            if not isinstance(inp[1], basestring):
                raise IOError("Bad input")
            field = inp[0]
            galid = inp[1]
        elif isinstance(inp, basestring):
            # Split at the first _
            under = inp.find('_')
            field = inp[:under]
            galid = inp[under+1:]
        else:
            raise IOError("Bad input")
        # Generate lists
        fields = np.array([cgm_abs.galaxy.field for cgm_abs in self.cgm_abs])
        galids = np.array([cgm_abs.galaxy.gal_id for cgm_abs in self.cgm_abs])
        #
        mt = np.where( (fields == field) & (galids == galid))[0]
        if len(mt) == 0:
            warnings.warn('CosHalos: CGM not found')
            return None
        elif len(mt) == 1:
            return self.cgm_abs[mt[0]]
        else:
            raise ValueError("Multiple hits.  Should not happen")


def update_cos_halos(v10=False, v11=True):
    """ Updates the JSON tarball
    Returns
    -------
    """
    print("See the COS-Halos_database notebook for details")

    # Generate v1.0
    if v10:
        print("Generate v1.0 of the JSON tarball")
        cos_halos = COSHalos()
        cos_halos.load_mega()
        tarfil = cos_halos.cdir+'/cos-halos_systems.v1.0.tar.gz'
        cos_halos.write_survey(tarfil)
        del cos_halos

    # Generate v1.1 which uses the NHI measurements from P+16 and
    #   modifies a few ions and transitions as limits or NG
    if v11:
        print("Generate v1.1 of the JSON tarball")
        cos_halos_v10 = COSHalos(load=False)
        tfile = pyigm.__path__[0]+'/data/CGM/COS_Halos/cos-halos_systems.v1.0.tar.gz'
        cos_halos_v10.load_sys(tfile=tfile)
        # NHI
        LLS_file = pyigm.__path__[0]+'/data/CGM/COS_Halos/COS_Halos_LLS.json'
        with open(LLS_file) as json_file:
            fdict = json.load(json_file)
        # Loop on systems
        names = cos_halos_v10.name
        for key in fdict.keys():
            # Match
            mt = np.where(names == key)[0]
            if len(mt) != 1:
                raise ValueError("No match?!")
            # Fill in
            if fdict[key]['flag_NHI'] == 1:
                cos_halos_v10.cgm_abs[mt[0]].igm_sys.NHI = fdict[key]['fit_NHI'][0]
                cos_halos_v10.cgm_abs[mt[0]].igm_sys.sig_NHI = [fdict[key]['fit_NHI'][0]-fdict[key]['fit_NHI'][1],
                                                                fdict[key]['fit_NHI'][2]-fdict[key]['fit_NHI'][0]]
                cos_halos_v10.cgm_abs[mt[0]].igm_sys.flag_NHI = 1
            elif fdict[key]['flag_NHI'] == 2:
                cos_halos_v10.cgm_abs[mt[0]].igm_sys.NHI = fdict[key]['fit_NHI'][1]
                cos_halos_v10.cgm_abs[mt[0]].igm_sys.flag_NHI = 2
            elif fdict[key]['flag_NHI'] == 3:
                cos_halos_v10.cgm_abs[mt[0]].igm_sys.NHI = fdict[key]['fit_NHI'][2]
                cos_halos_v10.cgm_abs[mt[0]].igm_sys.flag_NHI = 3
        # Bad ions
        filename = pyigm.__path__[0]+'/data/CGM/COS_Halos/COS_Halos_ions_updates_v1.0.yaml'
        with open(filename, 'r') as infile:
                up_ion_data = yaml.load(infile)
        for key in up_ion_data.keys():
            # Match
            mt = np.where(names == key)[0]
            if len(mt) != 1:
                raise ValueError("No match?!")
            igm_sys = cos_halos_v10.cgm_abs[mt[0]].igm_sys
            # Fill in
            for mod_type in up_ion_data[key].keys():
                if mod_type == 'ion':
                    for ionkey in up_ion_data[key][mod_type].keys():
                        Zion = tuple([int(ii) for ii in ionkey.split(',')])
                        #
                        Zicomp = [comp.Zion for comp in igm_sys._components]
                        mtZi = Zicomp.index(Zion)
                        # Set
                        for att_key in up_ion_data[key][mod_type][ionkey].keys():
                            if att_key == 'flag_N':
                                #cos_halos_v10.cgm_abs[mt[0]].igm_sys._components[mtZi].flag_N = up_ion_data[key][mod_type][ionkey][att_key]
                                igm_sys._components[mtZi].flag_N = up_ion_data[key][mod_type][ionkey][att_key]
                            else:
                                raise ValueError("Bad key for attribute")
                        print(cos_halos_v10.cgm_abs[mt[0]].igm_sys._components[mtZi])
                elif mod_type == 'trans':
                    for transkey in up_ion_data[key][mod_type].keys():
                        # Update AbsLine
                        lines = igm_sys.list_of_abslines()
                        trans = [iline.name for iline in lines]
                        aline = lines[trans.index(transkey)]
                        comp = igm_sys.get_comp_from_absline(aline)  # Grab it now before it changes
                        if att_key == 'flag_N':
                            aline.attrib['flag_N'] = up_ion_data[key][mod_type][transkey][att_key]
                        # Remake component
                        try:
                            comp.synthesize_colm(overwrite=True)
                        except ValueError:
                            pdb.set_trace()
                else:
                    raise ValueError("Bad mod_type")
        # Metallicity
        mtlfil = pyigm.__path__[0]+'/data/CGM/COS_Halos/COS_Halos_MTL_final.hdf5'
        cos_halos_v10.load_mtl_pdfs(mtlfil)
        #
        for cgm_abs in cos_halos_v10.cgm_abs:
            if hasattr(cgm_abs.igm_sys, 'metallicity'):
                cgm_abs.igm_sys.ZH = cgm_abs.igm_sys.metallicity.medianZH
                cgm_abs.igm_sys.sig_ZH = cgm_abs.igm_sys.metallicity.confidence_limits(0.68)
        # Write
        tarfil = pyigm.__path__[0]+'/data/CGM/COS_Halos/cos-halos_systems.v1.1.tar.gz'
        cos_halos_v10.write_survey(tarfil)


class COSDwarfs(COSHalos):
    """Inherits COS Halos Class

    Attributes:
    -----------
    fits_path: str, optional
      Path to the FITS data files for COS-Halos
    """
    # Initialize with a .dat file
    def __init__(self, tree=None, fits_path=None, kin_init_file=None, cdir=None):

        # Generate with type
        CGMAbsSurvey.__init__(self)
        self.survey = 'COS-Dwarfs'
        self.ref = 'Bordoloi+14'
        if cdir is None:
            self.cdir = os.environ.get('COSHALOS_DATA')
            #self.cdir = os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/'
        if fits_path is None:
            self.fits_path = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/Targets/FITS')
        else:
            self.fits_path = fits_path
        # Kinematics
        if kin_init_file is None:
            #self.kin_init_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/Kin/cosdwarfs_kin_driver.dat')
            self.kin_init_file = self.cdir+'/Kin/cosdwarfs_kin_driver.dat'
        else:
            self.kin_init_file = kin_init_file

