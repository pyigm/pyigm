"""  Module for the COS-Halos survey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob
import pdb
import warnings
import h5py

from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table import Table, Column

from linetools.spectra import io as lsio
from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.isgm.abscomponent import AbsComponent

from pyigm.metallicity.pdf import MetallicityPDF
from pyigm.cgm.cgmsurvey import CGMAbsSurvey
from pyigm.field.galaxy import Galaxy
from .cgm import CGMAbsSys
from pyigm.abssys.igmsys import IGMSystem

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
    def __init__(self, cdir=None, fits_path=None, kin_init_file=None):
        CGMAbsSurvey.__init__(self)
        self.survey = 'COS-Halos'
        self.ref = 'Tumlinson+11; Werk+12; Tumlinson+13; Werk+13'
        #
        if cdir is None:
            if os.environ.get('COSHALOS_DATA') is None:
                raise ValueError("Need to set COSHALOS_DATA variable")
            self.cdir = os.environ.get('COSHALOS_DATA')
        else:
            self.cdir = cdir
        # Summary Tables
        if fits_path is None:
            self.fits_path = self.cdir+'/Summary/'
        else:
            self.fits_path = fits_path
        try:
            self.cldy = Table.read(self.fits_path+'coshaloscloudysol_newphi.fits')
        except IOError:
            self.cldy = None
        # Kinematics
        self.kin_init_file = self.cdir+'/Kin/coshalo_kin_driver.dat'

    def load_single_fits(self, inp, skip_ions=False, verbose=True):
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
        gal.sfr = (galx['SFR_UPLIM'][0], galx['SFR'][0],
                                       galx['SFR_FLAG'][0]) # FLAG actually gives method used
        # Instantiate the IGM System
        igm_sys = IGMSystem('CGM',(galx['QSORA'][0], galx['QSODEC'][0]),
                            summ['ZFINAL'][0], [-600, 600.]*u.km/u.s)
        igm_sys.zqso = galx['ZQSO'][0]
        # Metallicity
        igm_sys.ZH = -99.
        if self.cldy is not None:
            mtc = np.where((self.cldy['GALID'] == gal.gal_id) &
                           (self.cldy['FIELD'] == gal.field))[0]
            if len(mtc) == 1:
                igm_sys.ZH = self.cldy[mtc]['ZBEST'][0]
        # Instantiate
        cgabs = CGMAbsSys(gal, igm_sys, name=gal.field+'_'+gal.gal_id)
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
        mm = len(self.cgm_abs)-1
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
            for kk in range(iont['NTRANS']):
                flg = iont['FLG'][0][kk]
                if ((flg % 2) == 0) or (flg == 15) or (flg == 13):
                    print('Skipping {:g} as NG for a line'.format(iont['LAMBDA'][0][kk]))
                    continue
                elif (flg == 1) or (flg == 3):
                    flgN = 1
                elif (flg == 5) or (flg == 7):
                    flgN = 3
                elif (flg == 9) or (flg == 11):
                    flgN = 2
                else:
                    pdb.set_trace()
                    raise ValueError("Bad flag!")
                # Fill in
                aline = AbsLine(iont['LAMBDA'][0][kk]*u.AA, closest=True)
                aline.attrib['EW'] = iont['WOBS'][0][kk]*u.AA/1e3  # Observed
                aline.attrib['sig_EW'] = iont['SIGWOBS'][0][kk]*u.AA/1e3
                aline.analy['vlim'] = [iont['VMIN'][0][kk],iont['VMAX'][0][kk]]*u.km/u.s
                aline.attrib['z'] = igm_sys.zabs
                aline.attrib['coord'] = igm_sys.coord
                # Check f
                if (np.abs(aline.data['f']-iont['FVAL'][0][kk])/aline.data['f']) > 0.01:
                    warnings.warn('COS-Halos f-value does not match linetools for {:g}.  Using COS-Halos for now'.format(aline.wrest))
                    aline.data['f'] = iont['FVAL'][0][kk]
                # Colm
                aline.attrib['logN'] = iont['LOGN'][0][kk]
                aline.attrib['sig_logN'] = iont['SIGLOGN'][0][kk]
                aline.attrib['flag_N'] = int(flgN)
                #pdb.set_trace()
                _,_ = ltaa.linear_clm(aline.attrib)
                # Append
                abslines.append(aline)
            # Component
            if len(abslines) == 0:
                comp = AbsComponent(cgabs.igm_sys.coord,
                                    (iont['ZION'][0][0],iont['ZION'][0][1]),
                                    igm_sys.zabs, igm_sys.vlim)

            else:
                comp = AbsComponent.from_abslines(abslines)
            comp.logN = float(iont['CLM'][0])
            comp.sig_logN = float(iont['SIG_CLM'][0])
            comp.flag_N = int(iont['FLG_CLM'][0])
            _,_ = ltaa.linear_clm(comp)
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
        self.cgm_abs[mm].igm_sys._ionN = dat_tab
        # NHI
        self.cgm_abs[mm].igm_sys.NHI = dat_tab[
            (dat_tab['Z']==1)&(dat_tab['ion']==1)]['logN'][0]
        self.cgm_abs[mm].igm_sys.flag_NHI = dat_tab[
            (dat_tab['Z']==1)&(dat_tab['ion']==1)]['flag_N'][0]

    def load_mega(self, data_file=None, cosh_dct=None, test=False, **kwargs):
        """ Load the data for COS-Halos from FITS files taken from the mega structure

        Parameters
        ----------
        data_file : string
          Name of data file
        pckl_fil : string
          Name of file for pickling
        """
        # Loop
        if test is True:
            cos_files = glob.glob(self.fits_path+'/J091*.fits.gz')  # For testing
        else:
            cos_files = glob.glob(self.fits_path+'/J*.fits.gz')
        # Read
        for fil in cos_files:
            self.load_single_fits(fil, **kwargs)

    def load_mtl_pdfs(self, ZH_fil):
        """ Load the metallicity PDFs from an input file (usually hdf5)

        Parameters
        ----------
        ZH_fil : str

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
            cgm_abs.igm_sys.metallicity = MetallicityPDF(fh5['met']['left_edge_bins']+
                                         fh5['met']['left_edge_bins'].attrs['BINSIZE']/2.,
                                         fh5['met'][mkeys[mt][0]])
            cgm_abs.igm_sys.metallicity.inputs = {}
            for key in fh5['inputs'][cgm_abs.name]:
                cgm_abs.igm_sys.metallicity.inputs[key] = fh5['inputs'][cgm_abs.name][key].value

    
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
                mt = np.where( (cgm_abs.field == kin_init['QSO']) &
                               (cgm_abs.gal_id == kin_init['Galaxy']) )[0]
                if len(mt) == 0:
                    print('load_kin: No kinematics for {:s}, {:s}'.format(cgm_abs.field,
                                                                          cgm_abs.gal_id))
                    continue
                mt = mt[0]

                # Metals
                if kin_init['flgL'][mt] > 0:
                    wrest = kin_init['mtl_wr'][mt]*u.AA 
                    if wrest.value <= 1:
                        pdb.set_trace()
                    spec = self.load_bg_cos_spec( qq, wrest )
                    vmnx = (kin_init['L_vmn'][mt]*u.km/u.s, kin_init['L_vmx'][mt]*u.km/u.s)
                    # Process
                    cgm_abs.abs_sys.kin['Metal'] = KinAbs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['Metal'].fill_kin(spec, per=0.07)
                    # Save spec
                    cgm_abs.abs_sys.kin['Metal'].spec = spec
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['Metal'] = KinAbs(0.*u.AA, (0., 0.))

                # HI
                if kin_init['flgH'][mt] > 0:
                    wrest = kin_init['HI_wrest'][mt]*u.AA 
                    if wrest.value <= 1:
                        pdb.set_trace()
                    spec = self.load_bg_cos_spec( qq, wrest )
                    vmnx = (kin_init['HIvmn'][mt]*u.km/u.s, kin_init['HIvmx'][mt]*u.km/u.s) 
                    # Process
                    cgm_abs.abs_sys.kin['HI'] = KinAbs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['HI'].fill_kin(spec, per=0.07)
                    cgm_abs.abs_sys.kin['HI'].spec = spec
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['HI'] = KinAbs(0.*u.AA, (0., 0.))


            #tmp = cos_halos.abs_kin('Metal')['Dv']
            #xdb.set_trace()
    # 
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

        JXP on 12 Oct 2015
        """
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
            spec = specb.splice(specr)
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
            spec = lsio.readspec(slicedir+slicename, flux_tags=['FNORM'], sig_tags=['ENORM'])
        except IOError:
            warnings.warn("File {:s} not found".format(slicedir+slicename))
            return None
        # Fill velocity
        spec.velo = spec.relative_vel((cgm_abs.galaxy.z+1)*wrest)
    
        #spec.qck_plot()
        return spec

    def stack_plot(self, inp, use_lines=None, ymnx=None, **kwargs):
        """ Generate a stack plot of the key lines for a given COS-Halos system
        Parameters
        ----------
        inp : int or tuple
          int -- Index of the cgm_abs list
          tuple -- (field,gal_id)
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

    def __getitem__(self, inp):
        """Grab CgmAbs Class from the list

        Parameters:
        -----------
        ion: tuple
          tuple:  (field,gal_id)

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
        else:
            raise IOError("Bad input")
        # Generate lists
        fields = np.array([cgm_abs.galaxy.field for cgm_abs in self.cgm_abs])
        galids = np.array([cgm_abs.galaxy.gal_id for cgm_abs in self.cgm_abs])
        #
        mt = np.where( (fields == inp[0]) & (galids == inp[1]))[0]
        if len(mt) != 1:
            warnings.warning('CosHalos: CGM not found')
            return None
        else:
            return self.cgm_abs[mt]


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
            self.cdir = os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/'
        if fits_path is None:
            self.fits_path = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/Targets/FITS')
        else:
            self.fits_path = fits_path
        # Kinematics
        if kin_init_file is None:
            self.kin_init_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/Kin/cosdwarfs_kin_driver.dat') 
        else:
            self.kin_init_file = kin_init_file

