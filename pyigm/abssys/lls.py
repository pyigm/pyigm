""" Subclasses for LLS IGMSystem and IGMSurvey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import warnings
import pdb

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.lists.linelist import LineList
from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm.abssystem import AbsSystem
from linetools.isgm import utils as ltiu

from pyigm.abssys.igmsys import IGMSystem, AbsSubSystem
from pyigm.abssys import utils as igmau
from .utils import dict_to_ions

class LLSSystem(IGMSystem):
    """
    Class for an LLS absorption system

    Attributes
    ----------
    tau_ll : float
      Opacity at the Lyman limit
    ZH : float, optional
      Mean metallicity (log10)
    metallicity : MetallicityPDF, optional
    """

    @classmethod
    def from_datfile(cls, dat_file, tree=None, **kwargs):
        """ Read from dat_file (historical JXP format)

        Parameters
        ----------
        dat_file : str
          dat file
        tree : str, optional
          Path to data files
        kwargs :
          Passed to __init__

        Returns
        -------
        _datdict : dict
          Fills this attribute
        """
        if tree is None:
            tree = ''
        # Read datfile
        datdict = igmau.read_dat_file(tree+dat_file)
        # Parse
        coord, zabs, name, NHI, sigNHI, clm_fil = igmau.parse_datdict(datdict)
        kwargs['NHI'] = NHI
        kwargs['sig_NHI'] = sigNHI
        # Generate with type
        vlim = None
        slf = cls(coord, zabs, vlim, **kwargs)

        # Fill files
        slf.tree = tree
        slf.dat_file = slf.tree+dat_file

        # Parse datdict
        #   Includes Sub systems
        slf._datdict = datdict
        slf.parse_dat_file()

        return slf

    @classmethod
    def from_dict(cls, idict, **kwargs):
        """ Generate an LLSSystem from a dict

        Parameters
        ----------
        idict : dict
          Usually read from the hard-drive
        """
        from linetools.isgm.abssystem import add_comps_from_dict, add_other_from_dict
        kwargs = dict(zem=idict['zem'], NHI=idict['NHI'],
                      sig_NHI=idict['sig_NHI'], name=idict['Name'])
        slf = cls(SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg),
                  idict['zabs'], idict['vlim']*u.km/u.s, **kwargs)
        #
        add_other_from_dict(slf, idict)
        add_comps_from_dict(slf, idict, **kwargs)

        # Subsystems
        if 'A' in idict.keys():
            lbls= map(chr, range(65, 91))
            for lbl in lbls:
                if lbl in idict.keys():
                    # Generate
                    subsys = AbsSubSystem.from_dict(slf, idict[lbl], lbl)
                    slf.subsys[lbl] = subsys
                else:
                    pass
            # Total them
            slf.nsub = len(slf.subsys.keys())

        # Return
        return slf

    def __init__(self, radec, zabs, vlim, **kwargs):
        """Standard init
        NHI keyword is required

        Parameters
        ----------
        radec : tuple or coordinate
            RA/Dec of the sightline or astropy.coordinate
        zabs : float
          Absorption redshift
        vlim : Quantity array (2)
          Velocity limits of the system
          Defaulted to +/- 500 km/s if None (see Prochaska et al. 2016 HDLLS)
        NHI= : float, required despite being a keyword
          log10 of HI column density
        **kwargs : keywords
          passed to IGMSystem.__init__
        """
        # NHI
        try:
            NHI = kwargs['NHI']
        except KeyError:
            raise ValueError("NHI must be specified for LLSSystem")
        else:
            kwargs.pop('NHI')
        # vlim
        if vlim is None:
            vlim = [-500.,500.]*u.km/u.s
        # Generate with type
        IGMSystem.__init__(self, radec, zabs, vlim, NHI=NHI, abs_type='LLS', **kwargs)

        # Set tau_LL
        self.tau_LL = (10.**self.NHI)*ltaa.photo_cross(1, 1, 1*u.Ry).to('cm**2').value

        # Other
        self.zpeak = None  # Optical depth weighted redshift
        self.ZH = 0.
        self.metallicity = None  # MetallicityPDF class usually

        # Subsystems
        self.nsub = 0
        self.subsys = {}

    def parse_dat_file(self, vlim=[-300.,300]*u.km/u.s):
        """ Parse the datdict read from the .dat file

        Parameters
        ----------
        vlim : Quantity array (2), optional
          Velocity limits of the subsystems
          Should be pulled from the .clm files
        """

        # LLS keys
        self.bgsrc = self._datdict['QSO name']
        self.zem = float(self._datdict['QSO zem'])  # Was zqso
        self.ZH = float(self._datdict['[M/H] ave'])
        self.nsub = int(self._datdict['N subsys'])
        self.cldyfil = self._datdict['Cloudy Grid File']

        # LLS Subsystems
        if self.nsub > 0:
            lbls= map(chr, range(65, 91))
            # Dict
            keys = (['zabs','NHI','NHIsig','NH','NHsig','log x','sigx','b','bsig','Abund file',
                     'U','Usig','flg_low','flg_alpha','[alpha/H]','sig[a/H]',
                     'flg_Fe','[Fe/H]','sig[Fe/H]','VPFIT file'])
            att = (['zabs','NHI','NHIsig','NH','NHsig','logx','sigx','bval','bsig','clm_file',
                     'U','Usig','flg_low','flg_alpha','alpha_H','sig_a_H',
                     'flg_Fe','Fe_H','sig_Fe_H','VPFIT_file'])
            values = ([0., 0., np.zeros(2), 0., np.zeros(2), 0., np.zeros(2), 0., 0.,
                    '', 0., np.zeros(2), 0, 0, 0., 0., 0, 0., 0., ''])
            null_dict = dict(zip(keys,values))
            # Loop on subsystems
            for i in range(self.nsub):
                # Generate
                zabs = float(self._datdict[lbls[i] + ' zabs'])
                self.subsys[lbls[i]] = AbsSubSystem(self, zabs, vlim, lbls[i])
                self.subsys[lbls[i]]._datdict = {}
                # Fill in dict
                for ii, key in enumerate(keys):
                    try:
                        tmpc = self._datdict[lbls[i]+' '+key]
                    except:
                        raise ValueError('lls_utils: Key "{:s}" not found in {:s}'
                                         .format(lbls[i]+key,self.dat_file))
                    else:  # Convert
                        val = null_dict[key]
                        if val.__class__ == np.ndarray:
                            self.subsys[lbls[i]]._datdict[att[ii]] = np.array(map(float,tmpc.split()))
                        else:  # Single value
                            self.subsys[lbls[i]]._datdict[att[ii]] = (map(type(val),[tmpc]))[0]
                # Set a few special ones as attributes
                self.subsys[lbls[i]].NHI = self.subsys[lbls[i]]._datdict['NHI']
                self.subsys[lbls[i]].sig_NHI = self.subsys[lbls[i]]._datdict['NHIsig']

    def get_ions(self, use_Nfile=False, idict=None, update_zvlim=True,
                 linelist=None, verbose=True):
        """Parse the ions for each Subsystem

        And put them together for the full system
        Fills ._ionN with a QTable

        Parameters
        ----------
        idict : dict, optional
          dict containing the IonClms info
        use_Nfile : bool, optional
          Parse ions from a .clm file (JXP historical)
          NOTE: This ignores velocity constraints on components (i.e. chk_vel=False)
        update_zvlim : bool, optional
          Update zvlim from lines in .clm (as applicable)
        linelist : LineList
        """
        if idict is not None:
            table = dict_to_ions(idict)
            self._ionN = table
        elif use_Nfile:
            # Subsystems
            if self.nsub > 0:  # This speeds things up (but is rarely used)
                linelist = LineList('ISM')
            for lbl in self.subsys.keys():
                clm_fil = self.tree+self.subsys[lbl]._datdict['clm_file']
                # Parse .clm file
                self.subsys[lbl]._clmdict = igmau.read_clmfile(clm_fil, linelist=linelist)
                # Build components from lines
                components = ltiu.build_components_from_dict(self.subsys[lbl]._clmdict,
                                                             coord=self.coord, chk_vel=False)
                self.subsys[lbl]._components = components
                # Update z, vlim
                if update_zvlim:
                    self.update_vlim(sub_system=lbl)
                    self.subsys[lbl].zabs = self.subsys[lbl]._clmdict['zsys']
                # Read .ion file and fill in components
                ion_fil = self.tree+self.subsys[lbl]._clmdict['ion_fil']
                self.subsys[lbl]._indiv_ionclms = igmau.read_ion_file(ion_fil, components)
                # Parse .all file
                all_file = ion_fil.split('.ion')[0]+'.all'
                self.subsys[lbl].all_file=all_file #MF: useful to have
                _ = igmau.read_all_file(all_file, components=components)
                # Build table
                self.subsys[lbl]._ionN = ltiu.iontable_from_components(components,ztbl=self.subsys[lbl].zabs)
                # Add to IGMSystem
                for comp in components:
                    self.add_component(comp)

            # Combine
            if self.nsub == 1:
                self._ionN = self.subsys['A']._ionN
                self._clmdict = self.subsys['A']._clmdict
                #xdb.set_trace()
            elif self.nsub == 0:
                raise ValueError('lls_utils.get_ions: Cannot have 0 subsystems..')
            else:
                self._ionN = self.subsys['A']._ionN
                self._clmdict = self.subsys['A']._clmdict
                warnings.warn('lls_utils.get_ions: Need to update multiple subsystems!! Taking A.')
        else:
            raise ValueError("Need an option in get_ions")

    def fill_lls_lines(self, bval=20.*u.km/u.s, do_analysis=1):
        """
        Generate an HI line list for an LLS.
        Goes into self.lls_lines 

        Now generates a component too.
        Should have it check for an existing HI component..

        Parameters
        ----------
        bval : float, optional
          Doppler parameter in km/s
        do_analysis : int, optional
          flag for analysis
        """
        from linetools.lists import linelist as lll

        # May be replaced by component class (as NT desires)
        HIlines = lll.LineList('HI')

        self.lls_lines = []
        Nval = 10**self.NHI / u.cm**2
        for lline in HIlines._data:
            aline = AbsLine(lline['wrest'], linelist=HIlines)
            # Attributes
            aline.attrib['N'] = Nval
            aline.attrib['b'] = bval
            aline.attrib['z'] = self.zabs
            aline.analy['vlim'] = self.vlim
            aline.analy['do_analysis'] = do_analysis
            aline.attrib['coord'] = self.coord
            self.lls_lines.append(aline)
        # Generate a component (should remove any previous HI)
        self.add_component(AbsComponent.from_abslines(self.lls_lines))

    def flux_model(self, spec, smooth=0):
        """ Generate a LLS model given an input spectrum

        Parameters
        ----------
        spec :  Spectrum1D
        smooth : int, optional
          Number of pixels to smooth by

        Returns
        -------
        model : XSpectrum1D
          Output model is passed back as a Spectrum 
        """
        from linetools.analysis import voigt as lav

        # Energies in LLS rest-frame
        wv_rest = spec.wavelength / (self.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())

        # Get photo_cross and calculate tau
        tau_LL = (10.**self.NHI / u.cm**2) * ltaa.photo_cross(1, 1, energy)

        # Check for lines
        if 'lls_lines' not in self.__dict__.keys():
            self.fill_lls_lines()

        tau_Lyman = lav.voigt_from_abslines(spec.wavelength, self.lls_lines, ret='tau')

        # Combine
        tau_model = tau_LL + tau_Lyman

        # Kludge around the limit
        pix_LL = np.argmin(np.fabs( wv_rest- 911.3*u.AA))
        pix_kludge = np.where((wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA))[0]
        tau_model[pix_kludge] = tau_model[pix_LL]
        
        # Fill in flux
        model = spec.copy()
        model.flux = np.exp(-1. * tau_model).value

        # Smooth?
        if smooth > 0:
            model.gauss_smooth(smooth)

        # Return
        return model

    def load_components(self, inp):
        """ Load components for subsystems from an input object

        May also update/create subsystems

        Parameters
        ----------
        inp : dict or ??
          Input object for loading the components
        """
        if isinstance(inp, dict):
            lbls= map(chr, range(65, 91))
            # Subsystems?
            if 'A' in inp.keys():
                for lbl in lbls:
                    if lbl in inp.keys():
                        if lbl not in self.subsys.keys():
                            self.subsys[lbl] = AbsSubSystem(self,
                                                            inp[lbl]['zsys'],
                                                            [-300., 300]*u.km/u.s,
                                                            lbl)
                        # Fill/update
                        self.subsys[lbl]._clmdict = inp[lbl]  # Not so necessary
                        components = ltiu.build_components_from_dict(self.subsys[lbl]._clmdict,
                                                                 coord=self.coord,
                                                                 chk_vel=True)
                        self.subsys[lbl]._components = components
                        # Update vlim
                        self.update_vlim(sub_system=lbl)
                    else:
                        pass
                self.nsub = len(self.subsys.keys())
            else:
                raise ValueError("Not sure what to do here")
        else:
            raise NotImplementedError("Not ready for this input")

    def get_zpeak(self):
        """ Measure zpeak from an ionic transition
        """
        if self.ions is None:
            print('get_zpeak: Need to fill ions with get_ions first.')
            return

        # Ions for analysis
        low_ions = [ (14,2), (6,2), (13,2), (26,2), (13,3)]  # SiII,CII,AlII,FeII,AlIII
        high_ions= [(14,4), (6,4)]  # SiIV, CIV

        for tt in range(4):
            if tt == 0:
                ions = low_ions
                iflg = 1 # Standard
            elif tt == 1:
                ions = low_ions
                iflg = 2 # Saturated
            elif tt == 2:
                ions = high_ions
                iflg = 1 # Standard
            elif tt == 3:
                ions = high_ions
                iflg = 2 # Standard
            else:
                raise ValueError('Bad value')

            # Search 
            for ion in ions:
                try:
                    t = self.ions[ion]
                except KeyError:
                    continue
                # Measurement?
                if t['flag_N'] == iflg:
                # Identify the transition
                    gdi = np.where( (self.ions.trans['Z'] == ion[0]) &
                                (self.ions.trans['ion'] == ion[1]) &
                                (self.ions.trans['flag_N'] <= iflg) )[0]
                    # Take the first one
                    gdt = self.ions.trans[gdi[0]]
                    wrest = gdt['wrest']
                    flgs = self.clm_analy.clm_lines[wrest].analy['FLAGS']
                    spec_file = self.clm_analy.fits_files[flgs[1] % 64]
                    # Generate an Abs_Line with spectrum
                    line = AbsLine(wrest, z=self.clm_analy.zsys, spec_file=spec_file)
                    # vpeak
                    from linetools import utils as ltu
                    vpeak = line.vpeak()
                    self.zpeak = ltu.z_from_v(self.clm_analy.zsys, vpeak)
                    if tt == 3:
                        print('zpeak WARNING: Using saturated high-ions!!')
                    break
            else:
                continue
            break

        # Error catching
        if self.zpeak is None:
            # Skip primordial LLS
            print('lls.zpeak: No transition in {:s}'.format(self.clm_analy.clm_fil))
            return (0,0), 0.
        # Return
        return ion, vpeak



    # Output
    def __repr__(self):
        return ('<{:s}: {:s} {:s}, zabs={:g}, logNHI={:g}, tau_LL={:g}, [Z/H]={:g} dex>'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                 self.zabs, self.NHI, self.tau_LL, self.ZH))

    def print_abs_type(self):
        """Return a string representing the type of vehicle this is."""
        return 'LLS'


def tau_multi_lls(wave, all_lls, **kwargs):
    """Calculate opacities on an input observed wavelength grid

    Parameters
    ----------
    wave : Quantity array
      Wavelengths
    all_lls : list
      List of LLS Class
    **kwargs : dict
      extra keywords go to lav.voigt_from_abslines

    Returns
    -------
    tau : ndarray
      Optical depth values at input wavelengths
    """
    from linetools.analysis import voigt as lav
    #
    all_tau_model = np.zeros(len(wave))
    # Loop on LLS
    for lls in all_lls:
        # LL
        wv_rest = wave / (lls.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())
        # Get photo_cross and calculate tau
        tau_LL = (10.**lls.NHI / u.cm**2) * ltaa.photo_cross(1,1,energy)

        # Lyman
        tau_Lyman = lav.voigt_from_abslines(wave, lls.lls_lines, ret='tau', **kwargs)
        tau_model = tau_LL + tau_Lyman

        # Kludge around the limit
        pix_LL = np.argmin(np.fabs( wv_rest- 911.3*u.AA ))
        pix_kludge = np.where((wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA))[0]
        tau_model[pix_kludge] = tau_model[pix_LL]
        # Add
        all_tau_model += tau_model
    # Return
    return all_tau_model

