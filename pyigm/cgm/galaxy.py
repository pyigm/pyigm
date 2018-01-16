"""  Module for CGM of the Milky Way
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from pkg_resources import resource_filename

from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord

from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm import utils as ltiu
from linetools.lists.linelist import LineList
from linetools.analysis.absline import linear_clm

from pyigm.abssys.igmsys import IGMSystem
from pyigm.field.galaxy import Galaxy
from pyigm.cgm.cgm import CGM, CGMAbsSys
from pyigm.cgm.cgmsurvey import CGMAbsSurvey
from pyigm.utils import calc_Galactic_rho

c_kms = const.c.to('km/s').value

class GalaxyCGM(CGM):
    """Inherits CGM class

    Parameters:
    -----------
    load : bool, optional
      Load datasets
    verbose : bool, optional
    """
    def __init__(self, load=True, verbose=False, debug=False, **kwargs):
        # Init
        self.verbose = verbose
        self.debug = debug
        # Generate the Milky Way
        milkyway = Galaxy((0.,0.), z=0.)
        self.galaxy = milkyway
        CGM.__init__(self, milkyway)
        self.refs = ''
        # Absorption
        self.abs = CGMAbsSurvey(survey='Galaxy')
        # Load Data
        if load:
            print("Loading data.  This takes ~20s to build it all...")
            self.load_coolgas()  # This needs to be first!
            self.load_hotgas()

    def load_coolgas(self):
        """ Load data on cool gas (CII, CIV, SiII, SiIII)
        Richter+17
        """
        llist = LineList('ISM')
        # Ricther+17
        print('Loading Richter+17 for CII, CIV, SiII, SiIII')
        r17_a1_file = resource_filename('pyigm','/data/CGM/Galaxy/richter17_A1.fits')
        r17_a1 = Table.read(r17_a1_file)
        r17_a2_file = resource_filename('pyigm','/data/CGM/Galaxy/richter17_A2.fits')
        r17_a2 = Table.read(r17_a2_file)
        # Coords
        coords = SkyCoord(ra=r17_a1['_RAJ2000'], dec=r17_a1['_DEJ2000'], unit='deg')
        gc = coords.transform_to('galactic')
        ra = np.zeros((len(r17_a2)))
        dec = np.zeros((len(r17_a2)))

        # Loop on Sightlines
        for kk,row in enumerate(r17_a1):
            if self.debug and (kk == 5):
                break
            a2_idx = np.where(r17_a2['Name'] == row['Name'])[0]
            if len(a2_idx) == 0:
                continue
            ra[a2_idx] = row['_RAJ2000']
            dec[a2_idx] = row['_DEJ2000']
            # Generate the components
            icoord = gc[kk]
            alines = []
            for jj,idx in enumerate(a2_idx):
                # Transition
                trans = '{:s} {:d}'.format(r17_a2['Ion'][idx].strip(), int(r17_a2['lambda0'][idx]))
                try:
                    aline = AbsLine(trans, linelist=llist)
                except ValueError:
                    pdb.set_trace()
                aline.attrib['coord'] = icoord
                # Velocity
                z = 0.
                aline.setz(z)
                vlim = np.array([r17_a2['vmin'][idx], r17_a2['vmax'][idx]]) * u.km / u.s
                aline.limits.set(vlim)
                # EW
                aline.attrib['flag_EW'] = 1
                aline.attrib['EW'] = r17_a2['W'][idx] / 1e3 * u.AA
                aline.attrib['sig_EW'] = r17_a2['e_W'][idx] / 1e3 * u.AA
                # Column
                if r17_a2['l_logN'][idx] == '>':
                    aline.attrib['flag_N'] = 2
                    aline.attrib['sig_logN'] = 99.99
                else:
                    aline.attrib['flag_N'] = 1
                    aline.attrib['sig_logN'] = r17_a2['e_logN'][idx]
                aline.attrib['logN'] = r17_a2['logN'][idx]
                # Fill linear
                _, _ = linear_clm(aline.attrib)
                alines.append(aline)
            # Generate components from abslines
            comps = ltiu.build_components_from_abslines(alines, chk_sep=False, chk_vel=False)
            # Limits
            vmin = np.min([icomp.limits.vmin.value for icomp in comps])
            vmax = np.max([icomp.limits.vmax.value for icomp in comps])
            # Instantiate
            s_kwargs = dict(name=row['Name'] + '_z0')
            c_kwargs = dict(chk_sep=False, chk_z=False)
            abssys = IGMSystem.from_components(comps, vlim=[vmin,vmax]*u.km/u.s, s_kwargs=s_kwargs, c_kwargs=c_kwargs)
            # CGM Abs
            rho, ang_sep = calc_Galactic_rho(abssys.coord)
            cgmabs = CGMAbsSys(self.galaxy, abssys, rho=rho, ang_sep=ang_sep, cosmo=self.cosmo)
            # Add to cgm_abs
            self.abs.cgm_abs.append(cgmabs)
        # Finish
        r17_a2['RA'] = ra
        r17_a2['DEC'] = dec
        self.richter17 = r17_a2
        # Reference
        if len(self.refs) > 0:
            self.refs += ','
        self.refs += 'Richter+17'


    def load_hotgas(self):
        """ Load data on hot gas (e.g. OVII, OVIII)
        Fang+15
        """

        # Init
        llist = LineList('EUV')
        ovii = AbsLine('OVII 21', linelist=llist)
        scoord = self.abs.scoord  # Sightline coordiantes

        # Fang+15  Table 1  [OVII]
        fang15_file = resource_filename('pyigm','/data/CGM/Galaxy/fang15_table1.dat')
        self.fang15 = Table.read(fang15_file, format='cds')
        print('Loading Fang+15 for OVII')
        # Reference
        if len(self.refs) > 0:
            self.refs += ','
        self.refs += 'Fang+15'
        # Generate the systems
        # # (should check to see if low-ion ones exist already)
        for row in self.fang15:
            # Coordinates
            gc = SkyCoord(l=row['GLON']*u.degree, b=row['GLAT']*u.degree, frame='galactic')
            # Limits
            # OVII line
            aline = ovii.copy()
            aline.attrib['coord'] = gc
            z = row['Vel']/c_kms
            try:
                aline.setz(z)
            except IOError:
                z = 0.
                vlim = np.array([-300,300]) * u.km/u.s
                aline.attrib['flag_EW'] = 3
                aline.attrib['EW'] = row['EW1'] / 1e3 * u.AA
                aline.attrib['sig_EW'] = 99. * u.AA
                #
                aline.attrib['flag_N'] = 0  # Might be able to set an upper limit
            else:
                aline.attrib['b'] = row['b'] * u.km / u.s
                aline.attrib['flag_EW'] = 1
                aline.attrib['EW'] = row['EW1'] / 1e3 * u.AA
                aline.attrib['sig_EW'] = row['e_EW1'] / 1e3 * u.AA
                vlim = np.array([-1,1]) * (2 * row['b'] + 2 * row['E_b']) * u.km/u.s
                # N_OVII
                aline.attrib['flag_N'] = 1
                aline.attrib['logN'] = row['logNO']
                aline.attrib['sig_logN'] = np.array([row['e_logNO'], row['E_logNO']])
                # Fill linear
                _,_ = linear_clm(aline.attrib)
            # OVII
            aline.limits.set(vlim)
            # Generate component and add
            comp = AbsComponent.from_abslines([aline])
            if aline.attrib['flag_N'] == 0: # Hack to merge later
                comp.attrib['sig_logN'] = np.array([0., 0.])
            else:
                pdb.set_trace()
            # Check for existing system
            minsep = np.min(comp.coord.separation(scoord).to('arcsec'))
            if minsep < 30*u.arcsec:  # Add component to existing system
                idx = np.argmin(comp.coord.separation(scoord).to('arcsec'))
                if self.verbose:
                    print("Adding OVII system to {}".format(self.abs.cgm_abs[idx].igm_sys))
                self.abs.cgm_abs[idx].igm_sys.add_component(comp, chk_sep=False, debug=True)
            else: # Instantiate
                abssys = IGMSystem(gc, z, vlim, name=row['Name']+'_z0', zem=row['z'])
                abssys.add_component(comp, chk_sep=False)
                # CGM Abs
                rho, ang_sep = calc_Galactic_rho(abssys.coord)
                cgmabs = CGMAbsSys(self.galaxy, abssys, rho=rho, ang_sep=ang_sep, cosmo=self.cosmo)
                # Add to cgm_abs
                self.abs.cgm_abs.append(cgmabs)

        scoord = self.abs.scoord  # Sightline coordiantes
        # Savage+03  Table 2  [OVI] -- Thick disk/halo only??
        print('Loading Savage+03 for OVI')
        savage03_file = resource_filename('pyigm', '/data/CGM/Galaxy/savage03_table2.fits')
        self.savage03 = Table.read(savage03_file)
        # Reference
        if len(self.refs) > 0:
            self.refs += ','
        self.refs += 'Savage+03'
        # Generate the systems
        # # (should check to see if low-ion ones exist already)
        for row in self.savage03:
            # Coordinates
            coord = SkyCoord(ra=row['_RA']*u.deg, dec=row['_DE']*u.deg, frame='icrs')
            gc = coord.transform_to('galactic')
            # Build the component
            vlim = np.array([row['V-'],row['V_']])*u.km/u.s
            comp = AbsComponent(gc, (8,6), 0., vlim)
            # Add attributes
            if row['b'] > 0.:
                comp.attrib['vcen'] = row['__v_obs']*u.km/u.s
                comp.attrib['sig_vcen'] = row['e__v_obs']*u.km/u.s
                comp.attrib['b'] = row['b']*u.km/u.s
                comp.attrib['sig_b'] = row['e_b']*u.km/u.s
                # Column
                comp.flag_N = 1
                comp.logN = row['logN_OVI_']
                comp.sig_logN = np.sqrt(row['e_sc']**2 + row['e_sys']**2)
            else: # Upper limit
                comp.flag_N = 3
                comp.logN = row['logN_OVI_']
                comp.sig_logN = 99.
            # Check for existing system
            minsep = np.min(comp.coord.separation(scoord).to('arcsec'))
            if minsep < 30*u.arcsec:
                idx = np.argmin(comp.coord.separation(scoord).to('arcsec'))
                self.abs.cgm_abs[idx].igm_sys.add_component(comp, chk_sep=False, debug=True,
                                                            update_vlim=True)
            else:  # New
                if row['RV'] > 0:
                    zem = row['RV']/c_kms
                else:
                    zem = row['z']
                abssys = IGMSystem(gc, comp.zcomp, vlim, name=row['Name']+'_z0', zem=zem)
                abssys.add_component(comp, chk_sep=False, debug=True)
                # CGM Abs
                rho, ang_sep = calc_Galactic_rho(abssys.coord)
                cgmabs = CGMAbsSys(self.galaxy, abssys, rho=rho, ang_sep=ang_sep, cosmo=self.cosmo)
                # Add to cgm_abs
                self.abs.cgm_abs.append(cgmabs)

    def write(self, outfil='MW_CGM.tar.gz'):
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
        return  # NOT YET IMPLEMENTED
        '''
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
        '''


