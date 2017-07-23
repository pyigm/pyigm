""" Subclasses for DLA AbsSystem and AbsSurvey
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import warnings
import pdb
import numpy as np

from astropy import units as u

from linetools.isgm import utils as ltiu

from pyigm.abssys.igmsys import IGMSystem
from pyigm.abssys import utils as igmau

from .utils import dict_to_ions

class DLASystem(IGMSystem):
    """
    Class for a DLA absorption system

    Parameters
    ----------
    radec : tuple or coordinate
        RA/Dec of the sightline or astropy.coordinate
    zabs : float
      Absorption redshift
    vlim : Quantity array (2)
      Velocity limits of the system
      Defaulted to +/- 500 km/s if None
    NHI : float, required despite being a keyword
      log10 of HI column density
      must be 20.3 or higher
    **kwargs : keywords
      passed to AbsSystem.__init__
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
        slf._datdict = datdict

        # QSO keys
        slf.qso = slf._datdict['QSO name']
        slf.zem = float(slf._datdict['QSO zem'])
        # Name
        slf.name = '{:s}_z{:0.3f}'.format(slf.qso,zabs)

        # Abund
        slf.flg_ZH = float(slf._datdict['flg_mtl'])
        slf.ZH = float(slf._datdict['[M/H]'])
        slf.sig_ZH = float(slf._datdict['sig([M/H])'])

        return slf

    '''
    @classmethod
    def from_dict(cls, idict):
        """ Generate a DLASystem from a dict

        Parameters
        ----------
        idict : dict
          Usually read from the hard-drive
        """
        kwargs = dict(zem=idict['zem'], sig_NHI=idict['sig_NHI'],
                      name=idict['Name'])
        slf = cls(SkyCoord(idict['RA'], idict['DEC'], unit='deg'),
                  idict['zabs'], idict['vlim']*u.km/u.s, idict['NHI'],
                  **kwargs)
        # Components
        components = ltiu.build_components_from_dict(idict)
        for component in components:
            # This is to insure the components follow the rules
            slf.add_component(component)
        # Return
        return slf
    '''

    def __init__(self, radec, zabs, vlim, NHI, **kwargs):
        """Standard init

        NHI keyword is required
        """
        # NHI
        if NHI < 20.3:
            raise ValueError("This is not a DLA!  Try an LLS (or SLLS)")
        # vlim
        if vlim is None:
            vlim = [-500., 500.]*u.km/u.s
        # Generate with type
        IGMSystem.__init__(self, radec, zabs, vlim, NHI=NHI, abs_type='DLA', **kwargs)

    def model_abs(self, spec, **kwargs):
        """ Generate a model of the absorption from the DLA on an input spectrum
        This is a simple wrapper to pyigm.abssys.utils.hi_model

        Parameters
        ----------
        spec : XSpectrum1D

        Returns
        -------
        dla_model : XSpectrum1D
          Model spectrum with same wavelength as input spectrum
          Assumes a normalized flux
        lyman_lines : list
          List of AbsLine's that contributed to the DLA model

        """
        from pyigm.abssys.utils import hi_model
        vmodel, lines = hi_model(self, spec, **kwargs)
        # Return
        return vmodel, lines


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
          NOTE: This ignores velocity constraints on components (i.e. skip_vel=True)
        update_zvlim : bool, optional
          Update zvlim from lines in .clm (as applicable)
        linelist : LineList
        """
        if use_Nfile:
            clm_fil = self.tree+self._datdict['Abund file']
            # Read
            self._clmdict = igmau.read_clmfile(clm_fil, linelist=linelist)
            #pdb.set_trace()
            # Build components
            components = ltiu.build_components_from_dict(self._clmdict,
                                                         coord=self.coord,
                                                         chk_vel=False)
            # Read .ion file and fill in components
            ion_fil = self.tree+self._clmdict['ion_fil']
            self._indiv_ionclms = igmau.read_ion_file(ion_fil, components)
            # Parse .all file
            all_file = ion_fil.split('.ion')[0]+'.all'
            self.all_file=all_file  #MF: useful
            _ = igmau.read_all_file(all_file, components=components)
            # Build table
            self._ionN = ltiu.iontable_from_components(components, ztbl=self.zabs)
            # Add to AbsSystem
            for comp in components:
                self.add_component(comp)
        elif idict is not None:
            table = dict_to_ions(idict)
            self._ionN = table
        else:
            raise IOError("Not ready for this")

    def load_components(self, inp):
        """ Load components from an input object

        Parameters
        ----------
        inp : dict or ??
          Input object for loading the components
        """
        if isinstance(inp, dict):
            components = ltiu.build_components_from_dict(inp, coord=self.coord,
                                                         skip_vel=True)
            # Add in
            for component in components:
                self.add_component(component)
        else:
            raise NotImplementedError("Not ready for this input")

    # Output
    def __repr__(self):
        return ('<{:s}: {:s} {:s}, {:g}, NHI={:g}, Z/H={:g}>'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                 self.coord.dec.to_string(sep=':', pad=True),
                 self.zabs, self.NHI, self.ZH))

