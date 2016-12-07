""" Utilities for IGMSystem and IGMSurvey
Best to keep these separate from the Class modules
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import warnings
import pdb
import numpy as np
from collections import OrderedDict

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table, Column

from linetools.abund import ions as ltai
from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.lists.linelist import LineList


def dict_to_ions(idict):
    """  Manipulate dict into an ion astropy Table

    Will likely be deprecated

    Parameters
    ----------
    idict : dict

    Returns
    -------
    table : astropy.Table

    """
    #  Could probably use add_row or dict instantiation
    table = None
    for ion in idict.keys():
        Zion = ltai.name_ion(ion)
        if table is None:
            tkeys = idict[ion].keys()
            lst = [[idict[ion][tkey]] for tkey in tkeys]
            table = Table(lst, names=tkeys)
            # Extra columns
            if 'Z' not in tkeys:
                table.add_column(Column([Zion[0]], name='Z'))
                table.add_column(Column([Zion[1]], name='ion'))
        else:
            tdict = idict[ion]
            tkeys = idict[ion].keys()
            if 'Z' not in tkeys:
                tdict['Z'] = Zion[0]
                tdict['ion'] = Zion[1]
            # Add
            table.add_row(tdict)
    # Finish
    try:  # Historical keys
        table.rename_column('clm', 'logN')
    except:
        pass
    else:
        table.rename_column('sig_clm', 'sig_logN')
        table.rename_column('flg_clm', 'flag_N')
    # Return
    return table


def hi_model(abssys, spec, lya_only=False, add_lls=False, ret_tau=False,
             ignore_abslines=False, bval=30*u.km/u.s, **kwargs):
    """ Generate a model of the absorption from the absorption system
    on an input spectrum.

    Parameters
    ----------
    abssys : AbsSystem or list
      If list, must be a list of AbsSystem's
    spec : XSpectrum1D
    lya_only : bool, optional
      Only generate Lya
    ignore_abslines : bool, optional
      Ignore any existing abslines in the object
      NHI tag must be set
    add_lls : bool, optional
      Add Lyman continuum absorption
    bval : Quantity, optional
      Doppler parameter to use if abslines not adopted
    ret_tau : bool, optional
      Return only the optical depth (used for multiple systems)
    kwargs :
      Passed to voigt_from_abslines

    Returns
    -------
    vmodel : XSpectrum1D or ndarray
      Model spectrum with same wavelength as input spectrum
      Assumes a normalized flux
      Or optical depth array [used for multiple systems]
    lyman_lines : list
      List of AbsLine's that contributed to the model

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    from linetools.spectralline import AbsLine
    from linetools.analysis.voigt import voigt_from_abslines
    from linetools.analysis.absline import photo_cross
    # Input
    if isinstance(abssys, list):
        tau = None
        all_lines = []
        for iabssys in abssys:
            itau, ly_lines = hi_model(iabssys, spec, lya_only=lya_only, add_lls=add_lls,
                     ignore_abslines=ignore_abslines, bval=bval, ret_tau=True, **kwargs)
            all_lines += ly_lines
            if tau is None:
                tau = itau
            else:
                tau += itau
        # Flux
        flux = np.exp(-1*tau)
        vmodel = XSpectrum1D.from_tuple((spec.wavelength, flux))
        return vmodel, all_lines
    else:
        # Scan abs lines
        if not ignore_abslines:
            alines = []
        else:
            alines = abssys.list_of_abslines()
        lyman_lines = []
        lya_lines = []
        logNHIs = []
        # Scan alines
        for aline in alines:
            # Lya
            if aline.name == 'HI 1215':
                lya_lines.append(aline)
                logNHIs.append(np.log10(aline.attrib['N'].value))
            # Any HI
            if 'HI' in aline.name:
                lyman_lines.append(aline)
        if len(lya_lines) > 0: # Use the lines
            # Check we have a DLA worth
            if np.log10(np.sum(10**np.array(logNHIs))) < abssys.NHI:
                raise ValueError("Total NHI of the Lya lines is less than NHI of the system!  Something is wrong..")
        else: # Generate one
            warnings.warn("Generating the absorption lines from the system info, not abslines")
            if lya_only:
                lya_line = AbsLine('HI 1215', z=abssys.zabs)
                lya_line.attrib['N'] = 10**abssys.NHI / u.cm**2
                lya_line.attrib['b'] = bval
                lyman_lines.append(lya_line)
            else:
                HIlines = LineList('HI')
                wrest = HIlines._data['wrest']
                for iwrest in wrest:
                    # On the spectrum?
                    if iwrest >= spec.wvmin/(1+abssys.zabs):
                        lyman_line = AbsLine(iwrest, linelist=HIlines, z=abssys.zabs)
                        lyman_line.attrib['N'] = 10**abssys.NHI / u.cm**2
                        lyman_line.attrib['b'] = bval
                        lyman_lines.append(lyman_line)
        # tau for abs lines
        if len(lyman_lines) == 0:
            pdb.set_trace()
        tau_Lyman = voigt_from_abslines(spec.wavelength,lyman_lines, ret='tau', **kwargs)
        # LLS?
        if add_lls:
            wv_rest = spec.wavelength / (1+abssys.zabs)
            energy = wv_rest.to(u.eV, equivalencies=u.spectral())
            # Get photo_cross and calculate tau
            tau_LL = (10.**abssys.NHI / u.cm**2) * photo_cross(1,1,energy)
            # Kludge
            pix_LL = np.argmin(np.fabs(wv_rest- 911.3*u.AA))
            pix_kludge = np.where((wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA))[0]
            tau_LL[pix_kludge] = tau_LL[pix_LL]
            # Generate the spectrum
        else:
            tau_LL = 0.
        # Sum
        tau_tot = tau_LL + tau_Lyman
        if ret_tau:
            vmodel = tau_tot
        else:
            flux = np.exp(-1*tau_tot)
            vmodel = XSpectrum1D.from_tuple((spec.wavelength, flux))
        # Return
        return vmodel, lyman_lines

def parse_datdict(datdict):
    """ Parse a datdict

    Parameters
    ----------
    datdict : OrderedDict
      dict from .dat file

    Returns
    -------
    tuple of coord, zabs, name, NHI, sigNHI, clmfil
    """
    # RA/DEC
    try:
        ras, decs = (datdict['RA (2000)'], datdict['DEC (2000)'])
    except KeyError:
        ras, decs = ('00 00 00', '+00 00 00')
    coord = SkyCoord(ras, decs, frame='icrs', unit=(u.hour, u.deg))

    # zabs
    try:
        zabs = float(datdict['zabs'])
    except KeyError:
        zabs = 0.

    # Name
    name = ('J' +
                 coord.ra.to_string(unit=u.hour, sep='', pad=True) +
                 coord.dec.to_string(sep='', pad=True, alwayssign=True) +
                 '_z{:0.3f}'.format(zabs))

    # NHI
    try:
        NHI = float(datdict['NHI'])  # DLA format
    except KeyError:
        try:
            NHI = float(datdict['NHI tot'])  # LLS format
        except KeyError:
            NHI = 0.

    # NHIsig
    try:
        key_sigNHI = datdict['sig(NHI)']  # DLA format
    except KeyError:
        try:
            key_sigNHI = datdict['NHI sig']  # LLS format
        except KeyError:
            key_sigNHI = '0.0 0.0'
    sigNHI = np.array(map(float,key_sigNHI.split()))

    # Abund file
    try:
        key_clmfil = datdict['Abund file']  # DLA format
    except KeyError:
        key_clmfil = ''
    clm_fil = key_clmfil.strip()

    # Return
    return coord, zabs, name, NHI, sigNHI, clm_fil


def read_all_file(all_file, components=None, verbose=False):
    """Read in JXP-style .all file in an appropriate manner

    NOTE: If program breaks in this function, check the all file
    to see if it is properly formatted.

    Fills components if inputted
    May need to worry more about CII*

    Parameters
    ----------
    all_file : str
      Full path to the .all file
    components : list, optional
      List of AbsComponent objects
    verbose : bool, optional
    """
    # Read
    if verbose:
        print('Reading {:s}'.format(all_file))
    names = ('Z', 'ion', 'logN', 'sig_logN', 'flag_N', 'flg_inst')  # was using flg_clm
    table = ascii.read(all_file, format='no_header', names=names)

    # Fill components
    if components is not None:
        allZ = np.array([comp.Zion[0] for comp in components])
        allion = np.array([comp.Zion[1] for comp in components])
        allEj = np.array([comp.Ej.value for comp in components])
        # Loop
        for row in table:
            mt = np.where((allZ == row['Z']) & (allion == row['ion']) & (allEj == 0.))[0]
            if len(mt) == 0:
                pass
            elif len(mt) == 1:
                # Fill
                components[mt[0]].flag_N = row['flag_N']
                components[mt[0]].logN = row['logN']
                components[mt[0]].sig_logN = row['sig_logN']
            else:
                raise ValueError("Found multiple component matches in read_all_file")
    # Write
    return table


def read_clmfile(clm_file, linelist=None, verbose=True):
    """ Read in a .CLM file in an appropriate manner

    NOTE: If program breaks in this function, check the clm to see if it is properly formatted.


    RETURNS two dictionaries CLM and LINEDIC. CLM contains the contents of CLM
    for the given DLA. THe LINEDIC that is passed (when not None) is updated appropriately.

    Keys in the CLM dictionary are:
      INST - Instrument used
      FITS - a list of fits files
      ZABS - absorption redshift
      ION - .ION file location
      HI - THe HI column and error; [HI, HIerr]
      FIX - Any abundances that need fixing from the ION file
      VELS - Dictioanry of velocity limits, which is keyed by
        FLAGS - Any measurment flags assosicated with VLIM
        VLIM - velocity limits in km/s [vmin,vmax]
        ELEM - ELement (from get_elem)

    See get_elem for properties of LINEDIC

    Parameters
    ----------
    clm_file : str
      Full path to the .clm file
    linelist : LineList
      can speed up performance
    """
    clm_dict = {}
    # Read file
    f=open(clm_file, 'r')
    arr=f.readlines()
    f.close()
    nline = len(arr)
    # Data files
    clm_dict['flg_data'] = int(arr[1][:-1])
    clm_dict['fits_files']={}
    ii = 2
    for jj in range(0,6):
        if (clm_dict['flg_data'] % (2**(jj+1))) > (2**jj - 1):
            clm_dict['fits_files'][2**jj] = arr[ii].strip()
            ii += 1

    # Redshift
    clm_dict['zsys']=float(arr[ii][:-1]) ; ii += 1
    clm_dict['ion_fil']=arr[ii].strip() ; ii += 1
    # NHI
    tmp = arr[ii].split(','); ii += 1
    if len(tmp) != 2:
        raise ValueError('ionic_clm: Bad formatting {:s} in {:s}'.format(arr[ii-1], clm_file))
    clm_dict['NHI']=float(tmp[0])
    clm_dict['sigNHI']=float(tmp[1])
    # Abundances by hand
    numhand=int(arr[ii][:-1]); ii += 1
    clm_dict['fixabund'] = {}
    if numhand > 0:
        for jj in range(numhand):
            # Atomic number
            atom = int(arr[ii][:-1]) ; ii += 1
            # Values
            tmp = arr[ii].strip().split(',') ; ii += 1
            clm_dict['fixabund'][atom] = float(tmp[0]), float(tmp[1]), int(tmp[2])
    # Loop on lines
    clm_dict['lines'] = {}
    while ii < (nline-1):
        # No empty lines allowed
        if len(arr[ii].strip()) == 0:
           break
        # Read flag
        ionflg = int(arr[ii].strip()); ii += 1
        # Read the rest
        tmp = arr[ii].split(',') ; ii += 1
        if len(tmp) != 4:
            raise ValueError('ionic_clm: Bad formatting {:s} in {:s}'.format(arr[ii-1],clm_file))
        vmin = float(tmp[1].strip())
        vmax = float(tmp[2].strip())
        key = float(tmp[0].strip()) # Using a numpy float not string!
        # Generate
        clm_dict['lines'][key] = AbsLine(key*u.AA,closest=True,linelist=linelist, verbose=verbose)
        clm_dict['lines'][key].setz(clm_dict['zsys'])
        clm_dict['lines'][key].analy['FLAGS'] = ionflg, int(tmp[3].strip())
        # By-hand
        if ionflg >= 8:
            clm_dict['lines'][key].attrib['N'] = 10.**vmin / u.cm**2
            clm_dict['lines'][key].attrib['sig_N'] = (10.**(vmin+vmax) - 10.**(vmin-vmax))/2/u.cm**2
        else:
            clm_dict['lines'][key].analy['vlim']= [vmin,vmax]*u.km/u.s
    # Return
    return clm_dict


def read_dat_file(dat_file):
    """ Read an ASCII ".dat" file from JXP format 'database'

    Parameters
    ----------
    dat_file : str
     filename

    Returns
    -------
    dat_dict : OrderedDict
      A dict containing the info in the .dat file
    """
    # Define
    datdict = OrderedDict()
    # Open
    f = open(dat_file, 'r')
    for line in f:
        tmp = line.split('! ')
        tkey = tmp[1].strip()
        key = tkey
        val = tmp[0].strip()
        datdict[key] = val
    f.close()

    return datdict


def read_ion_file(ion_fil, components, lines=None, linelist=None, tol=0.05*u.AA):
    """ Read in JXP-style .ion file in an appropriate manner

    NOTE: If program breaks in this function, check the .ion file
    to see if it is properly formatted.

    If components is passed in, these are filled as applicable.

    Parameters
    ----------
    ion_fil : str
      Full path to .ion file
    components : list
      List of AbsComponent objects
    lines : list, optional
      List of AbsLine objects [used for historical reasons, mainly]
    linelist : LineList
      May speed up performance
    tol : Quantity, optional
      Tolerance for matching wrest
    """
    # Read
    names = ('wrest', 'logN', 'sig_logN', 'flag_N', 'flg_inst')
    table = ascii.read(ion_fil, format='no_header', names=names)

    # Generate look-up table for quick searching
    all_wv = []
    all_idx = []
    for jj,comp in enumerate(components):
        for kk,iline in enumerate(comp._abslines):
            all_wv.append(iline.wrest)
            all_idx.append((jj,kk))
    all_wv = u.Quantity(all_wv)
    # Loop now
    for row in table:
        mt = np.where(np.abs(all_wv-row['wrest']*u.AA)<tol)[0]
        if len(mt) == 0:
            pass
        elif len(mt) == 1:
            # Fill
            jj = all_idx[mt[0]][0]
            kk = all_idx[mt[0]][1]
            components[jj]._abslines[kk].attrib['flag_N'] = row['flag_N']
            components[jj]._abslines[kk].attrib['logN'] = row['logN']
            components[jj]._abslines[kk].attrib['sig_logN'] = row['sig_logN']
            components[jj]._abslines[kk].analy['flg_inst'] = row['flg_inst']
        else:
            pdb.set_trace()
            raise ValueError("Matched multiple lines in read_ion_file")
    # Return
    return table


def sum_ionN(tbl1, tbl2):
    """ Sum two ion column density tables

    Parameters
    ----------
    tbl1 : Table
    tbl2 : Table

    Returns
    -------
    sum_tbl : Table

    """
    # Instantiate and use data form original as starting point
    sum_tbl = tbl1.copy()

    # Loop through other
    for row2 in tbl2:
        # New?
        Zion = (row2['Z'], row2['ion'])
        try:
            row1 = sum_tbl[(sum_tbl['Z'] == Zion[0])&(sum_tbl['ion'] == Zion[1])]
        except KeyError:
            # Add in the new row
            sum_tbl.add_row(row2)
        else:
            idx = np.where((sum_tbl['Z'] == Zion[0]) &
                           (sum_tbl['ion'] == Zion[1]))[0][0]
            # Clm
            flagN, logN, siglogN = ltaa.sum_logN(row1, row2)
            sum_tbl['logN'][idx] = logN
            # Error
            sum_tbl['sig_logN'][idx] = siglogN
            # Flag
            flags = [row1['flag_N'], row2['flag_N']]
            if 2 in flags:   # At least one saturated
                flag = 2
            elif 1 in flags:  # None saturated; at least one detection
                flag = 1
            else:            # Both upper limits
                flag = 3
            sum_tbl['flag_N'][idx] = flag
            # Instrument (assuming binary flag)
            if 'flg_inst' in row1.keys():
                binflg = [0]*10
                for jj in range(10):
                    if (row2['flg_inst'] % 2**(jj+1)) >= 2**jj:
                        binflg[jj] = 1
                    if (row1['flg_inst'] % 2**(jj+1)) >= 2**jj:
                        binflg[jj] = 1
                sum_tbl['flg_inst'][idx] = int(np.sum(
                    [2**kk for kk, ibinf in enumerate(binflg) if ibinf==1]))
    # Return
    return sum_tbl

def class_by_type(type):
    """ Enable IGMSystem init by type

    Parameters
    ----------
    type : str

    Returns
    -------

    """
    from .lls import LLSSystem
    from .dla import DLASystem

    if type == 'LLS':
        system = LLSSystem
    elif type == 'DLA':
        system = DLASystem
    else:
        raise IOError("Bad system type!")
    # Return
    return system
