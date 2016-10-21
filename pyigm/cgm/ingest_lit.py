""" Module to Ingest CGM samples
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent

import pyigm
from pyigm.field.galaxy import Galaxy
from pyigm.abssys.igmsys import IGMSystem
from pyigm.cgm.cgm import CGMAbsSys
from pyigm.cgm.cgmsurvey import CGMAbsSurvey

def p11():
    """ Ingest Prochaska et al. 2011 CGM survey
    """
    # Low z OVI summary
    ovi_file = pyigm.__path__[0]+'/data/CGM/P11/lowovidat.fits'
    ovidat = Table.read(ovi_file)
    qso_radec = SkyCoord(ra=ovidat['QSO_RA'], dec=ovidat['QSO_DEC'], unit=(u.hourangle, u.deg))
    qso_nms = np.array([row['QSO'].strip() for row in ovidat])

    # CGM Survey
    p11 = CGMAbsSurvey(survey='P11', ref='Prochaska+11')

    # Dwarfs
    cgm_dwarf_file = pyigm.__path__[0]+'/data/CGM/P11/dwarf_galabs_strct.fits'
    cgm_dwarfs = Table.read(cgm_dwarf_file)
    # sub L*
    cgm_subls_file = pyigm.__path__[0]+'/data/CGM/P11/subls_galabs_strct.fits'
    cgm_subls = Table.read(cgm_subls_file)
    # L*
    cgm_lstar_file = pyigm.__path__[0]+'/data/CGM/P11/lstar_galabs_strct.fits'
    cgm_lstar = Table.read(cgm_lstar_file)

    # Loop on subsets
    for subset in [cgm_dwarfs, cgm_subls, cgm_lstar]:
        for row in subset:
            # RA, DEC
            # Galaxy
            gal = Galaxy((row['RA'], row['DEC']), z=row['Z'])
            gal.Lstar = row['DDEC']
            gal.type = row['GAL_TYPE']
            # IGMSys
            mtqso = np.where(qso_nms == row['FIELD'].strip())[0]
            if len(mtqso) != 1:
                pdb.set_trace()
                raise ValueError("No Field match")
            igmsys = IGMSystem(qso_radec[mtqso[0]], row['Z'], (-400.,400.)*u.km/u.s)
            # HI
            if row['MAG'][2] > 0.:
                # Lya
                lya = AbsLine(1215.67*u.AA, z=float(row['MAG'][3]))
                lya.attrib['EW'] = row['MAG'][4]/1e3*u.AA
                if row['MAG'][5] >= 99.:
                    lya.attrib['flag_EW'] = 3
                else:
                    lya.attrib['flag_EW'] = 1
                lya.attrib['sig_EW'] = row['MAG'][5]/1e3*u.AA
                # Ref
                lya.attrib['Ref'] = int(row['MAG'][2])
                # HI component
                if row['MAG'][9] <= 0.:
                    flagN = 3
                elif row['MAG'][9] > 9.:
                    flagN = 2
                else:
                    flagN = 1
                HIcomp = AbsComponent(qso_radec[mtqso[0]],(1,1),float(row['MAG'][3]),
                                      (-400,400)*u.km/u.s,
                                      Ntup=(flagN, row['MAG'][8], row['MAG'][9]))
                HIcomp._abslines.append(lya)
                igmsys._components.append(HIcomp)
                # NHI
                igmsys.NHI = HIcomp.logN
                igmsys.flag_NHI = HIcomp.flag_N
                igmsys.sig_NHI = HIcomp.sig_N
            # OVI
            if row['MAGERR'][2] > 0.:
                # OVI 1031
                ovi1031 = None
                if row['MAGERR'][4] > 0.:
                    ovi1031 = AbsLine(1031.9261*u.AA, z=float(row['MAGERR'][3]))
                    if row['MAGERR'][5] >= 99.:
                        ovi1031.attrib['flag_EW'] = 3
                    else:
                        ovi1031.attrib['flag_EW'] = 1
                    ovi1031.attrib['EW'] = row['MAGERR'][4]/1e3*u.AA
                    ovi1031.attrib['sig_EW'] = row['MAGERR'][5]/1e3*u.AA
                # OVI component
                if row['MAGERR'][9] <= 0.:
                    flagN = 3
                elif row['MAGERR'][9] > 9.:
                    flagN = 2
                else:
                    flagN = 1
                OVIcomp = AbsComponent(qso_radec[mtqso[0]],(8,6),float(row['MAGERR'][3]),
                                      (-400,400)*u.km/u.s,
                                      Ntup=(flagN, row['MAGERR'][8], row['MAGERR'][9]))
                if ovi1031 is not None:
                    OVIcomp._abslines.append(ovi1031)
                # Ref
                OVIcomp.Ref = int(row['MAG'][2])
                igmsys._components.append(OVIcomp)
            # CGM
            cgmabs = CGMAbsSys(gal, igmsys, chk_lowz=False)
            p11.cgm_abs.append(cgmabs)
    # Write tarball
    out_file = pyigm.__path__[0]+'/data/CGM/P11/P11_sys.tar'
    p11.to_json_tarball(out_file)


def main(flg):

    if (flg % 2**1) >= 2**0:
        p11()  # Prochaska et al. 2011

# Command line execution
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:
        flg = 0
        flg += 2**0   # P11
    else:
        flg = sys.argv[1]

    main(flg)
