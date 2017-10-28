""" Module to Ingest CGM samples
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb

from pkg_resources import resource_filename

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent
from linetools.lists.linelist import LineList
from linetools import utils as ltu

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


def ingest_burchett16():
    """ Ingest Burchett+16
    """
    # Virial matching
    b16_vir_file = resource_filename('pyigm', 'data/CGM/z0/Burchett16_HI_virselect_sfr.fits')
    b16_vir = Table.read(b16_vir_file)

    # CGM Survey
    b16 = CGMAbsSurvey(survey='B16', ref='Burchett+16')

    # Linelist
    llist = LineList('ISM')

    for row in b16_vir:
        # RA, DEC
        # Galaxy
        gal = Galaxy((row['ra_gal'], row['dec_gal']), z=row['zgal'])
        gal.SFR = row['SFR']
        gal.sig_SFR = row['SFR_err']
        gal.Mstar = row['mstars']
        gal.field = row['field']
        gal.RRvir = row['rrvir']
        gal.NSAidx = row['NSAidx']
        #
        igmsys = IGMSystem((row['ra_qso'], row['dec_qso']), row['zgal'], (-400., 400.) * u.km / u.s)
        # HI
        # Lya
        lya = AbsLine(1215.67 * u.AA, z=row['zgal'], linelist=llist)
        lya.attrib['EW'] = row['EW'] / 1e3 * u.AA
        if row['h1colsig'] <= 0.:
            lya.attrib['flag_EW'] = 3
        else:
            lya.attrib['flag_EW'] = 1
        lya.attrib['sig_EW'] = row['sigEW']
        # Ref
        lya.attrib['Ref'] = 'Burchett+16'
        # HI component
        if row['h1colsig'] >= 99.:
            flagN = 2
        elif row['h1colsig'] <= 0.:
            flagN = 3
        else:
            flagN = 1
        HIcomp = AbsComponent((row['ra_qso'], row['dec_qso']),
                              (1, 1), row['zgal'],
                              (-400, 400) * u.km / u.s,
                              Ntup=(flagN, row['h1col'], row['h1colsig']))
        HIcomp._abslines.append(lya)
        igmsys._components.append(HIcomp)
        # NHI
        igmsys.NHI = HIcomp.logN
        igmsys.flag_NHI = HIcomp.flag_N
        igmsys.sig_NHI = HIcomp.sig_N
        # CGM
        cgmabs = CGMAbsSys(gal, igmsys, chk_lowz=False)
        b16.cgm_abs.append(cgmabs)
    # Write tarball
    out_file = resource_filename('pyigm', '/data/CGM/z0/B16_sys.tar')
    b16.to_json_tarball(out_file)


def ingest_johnson15():
    """ Ingest Johnson+15
    """
    # Dict for QSO coords
    qsos = {}
    qsos['1ES1028+511'] = ltu.radec_to_coord('J103118.52517+505335.8193')
    qsos['FBQS1010+3003'] = ltu.radec_to_coord((152.5029167,30.056111))
    qsos['HE0226-4110'] = ltu.radec_to_coord('J022815.252-405714.62')
    qsos['HS1102+3441'] = ltu.radec_to_coord('J110539.8189+342534.672')
    qsos['LBQS1435-0134'] = ltu.radec_to_coord((219.451183,-1.786328))
    qsos['PG0832+251'] = ltu.radec_to_coord('J083535.8048+245940.146')
    qsos['PG1522+101'] = ltu.radec_to_coord((231.1023075, 9.9749372))
    qsos['PKS0405-123'] = ltu.radec_to_coord('J040748.4376-121136.662')
    qsos['SBS1108+560'] = ltu.radec_to_coord((167.8841667,55.790556))
    qsos['SBS1122+594'] = ltu.radec_to_coord((171.4741250,59.172667))
    qsos['Ton236'] = ltu.radec_to_coord((232.1691746,28.424928))

    # Virial matching
    j15_file = resource_filename('pyigm', 'data/CGM/z0/johnson2015_table1.fits')
    j15_tbl = Table.read(j15_file)

    # Clip COS-Halos
    keep = j15_tbl['Survey'] != 'COS-Halos'
    j15_tbl = j15_tbl[keep]

    # CGM Survey
    j15 = CGMAbsSurvey(survey='J15', ref='Johnson+15')

    # Linelist
    llist = LineList('ISM')

    for row in j15_tbl:
        # RA, DEC
        # Galaxy
        gal = Galaxy((row['RAJ2000'], row['DEJ2000']), z=float(row['zgal']))
        gal.Class = row['Class']
        gal.Mstar = row['logM_']
        gal.field = row['Name']
        gal.Env = row['Env']
        gal.d_Rh = row['d_Rh']
        #
        igmsys = IGMSystem(qsos[row['Name']], float(row['zgal']), (-400., 400.) * u.km / u.s)
        # HI
        if np.isnan(row['logNHI']):
            pass
        else:
            # HI component
            if row['l_logNHI'] == '<':
                flagN = 3
                sigNHI = 99.
            elif np.isnan(row['e_logNHI']):
                flagN = 2
                sigNHI = 99.
            else:
                flagN = 1
                sigNHI = row['e_logNHI']
            HIcomp = AbsComponent(qsos[row['Name']], (1, 1), float(row['zgal']),
                                      (-400, 400)*u.km/u.s, Ntup=(flagN, row['logNHI'], sigNHI))
            igmsys._components.append(HIcomp)
            # NHI
            igmsys.NHI = HIcomp.logN
            igmsys.flag_NHI = HIcomp.flag_N
            igmsys.sig_NHI = HIcomp.sig_N
        # OVI
        if np.isnan(row['logNHOVI']):
            pass
        else:
            # OVI component
            if row['l_logNHOVI'] == '<':
                flagN = 3
                sigNHOVI = 99.
            elif np.isnan(row['e_logNHOVI']):
                flagN = 2
                sigNHOVI = 99.
            else:
                flagN = 1
                sigNHOVI = row['e_logNHOVI']
            OVIcomp = AbsComponent(qsos[row['Name']], (8, 6), float(row['zgal']),
                                      (-400, 400)*u.km/u.s, Ntup=(flagN, row['logNHOVI'], sigNHOVI))
            igmsys._components.append(OVIcomp)
        # CGM
        cgmabs = CGMAbsSys(gal, igmsys, chk_lowz=False)
        j15.cgm_abs.append(cgmabs)
    # Write tarball
    out_file = resource_filename('pyigm', '/data/CGM/z0/J15_sys.tar')
    j15.to_json_tarball(out_file)

def main(flg):

    if (flg % 2**1) >= 2**0:
        p11()  # Prochaska et al. 2011
    if flg & 2**1:
        ingest_burchett16()
    if flg & 2**2:
        ingest_johnson15()

# Command line execution
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:
        flg = 0
        #flg += 2**0   # P11
        #flg += 2**1   # Burchett+16
        flg += 2**2   # Johnson+15
    else:
        flg = sys.argv[1]

    main(flg)
