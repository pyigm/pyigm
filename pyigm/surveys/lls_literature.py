"""
#;+ 
#; NAME:
#; lls_literature
#;   Ordered by publication date
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for loading up literature data on Lyman Limit Systems
#;   29-Jun-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import imp, glob
import numpy as np
import urllib2
import pdb

from astropy import units as u

from linetools.lists.linelist import LineList
from linetools.spectralline import AbsLine
from linetools.abund import ions as ltai
from linetools.analysis import absline as ltaa

from pyigm.abssys.igmsys import AbsSubSystem
from pyigm.abssys import utils as pyiau
from pyigm.abssys.lls import LLSSystem
from .llssurvey import LLSSurvey

pyigm_path = imp.find_module('pyigm')[1]


def zonak2004():
    """Zonak, S. et al. 2004, ApJ, 2004, 606, 196

    PG1634+706
    HST+Keck spectra
    MgII, SiIV, SiIII from Table 2.  Summing Subsystems A (Model 2) and B
       Errors estimated by JXP (not reported)
       SiIII in A may be a model
       SiIV in B may be a model
    Total NHI from LL. Taken from Fig 3 caption.  
       Error estimated by JXP 
    Not all EWs in Table 1 included
    Adopting their M/H
    """
    # Setup
    radec = '163428.9897+703132.422'  # SIMBAD
    lls = LLSSystem(name='PG1634+706_z1.041', radec=radec, zem=1.337,
        zabs=1.0414, vlim=[-250., 100.]*u.km/u.s, NHI=17.23, ZH=-1.4,
        sig_NHI=np.array([0.15,0.15]))
    # SubSystems
    lls.nsub = 2
    # Abundances
    adict = dict(MgII={'clm': log_sum([11.45,11.90,12.02,11.68]), 'sig_clm': 0.05, 'flg_clm': 1},
        SiIII={'clm': log_sum([12.5,12.5,12.8,12.7]), 'sig_clm': 0.25, 'flg_clm': 1},
        SiIV={'clm': log_sum([10.9,10.8,11.2,11.1]), 'sig_clm': 0.15, 'flg_clm': 1} )
    lls.subsys['A'] = AbsSubSystem(lls, 1.0414, [-80, 100]*u.km/u.s, 'A')
    lls.subsys['A']._ionN = pyiau.dict_to_ions(adict)

    bdict = dict(SiIII={'clm': log_sum([11.8,12.8,12.4]), 'sig_clm': 0.15, 'flg_clm': 1},
        SiIV={'clm': log_sum([11.2,12.2,11.8]), 'sig_clm': 0.15, 'flg_clm': 1} )
    lls.subsys['B'] = AbsSubSystem(lls, 1.0414, [-240, -80]*u.km/u.s, 'B')
    lls.subsys['B']._ionN = pyiau.dict_to_ions(bdict)
    # Total
    lls._ionN = pyiau.sum_ionN(lls.subsys['A']._ionN, lls.subsys['B']._ionN)
    lls.Refs.append('Zon04')
    # Return
    return lls

def jenkins2005():
    """Jenkins, E. et al. 2005, ApJ, 2005, 623, 767
    PHL 1811
    HST/STIS, FUSE
    Metals parsed from Table 1
      OI taken from text
      Had to input error on columns by hand (JXP)
    Total NHI from Lyman series. see Fig 3
    M/H from O/H
    """
    # Grab ASCII file from ApJ
    tab_fil = pyigm_path+"/data/LLS/Literature/jenkins2005.tb1.ascii"
    chk_fil = glob.glob(tab_fil)
    if len(chk_fil) > 0:
        tab_fil = chk_fil[0]
    else:
        url = 'http://iopscience.iop.org/0004-637X/623/2/767/fulltext/61520.tb1.txt'
        print('LLSSurvey: Grabbing table file from {:s}'.format(url))
        f = urllib2.urlopen(url)
        with open(tab_fil, "wb") as code:
            code.write(f.read())
    # Setup
    radec = '215501.5152-092224.688'  # SIMBAD
    lls = LLSSystem(name='PHL1811_z0.081', radec=radec, zem=0.192,
        zabs=0.080923, vlim=[-100., 100.]*u.km/u.s, NHI=17.98, ZH=-0.19,
        sig_NHI=np.array([0.05,0.05]))
    lls.lines = []  # Probably not used

    # AbsLines
    ism = LineList('ISM')
    Nsig = {'C IV': 0.4, 'N II': 0.4, 'Si II': 0.05, 'Si IV': 0.25, 
        'S II': 0.2, 'Fe II': 0.12, 'H I': 0.05, 'S III': 0.06}

    # Parse Table
    with open(tab_fil,'r') as f:
        flines = f.readlines()
    ion_dict = {}
    for iline in flines:
        iline = iline.strip()
        if (len(iline) == 0): 
            continue
        # Split on tabs
        isplit = iline.split('\t')
        # Offset?
        ioff = 0
        if isplit[0][0] in ['1','2']:
            ioff = -1
        # Catch bad lines
        if (isplit[1+ioff][0:6] in ['1442.0','1443.7','1120.9']): # Skip goofy CII line and CII*
            continue
        if len(isplit[2+ioff]) == 0: 
            continue
        # Ion
        if (len(isplit[0].strip()) > 0) & (isplit[0][0] not in ['1','2']):
            ionc = isplit[0].strip()
            try:
                Zion = ltai.name_ion(ionc)
            except KeyError:
                pdb.set_trace()
        # Generate the Line
        try:
            newline = AbsLine(float(isplit[2+ioff])*u.AA,linelist=ism, closest=True)
        except ValueError:
            pdb.set_trace()
        newline.attrib['z'] = lls.zabs
        # Spectrum
        newline.analy['datafile'] = 'STIS' if 'S' in isplit[1] else 'FUSE'
        # EW
        try:
            EWvals = isplit[4+ioff].split(' ')
        except IndexError:
            pdb.set_trace()
        newline.attrib['EW'] = float(EWvals[0])*u.AA/1e3
        newline.attrib['sig_EW'] = float(EWvals[2])*u.AA/1e3
        newline.attrib['flag_EW'] = 1
        if len(isplit) < (5+ioff+1):
            continue
        # Colm?
        #xdb.set_trace()
        newline.attrib['sig_logN'] = 0.
        if (len(isplit[5+ioff].strip()) > 0) & (isplit[5+ioff].strip() != '\\ldots'):
            if isplit[5+ioff][0] == '\\':
                ipos = isplit[5+ioff].find(' ')
                newline.attrib['logN'] = float(isplit[5+ioff][ipos+1:])
                newline.attrib['flag_N'] = 2
            elif isplit[5+ioff][0] == '<':
                ipos = 0
                newline.attrib['logN'] = float(isplit[5+ioff][ipos+1:])
                newline.attrib['flag_N'] = 3
            elif isplit[5+ioff][0] == '1':
                try:
                    newline.attrib['logN'] = float(isplit[5+ioff][0:5])
                except ValueError:
                    pdb.set_trace()
                newline.attrib['flag_N'] = 1
                try:
                    newline.attrib['sig_logN'] = Nsig[ionc]
                except KeyError:
                    print('No error for {:s}'.format(ionc))
            else:
                raise ValueError('Bad character')
            # ion_dict
            ion_dict[ionc] = dict(clm=newline.attrib['logN'], sig_clm=newline.attrib['sig_logN'],
                flg_clm=newline.attrib['flag_N'], Z=Zion[0], ion=Zion[1])
        # Append
        lls.lines.append(newline)
    # Fix NI, OI
    ion_dict['O I']['clm'] = 14.47
    ion_dict['O I']['sig_clm'] = 0.05
    ion_dict['N I']['flg_clm'] = 3
    lls._ionN = pyiau.dict_to_ions(ion_dict)

    lls.Refs.append('Jen05')
    # Return
    return lls

def tripp2005():
    '''Tripp, T. et al. 2005, ApJ, 2005, 619, 714
    PG 1216+069 (LLS in Virgo)
    HST/STIS, FUSE
    Metal columns parsed from Tables 2 and 3
    Total NHI from damping wings
    M/H from O/H
    '''
    # Grab ASCII files from ApJ
    tab_fils = [pyigm_path+"/data/LLS/tripp2005.tb3.ascii",
                pyigm_path+"/data/LLS/tripp2005.tb2.ascii"]
    urls = ['http://iopscience.iop.org/0004-637X/619/2/714/fulltext/60797.tb3.txt',
        'http://iopscience.iop.org/0004-637X/619/2/714/fulltext/60797.tb2.txt']
    for jj,tab_fil in enumerate(tab_fils):
        chk_fil = glob.glob(tab_fil)
        if len(chk_fil) > 0:
            tab_fil = chk_fil[0]
        else:
            url = urls[jj]
            print('LLSSurvey: Grabbing table file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            with open(tab_fil, "wb") as code:
                code.write(f.read())
    # Setup
    radec = '121920.9320+063838.476'  # SIMBAD
    lls = LLSSystem(name='PG1216+069_z0.006', radec=radec, zem=0.3313,
        zabs=0.00632, vlim=[-100., 100.]*u.km/u.s, NHI=19.32, ZH=-1.6,
        sig_NHI=np.array([0.03, 0.03]))

    # Columns
    # Start with Table 3 (VPFIT)
    with open(tab_fils[0],'r') as f:
        flines3 = f.readlines()
    ion_dict = {}
    for iline in flines3:
        if (len(iline.strip()) == 0): 
            continue
        isplit = iline.split('\t')
        # Ion
        flg = 2
        if (len(isplit[0].strip()) > 0):# & (isplit[0][0] not in ['1','2']):
            ipos = isplit[0].find('1')
            ionc = isplit[0][0:ipos-1].strip()
            try:
                Zion = ltai.name_ion(ionc)
            except KeyError:
                pdb.set_trace()
            flg = 1
        # Column
        csplit = isplit[3].split(' ')
        clm = float(csplit[0])
        sig = float(csplit[2])
        if flg == 1:
            ion_dict[ionc] = dict(logN=clm, sig_logN=sig, flag_N=1,
                                  Z=Zion[0], ion=Zion[1])
        else:  # Add it in
            tmp_dict = dict(logN=clm, sig_logN=sig, flag_N=1, Z=Zion[0],
                            ion=Zion[1])
            flagN, logN, siglogN = ltaa.sum_logN(ion_dict[ionc], tmp_dict)
            ion_dict[ionc]['logN'] = logN
            ion_dict[ionc]['sig_logN'] = siglogN
    ions = ion_dict.keys()

    # Now Table 2 for the extras
    with open(tab_fils[1],'r') as f:
        flines2 = f.readlines()
    # Trim the first 10 lines
    flines2 = flines2[10:]
    # Loop
    for iline in flines2:
        isplit = iline.split('\t')
        #
        ionc = isplit[0].strip()
        if (len(ionc) == 0) or (ionc in ions):
            continue
        #
        Zion = ltai.name_ion(ionc)
        ion_dict[ionc] = dict(Z=Zion[0], ion=Zion[1], sig_logN=0.)
        if isplit[4][0] == '<':
            ion_dict[ionc]['logN'] = float(isplit[4][1:])
            ion_dict[ionc]['flag_N'] = 3
        else:
            raise ValueError('Should not get here')


    # Finish
    lls._ionN = pyiau.dict_to_ions(ion_dict)
    lls.Refs.append('Tri05')
    return lls

def peroux06a():
    """Peroux, C. et al. 2006a, MNRAS, 372, 369
    SDSS J0134+0051
    One of her sample
    Metal columns taken by JXP from Table 2 (no online data)
    Total NHI from damping wings
    """
    # Setup
    radec = '013405.75+005109.4'  # SDSS Name
    lls = LLSSystem(name='SDSSJ0134+0051_z0.842', radec=radec, zem=1.522,
        zabs=0.842, vlim=[-150., 150.]*u.km/u.s, NHI=19.93, sig_NHI=np.array([0.15,0.15]))
    #  Table 2
    ion_dict = {}
    N = np.sum(np.array([5.56,12.6,13.7,23.5,61.4,39.8,6,9.14])*1e10)
    sig = np.sqrt(np.sum((np.array([2.32,3.1,3.68,4.13,8.02,6.65,3.37,2.82])*1e10)**2))
    ion_dict['Mg I'] = dict(clm=np.log10(N), sig_clm=sig/N/np.log(10),flg_clm=1,Z=12,ion=1)
    ion_dict['Mg II'] = dict(clm=np.log10(5e13), sig_clm=0.,flg_clm=2,Z=12,ion=2)
    N = np.sum(np.array([8.17,4.28,32.1,125,710,301,893,600,263,65.7])*1e11)
    sig = np.sqrt(np.sum((np.array([2.63,1.40,2.37,8.6,53.2,28.4,73.5,61.7,14.0,2.95])*1e11)**2))
    ion_dict['Fe II'] = dict(clm=np.log10(N), sig_clm=sig/N/np.log(10),flg_clm=1,Z=26,ion=2)
    sig = np.sqrt(np.sum((np.array([3.72,1.84,2.36,3.83])*1e11)**2))
    ion_dict['Zn II'] = dict(clm=np.log10(2*sig), sig_clm=0.,flg_clm=3,Z=30,ion=2)
    sig = np.sqrt(np.sum((np.array([19.4,9.79])*1e11)**2))
    ion_dict['Cr II'] = dict(clm=np.log10(2*sig), sig_clm=0.,flg_clm=3,Z=24,ion=2)
    # Not including MnII.  Appears as a detection but also given as a limit..

    # Finish
    lls._ionN = pyiau.dict_to_ions(ion_dict)
    lls.Refs.append('Prx06a')
    return lls
 
def peroux06b():
    """Peroux, C. et al. 2006b, A&A, 450, 53
    SDSS J1323-0021
    Metal rich
    Metal columns copied by JXP from Table 1 
    Total NHI from damping wings
    """
    # Setup
    radec = '132323.78-002155.2'  # SDSS Name
    lls = LLSSystem(name='SDSSJ1323-0021_z0.716', radec=radec, zem=1.390,
        zabs=0.716, vlim=[-200., 200.]*u.km/u.s, NHI=20.21,
                    sig_NHI=np.array([0.20,0.20]))
    # Parse table file
    tab_fil = pyigm_path+"/data/LLS/Literature/peroux06b.tb1.ascii"
    with open(tab_fil,'r') as f:
        flines = f.readlines()
    ion_dict = {}
    for iline in flines:
        isplit = iline.split('\t')
        if len(isplit[0]) == 0:
            # Grab ions and init
            ions = isplit[3:10]
            for ion in ions:
                Zion = ltai.name_ion(ion)
                ion_dict[ion] = dict(clm=0., sig_clm=0.,flg_clm=1,Z=Zion[0],ion=Zion[1])
            continue
        # Column or sigma?
        if isplit[0][0] == 'N': # Column
            for kk,iis in enumerate(isplit[3:10]):
                ion = ions[kk]
                if iis[0] == '>':
                    ion_dict[ion]['flg_clm'] = 2
                    ion_dict[ion]['clm'] += float(iis[1:])
                elif iis[0] == '<':
                    pass
                elif iis[0] == '.':
                    pass
                else:
                    ion_dict[ion]['clm'] += float(iis)
        else: # Sigma
            for kk,iis in enumerate(isplit[3:10]):
                ion = ions[kk]
                if iis[0] == '.':
                    pass
                else:
                    ion_dict[ion]['sig_clm'] += float(iis)**2
    # Convert to log
    for ion in ions:
        N = ion_dict[ion]['clm']
        sig = np.sqrt(ion_dict[ion]['sig_clm'])
        #
        ion_dict[ion]['clm'] = np.log10(N)
        if ion_dict[ion]['flg_clm'] == 2:
            ion_dict[ion]['sig_clm'] = 0.
        else:
            ion_dict[ion]['sig_clm'] = sig/N/np.log(10)
    # Finish
    lls._ionN = pyiau.dict_to_ions(ion_dict)
    lls.Refs.append('Prx06b')
    return lls

def meiring06():
    """Meiring et al. 2006, MNRAS, 370, 43
    Q1107+0003 
    Taken from Table 4 by JXP
    NHI from RTN06 (damping wings)
    RA/DEC from STIS header
    """
    # Setup
    lls = LLSSystem(name='SDSSJ1107+0003_z0.954', radec=(166.90273*u.deg,
                                                         0.05795000*u.deg),
                    zem=1.726, zabs=0.9542,
                    vlim=[-300., 300.]*u.km/u.s, NHI=20.26,
                    sig_NHI=np.array([0.14,0.09]))
    #  Meiring06, Table 4
    ion_dict = {}
    ion_dict['Zn II'] = dict(clm=12.08, sig_clm=0.,flg_clm=3,Z=30,ion=2)
    ion_dict['Ti II'] = dict(clm=13.01, sig_clm=0.,flg_clm=3,Z=22,ion=2)
    ion_dict['Cr II'] = dict(clm=12.76, sig_clm=0.,flg_clm=3,Z=24,ion=2)
    # Finish
    lls._ionN = pyiau.dict_to_ions(ion_dict)
    lls.Refs.append('Mei06')
    return lls

def meiring07():
    """Meiring et al. 2007, MNRAS, 376, 557
    SLLS with Magellan
    Abundances from Table 11 from astro-ph (LateX) by JXP [AODM]
    RA/DEC from Table 1
    """
    all_lls = []
    # Table 1
    tab_fil = pyigm_path+"/data/LLS/Literature/meiring07.tb1.ascii"
    with open(tab_fil,'r') as f:
        flines1 = f.readlines()
    # Grab RA/DEC
    qso_dict = {}
    for iline in flines1:
        if iline[0:2] in ['QS','\h','$\\', 'J2']:
            continue
        # Parse
        isplit = iline.split('&')
        if '-' not in isplit[3]:
            sgn = '+'
        else:
            sgn = ''
        radec = isplit[2].strip()+sgn+isplit[3].strip()
        radec = radec.replace(':','')
        # zem
        if isplit[0].strip() != 'Q0826-2230':
            zem = float(isplit[5].strip())
        else:
            zem = 0.911
        # Save
        qso_dict[isplit[0].strip()] = dict(radec=radec, zem=zem,
                                           vlim=[-500.,500]*u.km/u.s)
    # Abundances (AODM)
    # Table 11
    tab_fil = pyigm_path+"/data/LLS/Literature/meiring07.tb11.ascii"
    with open(tab_fil,'r') as f:
        flines11 = f.readlines()
    #
    for iline in flines11:
        if iline[0:2] in ['\h','  ']:
            continue
        # Parse
        isplit = iline.split('&')
        # Ions
        if iline[0:2] == 'QS':
            ioncs = []
            Zions = []
            for iis in isplit[3:-1]: # Skipping HI
                # Parse
                is2 = iis.split('\\')
                ip2 = is2[2].find('}')
                ionc = is2[1][2:].strip()+' '+is2[2][0:ip2].strip()
                # Zion
                Zion = ltai.name_ion(ionc)
                # Append
                ioncs.append(ionc)
                Zions.append(Zion)
            continue
        if iline[0] == 'Q':
            # QSO
            qso = isplit[0].strip()
            # zabs and name
            zabs = float(isplit[1].strip())
            qso_dict[qso]['name']=qso+'z_{:.3f}'.format(zabs) 
            qso_dict[qso]['zabs']=zabs
            # NHI
            is2 = isplit[2].strip()
            qso_dict[qso]['NHI'] = float(is2[0:5])
            #if qso_dict[qso]['NHI'] >= 20.3:
            #    print('Uh oh.  DLA')
            qso_dict[qso]['sig_NHI'] = np.array([float(is2[10:])]*2)
            # Generate LLS
            lls = LLSSystem(**qso_dict[qso])
            continue
        else:
            # ADOM Columns
            ion_dict = {}
            for kk,iis in enumerate(isplit[3:-1]):
                is2 = iis.strip()
                if is2[0:3] == '$>$':
                    ion_dict[ioncs[kk]] = dict(sig_clm=0.,flg_clm=2,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[3:])
                elif is2[0:3] == '$<$':
                    ion_dict[ioncs[kk]] = dict(sig_clm=0.,flg_clm=3,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[3:])
                elif len(is2) == 0:
                    pass
                else:
                    ion_dict[ioncs[kk]] = dict(flg_clm=1,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[0:5])
                    ion_dict[ioncs[kk]]['sig_clm'] = float(is2[10:])
            # Finish
            lls._ionN = pyiau.dict_to_ions(ion_dict)
            lls.Refs.append('Mei07')
            all_lls.append(lls)

    # Return SLLS only
    fin_slls = [ills for ills in all_lls if ills.NHI < 20.3]
    return fin_slls

def meiring08():
    '''Meiring et al. 2008, MNRAS, 384, 1015
    SLLS with Magellan
    Abundances from Table 3 from astro-ph (LateX) by JXP [AODM]
    RA/DEC from Table 1
    Threw out Q1436-0051B given NHI < 18.8 [i.e. unknown]
    ''' 
    all_lls = []
    # Table 1
    tab_fil = pyigm_path+"/data/LLS/Literature/meiring08.tb1.ascii"
    with open(tab_fil,'r') as f:
        flines1 = f.readlines()
    # Grab RA/DEC
    qso_dict = {}
    for iline in flines1:
        if iline[0:2] in ['QS','\h','$\\', 'J2', '  ']:
            continue
        # Parse
        isplit = iline.split('&')
        radec = isplit[3].strip()+isplit[4].strip()
        radec = radec.replace(':','')
        # zem
        zem = float(isplit[5].strip())
        # Save
        qso_dict[isplit[0].strip()] = dict(radec=radec,zem=zem,
                                           vlim=[-500,500.]*u.km/u.s)

    # Abundances (AODM)
    # Table 3
    tab_fil = pyigm_path+"/data/LLS/Literature/meiring08.tb3.ascii"
    with open(tab_fil,'r') as f:
        flines3 = f.readlines()
    #
    for iline in flines3:
        if iline[0:2] in ['\h','  ']:
            continue
        # Parse
        isplit = iline.split('&')
        # Ions
        if iline[0:3] == ' QS':
            ioncs = []
            Zions = []
            for iis in isplit[3:-1]: # Skipping HI
                # Parse
                #is2 = iis.split('\\')
                #ip2 = is2[2].find('}')
                ionc = iis.strip()
                # Zion
                Zion = ltai.name_ion(ionc)
                # Append
                ioncs.append(ionc)
                Zions.append(Zion)
            continue
        if iline[0:2] == ' Q':
            # QSO
            qso = isplit[0].strip()
            if qso[-1] in ['A','B']:
                qso = qso[0:-1]
            # zabs and name
            zabs = float(isplit[1].strip())
            qso_dict[qso]['name']=qso+'z_{:.3f}'.format(zabs) 
            qso_dict[qso]['zabs']=zabs
            # NHI
            is2 = isplit[2].strip()
            if is2[0] == '$':
                qso_dict[qso]['NHI'] = 99.99 # THROW OUT Q1436-0051B
                qso_dict[qso]['sig_NHI'] = np.array([0.,0.])
            else:
                qso_dict[qso]['NHI'] = float(is2[0:5])
                qso_dict[qso]['sig_NHI'] = np.array([float(is2[10:])]*2)
            #if qso_dict[qso]['NHI'] >= 20.3:
            #    print('Uh oh.  DLA')
            # Generate LLS
            lls = LLSSystem(**qso_dict[qso])
            continue
        else:
            # ADOM Columns
            ion_dict = {}
            for kk,iis in enumerate(isplit[3:-1]):
                is2 = iis.strip()
                if is2[0:3] == '$>$':
                    ion_dict[ioncs[kk]] = dict(sig_clm=0.,flg_clm=2,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[3:])
                elif is2[0:3] == '$<$':
                    ion_dict[ioncs[kk]] = dict(sig_clm=0.,flg_clm=3,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[3:])
                elif len(is2) == 0:
                    pass
                else:
                    ion_dict[ioncs[kk]] = dict(flg_clm=1,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[0:5])
                    ion_dict[ioncs[kk]]['sig_clm'] = float(is2[10:])
            # Finish
            lls._ionN = pyiau.dict_to_ions(ion_dict)
            lls.Refs.append('Mei08')
            all_lls.append(lls)

    # Return SLLS only
    fin_slls = [ills for ills in all_lls if ills.NHI < 20.3]
    return fin_slls

def nestor08():
    '''Nestor, D. et al. 2008, MNRAS, 390, 1670-1682
    Q2149+212
    Taken from Table 1 by JXP
    NHI from RTN06 (damping wings)
    RA/DEC from STIS header
    ''' 
    # Setup
    lls = LLSSystem(name='SDSSJ2151+2130_z1.002', radec=(327.94096*u.deg,
                                                         21.503750*u.deg),
                    zem=1.534, zabs=1.0023,
                    vlim=[-300., 300.]*u.km/u.s, NHI=19.30,
                    sig_NHI=np.array([0.10,0.10]))
    #  Meiring06, Table 4
    ion_dict = {}
    ion_dict['Zn II'] = dict(clm=12.13, sig_clm=0.,flg_clm=3,Z=30,ion=2)
    ion_dict['Cr II'] = dict(clm=12.59, sig_clm=0.,flg_clm=3,Z=24,ion=2)
    # Finish
    lls._ionN = pyiau.dict_to_ions(ion_dict)
    lls.Refs.append('Nes08')
    return lls


def meiring09():
    '''Meiring et al. 2009, MNRAS, 393, 1513
    SLLS with Magellan
    Abundances from Table 3 from astro-ph (LateX) by JXP [AODM]
    RA/DEC from Table 1
    ''' 
    all_lls = []
    # Table 1
    tab_fil = pyigm_path+"/data/LLS/Literature/meiring09.tb1.ascii"
    with open(tab_fil,'r') as f:
        flines1 = f.readlines()
    # Grab RA/DEC
    qso_dict = {}
    for iline in flines1:
        if iline[0:3] in [' QS','\hl','$\\c', ' J2', '   ']:
            continue
        # Parse
        isplit = iline.split('&')
        #xdb.set_trace()
        if '$' in isplit[3].strip():
            isplit[3] = '-'+(isplit[3].strip())[3:]
        radec = isplit[2].strip()+isplit[3].strip()
        radec = radec.replace(':','')
        # zem
        zem = float(isplit[5].strip())
        # Save
        qso_dict[isplit[0].strip()] = dict(radec=radec,zem=zem,
                                           vlim=[-500,500.]*u.km/u.s)

    # Abundances (AODM)
    # Table 3
    tab_fil = pyigm_path+"/data/LLS/Literature/meiring09.tb3.ascii"
    with open(tab_fil,'r') as f:
        flines3 = f.readlines()
    #
    for iline in flines3:
        if iline[0:2] in ['\h','  ']:
            continue
        # Parse
        isplit = iline.split('&')
        # Ions
        if iline[0:2] == 'QS':
            ioncs = []
            Zions = []
            for iis in isplit[3:-1]: # Skipping HI
                # Parse
                #is2 = iis.split('\\')
                #ip2 = is2[2].find('}')
                ionc = iis.strip()
                # Zion
                Zion = ltai.name_ion(ionc)
                # Append
                ioncs.append(ionc)
                Zions.append(Zion)
            continue
        if iline[0] == 'Q':
            # QSO
            qso = isplit[0].strip()
            if qso[-1] in ['A','B','C']:
                qso = qso[0:-1]
            # zabs and name
            zabs = float(isplit[1].strip())
            qso_dict[qso]['name']=qso+'z_{:.3f}'.format(zabs) 
            qso_dict[qso]['zabs']=zabs
            # NHI
            is2 = isplit[2].strip()
            if is2[0] == '$':
                qso_dict[qso]['NHI'] = 99.99 # THROW OUT Q1436-0051B
                qso_dict[qso]['sig_NHI'] = np.array([0.,0.])
            else:
                qso_dict[qso]['NHI'] = float(is2[0:5])
                qso_dict[qso]['sig_NHI'] = np.array([float(is2[10:])]*2)
            #if qso_dict[qso]['NHI'] >= 20.3:
            #    print('Uh oh.  DLA')
            # Generate LLS
            lls = LLSSystem(**qso_dict[qso])
            continue
        else:
            # ADOM Columns
            ion_dict = {}
            for kk,iis in enumerate(isplit[3:-1]):
                is2 = iis.strip()
                if is2[0:3] == '$>$':
                    ion_dict[ioncs[kk]] = dict(sig_clm=0.,flg_clm=2,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[3:])
                elif is2[0:3] == '$<$':
                    ion_dict[ioncs[kk]] = dict(sig_clm=0.,flg_clm=3,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[3:])
                elif len(is2) == 0:
                    pass
                else:
                    ion_dict[ioncs[kk]] = dict(flg_clm=1,Z=Zions[kk][0],ion=Zions[kk][1])
                    ion_dict[ioncs[kk]]['clm'] = float(is2[0:5])
                    ion_dict[ioncs[kk]]['sig_clm'] = float(is2[10:])
            # Finish
            lls._ionN = pyiau.dict_to_ions(ion_dict)
            lls.Refs.append('Mei09')
            all_lls.append(lls)

    # Return SLLS only
    fin_slls = [ills for ills in all_lls if ills.NHI < 20.3]
    return fin_slls

def dessauges09():
    '''Dessauges-Zavadsky et al. 2009, MNRAS, 396, L96
    SLLS with UVES
    Zn,Fe abundances from Table 1 from astro-ph (LateX) by JXP [AODM]
     Taken from the Zn/H and Fe/H assuming *no* ionization corrections
    RA/DEC from the 'other' name 
    ''' 
    # Solar abundances
    eZn = 4.63
    eFe = 7.45
    sol = [eFe,eZn]
    #
    all_lls = []
    # Table 1
    tab_fil = pyigm_path+"/data/LLS/Literature/dessauges09.tb1.ascii"
    with open(tab_fil,'r') as f:
        flines1 = f.readlines()
    # Trim the first few lines
    flines1 = flines1[3:]
    for iline in flines1:
        # Parse
        isplit = iline.split('&')
        # QSO
        if iline[0:2] == 'QS':
            # QSO, RA/DEC, zem
            qso = isplit[0][4:].strip()
            radec = isplit[1].strip()[1:].replace('$','')
            zem = float(isplit[3].strip())
        # NHI, zabs
        zabs = float(isplit[4].strip())
        is2 = isplit[6].strip()
        NHI = float(is2[1:6])
        sigNHI = np.array([float(is2[10:14])]*2)
        # name
        name = qso+'z_{:.3f}'.format(zabs) 
        lls = LLSSystem(name=name, radec=radec, vlim=[-500,500]*u.km/u.s,
            zem=zem, zabs=zabs, NHI=NHI, sig_NHI=sigNHI)
        # ADOM Columns
        ion_dict = {}
        for kk,ion in enumerate(['Fe II','Zn II']):
            Zion = ltai.name_ion(ion)
            is2 = isplit[7+kk].strip()
            if is2[0:2] == '$>':
                ion_dict[ion] = dict(sig_clm=0.,flg_clm=2,Z=Zion[0],ion=Zion[1])
                ion_dict[ion]['clm'] = float(is2[2:7]) + NHI - 12 + sol[kk]
            elif is2[0:2] == '$<':
                ion_dict[ion] = dict(sig_clm=0.,flg_clm=3,Z=Zion[0],ion=Zion[1])
                ion_dict[ion]['clm'] = float(is2[2:7]) + NHI - 12 + sol[kk]
            elif is2[0:2] == '..':
                pass
            else:
                ion_dict[ion] = dict(flg_clm=1,Z=Zion[0],ion=Zion[1])
                ion_dict[ion]['clm'] = float(is2[1:6]) + NHI - 12 + sol[kk]
                ion_dict[ion]['sig_clm'] = float(is2[10:14])
        #xdb.set_trace()
        # Finish
        lls._ionN = pyiau.dict_to_ions(ion_dict)
        lls.Refs.append('DZ09')
        all_lls.append(lls)

    # Return SLLS only
    fin_slls = [ills for ills in all_lls if ills.NHI < 20.3]
    return fin_slls

def tumlinson11():
    """Tumlinson, J. et al. 2011, ApJ, 733, 111
    J1009+0713
    HST/COS
    Metal columns parsed from Table 1
    NHI from LL+Lyman series (uncertain)
    """
    # Grab ASCII file from ApJ
    tab_fil = pyigm_path+"/data/LLS/Literature/tumlinson11.tb1.ascii"
    url = 'http://iopscience.iop.org/0004-637X/733/2/111/suppdata/apj388927t1_ascii.txt'
    chk_fil = glob.glob(tab_fil)
    if len(chk_fil) > 0:
        tab_fil = chk_fil[0]
    else:
        print('LLSSurvey: Grabbing table file from {:s}'.format(url))
        f = urllib2.urlopen(url)
        with open(tab_fil, "wb") as code:
            code.write(f.read())
    # Setup
    radec = '100902.06+071343.8'  # From paper
    lls = LLSSystem(name='J1009+0713_z0.356', radec=radec, zem=0.456,
        zabs=0.3558, vlim=[-200., 250.]*u.km/u.s, NHI=18.4, 
        sig_NHI=np.array([0.41,0.41]))

    # Columns
    # Start with Table 3 (VPFIT)
    with open(tab_fil,'r') as f:
        flines1 = f.readlines()
    # Trim
    flines1 = flines1[18:]
    #
    ion_dict = {}
    line_dict = dict(OI='1302',OVI='1038',MgII='2803^b',SiII='1190',
        CaII='3934',FeII='2586')
    ion = None
    for iline in flines1:
        isplit = iline.split('\t')
        if ion=='FeIII': # Last line
            break
        # Ion
        is2 = isplit[0].split(' ')
        ion = is2[0]+is2[1]
        try:
            gdl = line_dict[ion]
        except:
            pass
            #print('Taking {:s}'.format(isplit[0]))
        else:
            if is2[2] != gdl:
                continue
        Zion = ltai.name_ion(ion)
        ion_dict[ion] = dict(logN=0., sig_logN=0., flag_N=0, Z=Zion[0],ion=Zion[1])
        # Combine components [could replace with SubSystems some day]
        for iis in isplit[1:-1]:
            # Upper limit
            if (iis.strip()[0] == '<') & (ion_dict[ion]['flag_N']==0):
                ion_dict[ion]['flag_N']=3
                ion_dict[ion]['logN']=float(iis[1:])
            elif (iis.strip()[0] == '>'): # Saturated
                ion_dict[ion]['flag_N']=2
                ion_dict[ion]['logN']=log_sum([ion_dict[ion]['logN'],float(iis[1:5])])
            elif iis.strip()[0] in ['.','<']:
                pass
            else:
                if ion_dict[ion]['flag_N']==2: # Add to saturated
                    ion_dict[ion]['logN']=log_sum([ion_dict[ion]['logN'],float(iis[0:4])])
                else:
                    ion_dict[ion]['flag_N']=1
                    obj = dict(logN=float(iis[0:4]),sig_logN=float(iis[-4:]),
                               flag_N=1)
                    # Add
                    flag,N,sig = ltaa.sum_logN(ion_dict[ion],obj)
                    ion_dict[ion]['logN']=N
                    ion_dict[ion]['sig_logN']=sig
    # Finish
    lls._ionN = pyiau.dict_to_ions(ion_dict)
    lls.Refs.append('Tum11')
    return lls

def kacprzak12():
    '''Kacprzak, G. et al. 2012, MNRAS, 427, 3029-3043
    TON 153
    Taken from Table 1 by JXP
    NHI from Churchill+2007
    RA/DEC from Simbad
    ''' 
    # Setup
    radec = '131956.2209+272808.271'
    lls = LLSSystem(name='TON153_z1.002', radec=radec, zem=0.6610,
        zabs=1.0023, vlim=[-250., 200.]*u.km/u.s, NHI=18.30, sig_NHI=np.array([0.30,0.30]))
    # Table 1 (total)
    ion_dict = {}
    ion_dict['Mg II'] = dict(clm=13.11, sig_clm=0.07,flg_clm=1,Z=12,ion=2)
    ion_dict['Mg I'] = dict(clm=11.54, sig_clm=0.06,flg_clm=1,Z=12,ion=1)
    ion_dict['Si I'] = dict(clm=11.8, sig_clm=0.00,flg_clm=3,Z=14,ion=1)
    ion_dict['Si II'] = dict(clm=13.16, sig_clm=0.11,flg_clm=1,Z=14,ion=2)
    ion_dict['Si IV'] = dict(clm=12.4, sig_clm=0.0,flg_clm=3,Z=14,ion=4)
    ion_dict['C II'] = dict(clm=13.39, sig_clm=0.0,flg_clm=2,Z=6,ion=2)
    ion_dict['C III'] = dict(clm=14.20, sig_clm=0.05,flg_clm=1,Z=6,ion=3)
    ion_dict['C III'] = dict(clm=14.41, sig_clm=0.05,flg_clm=1,Z=6,ion=4)
    ion_dict['O VI'] = dict(clm=14.49, sig_clm=0.05,flg_clm=1,Z=8,ion=6)
    # Finish
    lls._ionN = pyiau.dict_to_ions(ion_dict)
    lls.Refs.append('Kcz12')
    return lls

def battisti12():
    '''Battisti, A. et al. 2012, ApJ, 744, 93
    HST/COS
    QSO info from Table 1
    Metal columns parsed from Table 3
    NHI from Lya
    '''
    all_lls = []
    # Grab ASCII files from ApJ
    tab_fils = [pyigm_path+"/data/LLS/Literature/battisti12.tb1.ascii",
                pyigm_path+"/data/LLS/Literature/battisti12.tb3.ascii"]
    urls = ['http://iopscience.iop.org/0004-637X/744/2/93/suppdata/apj413924t1_ascii.txt',
        'http://iopscience.iop.org/0004-637X/744/2/93/suppdata/apj413924t3_ascii.txt']
    for jj,tab_fil in enumerate(tab_fils):
        chk_fil = glob.glob(tab_fil)
        if len(chk_fil) > 0:
            tab_fil = chk_fil[0]
        else:
            url = urls[jj]
            print('LLSSurvey: Grabbing table file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            with open(tab_fil, "wb") as code:
                code.write(f.read())
    # QSO info 
    with open(tab_fils[0],'r') as f:
        flines1 = f.readlines()
    # Grab RA/DEC
    all_idict = []
    for iline in flines1:
        if iline[0:2] != 'SD':
            continue
        # Parse
        isplit = iline.split('\t')
        name = isplit[0].split(' ')[1]
        radec = name[1:]
        zem = float(isplit[1].strip())
        zabs = float(isplit[2].strip())
        NHI = float(isplit[3].strip()[0:4])
        sigNHI = np.array([float(isplit[3].strip()[11:])]*2)
        # Save
        lls = LLSSystem(name=name,radec=radec,zem=zem,
            zabs=zabs,NHI=NHI,sig_NHI=sigNHI, vlim=[-500,500]*u.km/u.s)
        #
        all_lls.append(lls)
        all_idict.append({})

    # Abundances
    with open(tab_fils[1],'r') as f:
        flines3 = f.readlines()
    flines3 = flines3[5:]
    ion = None
    for iline in flines3:
        if ion == 'Ni II':
            break
        isplit = iline.split('\t')
        if isplit[0] == 'C II*': # Skipping CII*
            continue
        # ion
        ipos = -1
        while (isplit[0][ipos] not in ['I','V']):
            ipos -= 1
        ion = isplit[0][0:ipos+1+len(isplit[0])]
        Zion = ltai.name_ion(ion)
        # Loop on systems
        for kk,iis in enumerate(isplit[1:-1]):
            if iis.strip()[0] == '.':
                continue
            all_idict[kk][ion] = dict(Z=Zion[0], ion=Zion[1],sig_clm=0.)
            if iis[0] == '>':
                all_idict[kk][ion]['flg_clm'] = 2
                all_idict[kk][ion]['clm'] = float(iis[1:6])
            elif iis[0] == '<':
                all_idict[kk][ion]['flg_clm'] = 3
                all_idict[kk][ion]['clm'] = float(iis[1:])
            else:
                all_idict[kk][ion]['flg_clm'] = 1
                all_idict[kk][ion]['clm'] = float(iis[0:5])
                all_idict[kk][ion]['sig_clm'] = float(iis[-4:])

    # Return SLLS only
    for kk,lls in enumerate(all_lls):
        try:
            lls._ionN = pyiau.dict_to_ions(all_idict[kk])
        except ValueError:
            pdb.set_trace()
        lls.Refs.append('Bat12')
    fin_slls = [ills for ills in all_lls if ills.NHI < 20.3]
    return fin_slls

def lehner13():
    """Lenher et al. 2013
    Low z LLS
    HST
    Metals parsed from Table 2
    """
    # Grab ASCII file from ApJ
    tab_fil = pyigm_path+"/data/LLS/Literature/lehner13_tab2.ascii"
    chk_fil = glob.glob(tab_fil)
    if len(chk_fil) > 0:
        tab_fil = chk_fil[0]
    else:
        url = 'http://iopscience.iop.org/0004-637X/770/2/138/suppdata/apj472363t2_mrt.txt'
        print('LLSSurvey: Grabbing table file from {:s}'.format(url))
        f = urllib2.urlopen(url)
        with open(tab_fil, "wb") as code:
            code.write(f.read())

def wotta16():
    """ Generate sys files from IDL save files
    Returns
    -------
    """
    from scipy.io import readsav
    # Non-excluded
    all = readsav(pyigm_path+'/data/LLS/Literature/wotta16_final.save')

    # Build Lehner+13
    assert False # Need RA/DEC


#####
def log_sum(logN):
    '''Sum up logN values return the log
    '''
    Nsum = np.sum(10.**np.array(logN))
    return np.log10(Nsum)


def load_lls_lit():
    # TODO
    #  Handle duplicates (e.g. DZ vs M08,M09)

    # Begin
    lls_lit = LLSSurvey()
    mask = []

    # Zonak 2004
    Zon04 = zonak2004()
    lls_lit._abs_sys.append(Zon04)
    lls_lit.ref = Zon04.Refs[0]
    mask += [True]
    # Jenkins 2005
    Jen05 = jenkins2005()
    lls_lit._abs_sys.append(Jen05)
    lls_lit.ref = lls_lit.ref + ',' + Jen05.Refs[0]
    mask += [True]
    # Tripp 2005
    Tri05 = tripp2005()
    lls_lit._abs_sys.append(Tri05)
    lls_lit.ref = lls_lit.ref + ',' + Tri05.Refs[0]
    mask += [True]
    # Peroux 2006a
    Prx06a = peroux06a()
    lls_lit._abs_sys.append(Prx06a)
    lls_lit.ref = lls_lit.ref + ',' + Prx06a.Refs[0]
    mask += [True]
    # Peroux 2006b
    Prx06b = peroux06b()
    lls_lit._abs_sys.append(Prx06b)
    lls_lit.ref = lls_lit.ref + ',' + Prx06b.Refs[0]
    mask += [True]
    # Meiring 2006
    Mei06 = meiring06()
    lls_lit._abs_sys.append(Mei06)
    lls_lit.ref = lls_lit.ref + ',' + Mei06.Refs[0]
    mask += [True]
    # Meiring 2007
    Mei07 = meiring07()
    for ills in Mei07:
        lls_lit._abs_sys.append(ills)
        mask += [True]
    lls_lit.ref = lls_lit.ref + ',' + ills.Refs[0]
    # Meiring 2008
    Mei08 = meiring08()
    for ills in Mei08:
        lls_lit._abs_sys.append(ills)
        mask += [True]
    lls_lit.ref = lls_lit.ref + ',' + ills.Refs[0]
    # Nestor 2008
    Nes08 = nestor08()
    lls_lit._abs_sys.append(Nes08)
    lls_lit.ref = lls_lit.ref + ',' + Nes08.Refs[0]
    mask += [True]
    # Meiring 2009
    Mei09 = meiring09()
    for ills in Mei09:
        lls_lit._abs_sys.append(ills)
        mask += [True]
    lls_lit.ref = lls_lit.ref + ',' + ills.Refs[0]
    # Dessauges-Zavadsky 2009
    DZ09 = dessauges09()
    for ills in DZ09:
        lls_lit._abs_sys.append(ills)
        mask += [True]
    lls_lit.ref = lls_lit.ref + ',' + ills.Refs[0]
    # Tumlinson 2011
    Tum11 = tumlinson11()
    lls_lit._abs_sys.append(Tum11)
    lls_lit.ref = lls_lit.ref + ',' + Tum11.Refs[0]
    mask += [True]
    # Kacprzak 2012
    Kcz12 = kacprzak12()
    lls_lit._abs_sys.append(Kcz12)
    mask += [True]
    lls_lit.ref = lls_lit.ref + ',' + Kcz12.Refs[0]
     # Dessauges-Zavadsky 2009
    Bat12 = battisti12()
    for ills in Bat12:
        lls_lit._abs_sys.append(ills)
        mask += [True]
    lls_lit.ref = lls_lit.ref + ',' + ills.Refs[0]

    # Final Mask
    lls_lit.mask = np.array(mask)
    # Return
    return lls_lit
