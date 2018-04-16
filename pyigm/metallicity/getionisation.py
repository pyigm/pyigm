from COSCGMLegacyPath_define import COSCGMLegacyPath
import glob
import pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
##CBW
import h5py
from matplotlib.ticker import NullFormatter
matplotlib.rcParams.update({'font.size': 16})
from scipy import interpolate
import math

mblue='#1E90FF'
mred='#DC143C'

from astropy import constants as const
from scipy import integrate 
from scipy.interpolate import interp1d


#import UVB interpolator
from pyigm.euvb import cloudyiongrid as cld

# print("...Done")




"""

Compute the ionization parameter along the chains - assumes HM05

"""

def getions_uparam(result,base):
    
    #print('Get ionization parameter for ', result['info']['name'])
    
    #init UVB at given redshift
    uvb=cld.Cldyion(uvb='HM05logU')
    uvb.redshift=result['info']['z']
    
    ##This should be taken care of when calling cld.Cldyion(uvb='HM05logU')
    ##  but I'm going to do this again, just in case.
    uvb.uvbtype='HM05logU'
    #Get the freqeunce Log of Hz and UVB Log of erg/s/cm2/Hz
    uvb.inituvb()
    
    #now integrate 4pi Jnu / h  in d log nu
    sort=np.argsort(uvb.source['lgnu'])
    uvb.source['lgnu']=uvb.source['lgnu'][sort]
    uvb.source['lgfpiJ']=uvb.source['lgfpiJ'][sort]
    
    #define integal quantity (integrate in d log nu)
    lognu=np.log(10**uvb.source['lgnu'])
    hplanck = 6.6260755e-27 # erg/s
    integrand=10**uvb.source['lgfpiJ']/hplanck
    
    #Define min-max in log_natural Hz
    maxnu=np.max(lognu)
    ionnu=np.log((const.c.to('cm/s')/(912*1e-8)).value)
    
    #integrate [checked against cloudy output]
    fint = interp1d(lognu,integrand)
    phi,err = integrate.quad(fint,ionnu,maxnu)
    
    #now compute the ionization paramater
    den=result['tags'].index('dens')
    Uparam=np.log10(phi)-result['pdfs'][:,den]-np.log10(const.c.to('cm/s').value)
    
    hist,edge=np.histogram(Uparam,bins=base)
    hist=hist/(1.*len(result['pdfs'][:,den]))
    
    return (hist, np.median(Uparam))



def getions_uparam_unbinned(result):
    ##Added by CBW
    
    #print('Get ionization parameter for ', result['info']['name'])
    
    #init UVB at given redshift
    uvb=cld.Cldyion(uvb='HM05logU')
    uvb.redshift=result['info']['z']
    
    ##This should be taken care of when calling cld.Cldyion(uvb='HM05logU')
    ##  but I'm going to do this again, just in case.
    uvb.uvbtype='HM05logU'
    #Get the freqeunce Log of Hz and UVB Log of erg/s/cm2/Hz
    uvb.inituvb()
    
    #now integrate 4pi Jnu / h  in d log nu
    sort=np.argsort(uvb.source['lgnu'])
    uvb.source['lgnu']=uvb.source['lgnu'][sort]
    uvb.source['lgfpiJ']=uvb.source['lgfpiJ'][sort]
    
    #define integal quantity (integrate in d log nu)
    lognu=np.log(10**uvb.source['lgnu'])
    hplanck = 6.6260755e-27 # erg/s
    integrand=10**uvb.source['lgfpiJ']/hplanck
    
    #Define min-max in log_natural Hz
    maxnu=np.max(lognu)
    ionnu=np.log((const.c.to('cm/s')/(912*1e-8)).value)
    
    #integrate [checked against cloudy output]
    fint = interp1d(lognu,integrand)
    phi,err = integrate.quad(fint,ionnu,maxnu)
    
    #now compute the ionization paramater
    den=result['tags'].index('dens')
    Uparam=np.log10(phi)-result['pdfs'][:,den]-np.log10(const.c.to('cm/s').value)
    
    # hist,edge=np.histogram(Uparam,bins=base)
    # hist=hist/(1.*len(result['pdfs'][:,den]))
    
    # return (hist, np.median(Uparam))
    return Uparam



"""

Compute the XHI along the chains

"""

def getions_xhi(result,base,interp):

    #print('Get XHI  for ', result['info']['name'])

    #grab XHI [tested against cloudy]
    xhichain=interp['HI'](result['pdfs'])
    hist,edge=np.histogram(xhichain,bins=base)
    hist=hist/(1.*len(xhichain))

    return (hist,np.median(xhichain))


"""
Compute the X_MgII ionization correction along the chains

"""

def getions_xmg2(result,base,interp):

    #print('Get XSiII  for ', result['info']['name'])
    xhichain=interp['HI'](result['pdfs'])
    xsi2chain=interp['MgII'](result['pdfs'])
    diff=xsi2chain-xhichain
    hist,edge=np.histogram(diff,bins=base)
    hist=hist/(1.*len(diff))

    return (hist,np.median(diff))

"""
Compute the X_SiII ionization correction along the chains

"""

def getions_xsi2(result,base,interp):

    #print('Get XSiII  for ', result['info']['name'])
    xhichain=interp['HI'](result['pdfs'])
    xsi2chain=interp['SiII'](result['pdfs'])
    diff=xsi2chain-xhichain
    hist,edge=np.histogram(diff,bins=base)
    hist=hist/(1.*len(diff))

    return (hist,np.median(diff))

"""
Compute the X_SiIII along the chains

"""

def getions_xsi3(result,base,interp):

    #print('Get XSiII  for ', result['info']['name'])
    xsi3chain=interp['SiIII'](result['pdfs'])
    hist,edge=np.histogram(xsi3chain,bins=base)
    hist=hist/(1.*len(xsi3chain))
    return (hist,np.median(xsi3chain))



"""
Compute the X_CII along the chains

"""

def getions_xc2(result,base,interp):

    xC2chain=interp['CII'](result['pdfs'])
    hist,edge=np.histogram(xC2chain,bins=base)
    hist=hist/(1.*len(xC2chain))

    return (hist,np.median(xC2chain))

"""
Compute the X_CIII along the chains

"""

def getions_xc3(result,base,interp):

    xC3chain=interp['CIII'](result['pdfs'])
    hist,edge=np.histogram(xC3chain,bins=base)
    hist=hist/(1.*len(xC3chain))

    return (hist,np.median(xC3chain))


"""
Compute the X_CIV along the chains

"""

def getions_xc4(result,base,interp):

    xC4chain=interp['CIV'](result['pdfs'])
    hist,edge=np.histogram(xC4chain,bins=base)
    hist=hist/(1.*len(xC4chain))

    return (hist,np.median(xC4chain))


"""
Compute the temperature along the chains

"""

def gettemperature(result,base,interp):

    temperaturechain=interp['temperature'](result['pdfs'])
    hist,edge=np.histogram(temperaturechain,bins=base)
    hist=hist/(1.*len(temperaturechain))

    return (hist,np.median(temperaturechain))


"""
Initialise the interpolator over the grid

"""

def init_interpolator(grid, ions=['HI','SiII','SiIII','CII','CIII','CIV']):
    
    print("Loading Cloudy grid...")
    
    #open the grid 
    try:
        ##Python2
        fil=open(grid)
        modl=pickle.load(fil)
    except:
        ##Python3
        fil=open(grid,'rb')
        modl=pickle.load(fil, encoding='latin1')
    
    fil.close()
    
    print("...Done")
    
    #define the dimension of the problem
    mod_axistag=modl[0]
    mod_axisval=[]
    ndim=len(mod_axistag)
    nmodels=1
    for tt in mod_axistag:
        nmodels=nmodels*(modl[1][tt]).size
        #append axis value in a list
        mod_axisval.append(modl[1][tt])
    
    
    interp = {}
    for ion in ions:
        try:
            xiongrid = modl[2][ion]
            x_interp=interpolate.RegularGridInterpolator(mod_axisval,xiongrid,method='linear',bounds_error=False,fill_value=-np.inf)
            interp[ion] = x_interp
        except:
            print("{} not in grid. Skipping.\n".format(ion))
    
    try:
        ##CBW: extract temperatures
        temperaturegrid=modl[4]
        #set up interpolator for this
        temperature_interp=interpolate.RegularGridInterpolator(mod_axisval,temperaturegrid,method='linear',bounds_error=False,fill_value=-np.inf)
        
        #return
        interp['temperature'] = temperature_interp
    except:
        pass
    
    return interp


################################################
################################################


def getionization(sample,grid_interp):
    
    print("Computing log U for "+sample+" sample...")
    
    ##CBW adjust
    #define the model
    mcmcoutputlist='{}/Strong_LLSs/Figures/plot-logU-{}_MCMC.list'.format(COSCGMLegacyPath, sample)
    outfile='{}/Strong_LLSs/Figures/plot-logU-{}_MCMC_ionization.pkl'.format(COSCGMLegacyPath, sample)

    interp = grid_interp

    #load the models
    # figname=mcmcoutputlist
    ##CBW
    # listin=glob.glob(mcmcoutputlist+'/*emcee*pkl')
    listin=[]
    for line in open(mcmcoutputlist):
        li=line.strip()
        ##skip blank lines and skip commented-out lines
        if li and not li.startswith("#"):
            listin.append(COSCGMLegacyPath+'/'+line.rstrip())
    print('Found {} fits'.format(len(listin)))
    
    #store arrays for distributions
    allname=[]
    allz=[]
    allnhi=[]

    #1D histo
    wdt=0.15
    
    U_hist=[]
    U_medians=[]
    U_hist_bin=np.arange(-7.,1.,wdt)
    
    X_hist=[]
    X_medians=[]
    X_hist_bin=np.arange(-6.,0.4,wdt)
    
    Xsi2_hist=[]
    Xsi2_medians=[]
    Xsi2_hist_bin=np.arange(-3.,3.,wdt)
    
    Xsi3_hist=[]
    Xsi3_medians=[]
    Xsi3_hist_bin=np.arange(-6.,0.4,wdt)
    
    XC2_hist=[]
    XC2_medians=[]
    XC2_hist_bin=np.arange(-6.,0.4,wdt)
    
    XC3_hist=[]
    XC3_medians=[]
    XC3_hist_bin=np.arange(-6.,0.4,wdt)
    
    XC4_hist=[]
    XC4_medians=[]
    XC4_hist_bin=np.arange(-6.,0.4,wdt)

    temperature_hist=[]
    temperature_medians=[]
    temperature_hist_bin=np.arange(1.,5.,wdt)

    #compress results
    for fit in listin:

        print('Process ', fit)

        #open the grid 
        try:
            ##Python2
            fl=open(fit)
            modl=pickle.load(fl)
        except:
            ##Python3
            fl=open(fit,'rb')
            modl=pickle.load(fl, encoding='latin1')

        fl.close()

        # result=h5py.File(fit, 'r')

        #grab the indexes of percentile and quantities
        med=result['percent'].index(50)
        p75=result['percent'].index(75)
        p25=result['percent'].index(25)
        col=result['tags'].index('col')
        met=result['tags'].index('met')
        red=result['tags'].index('red')
        den=result['tags'].index('dens')
        
        ##CBW: account for "some allow carbalpha to vary, some do not" problem
        ##Just append a column of "0.0"
        try:
            carbalph=result['tags'].index('carbalpha')
        except:
            carbalph=-999
            old = result['pdfs'].copy()
            new = np.empty([len(old), len(old[0])+1])
            for idx in range(len(new)):
                new[idx] = np.append(old[idx], 0.0)
            
            result['pdfs'] = new.copy()
        
        try:
            jqso=result['tags'].index('jqso')
        except:
            jqso=False
        try:
            jgal=result['tags'].index('jgal')
        except:
            jgal=False
        try: 
            fstar=result['tags'].index('fstar')
        except:
            fstar=False
        try:
            temp=result['tags'].index('temp')
        except: 
            temp=False

    
        ##CBW adjust
        #check condition if needed
        if((result['info']['z'] >= 0.0) & (result['info']['z'] <= 1.5) & (result['info']['NHI'] < 19.0) ):
            # if('LLS' not in figname):
            #     figname=figname+'LLSz'

            #store values of NHI
            allname.append(result['info']['name'])
            allz.append(result['info']['z'])
            allnhi.append(result['info']['NHI'])
    
            #grab the U histogram
            this_U,med_U=getions_uparam(result,U_hist_bin)
            U_hist.append(this_U/wdt)
            # U_hist.append(this_U)
            U_medians.append(med_U)
            
            #grab the X HI histogram 
            this_X,med_X=getions_xhi(result,X_hist_bin,interp)
            X_hist.append(this_X/wdt)
            # X_hist.append(this_X)
            X_medians.append(med_X)
            
            #grab the X Mg2 ionization correction histogram 
            this_Xmg2,med_Xmg2=getions_xmg2(result,Xmg2_hist_bin,interp)
            Xmg2_hist.append(this_Xmg2/wdt)
            # Xmg2_hist.append(this_Xmg2)
            Xmg2_medians.append(med_Xmg2)
            
            #grab the X Si2 ionization correction histogram 
            this_Xsi2,med_Xsi2=getions_xsi2(result,Xsi2_hist_bin,interp)
            Xsi2_hist.append(this_Xsi2/wdt)
            # Xsi2_hist.append(this_Xsi2)
            Xsi2_medians.append(med_Xsi2)
            
            #grab the X Si3 fraction histogram 
            this_Xsi3,med_Xsi3=getions_xsi3(result,Xsi3_hist_bin,interp)
            Xsi3_hist.append(this_Xsi3/wdt)
            # Xsi3_hist.append(this_Xsi3)
            Xsi3_medians.append(med_Xsi3)

            #grab the XC2 histogram 
            this_XC2,med_XC2=getions_xc2(result,XC2_hist_bin,interp)
            XC2_hist.append(this_XC2/wdt)
            # XC2_hist.append(this_XC2)
            XC2_medians.append(med_XC2)
            
            #grab the XC3 histogram 
            this_XC3,med_XC3=getions_xc3(result,XC3_hist_bin,interp)
            XC3_hist.append(this_XC3/wdt)
            # XC3_hist.append(this_XC3)
            XC3_medians.append(med_XC3)

            #grab the XC4 histogram 
            this_XC4,med_XC4=getions_xc4(result,XC4_hist_bin,interp)
            XC4_hist.append(this_XC4/wdt)
            # XC4_hist.append(this_XC4)
            XC4_medians.append(med_XC4)
        
            try:
                #grab the temperature histogram 
                this_temperature,med_temperature=gettemperature(result,temperature_hist_bin,interp)
                temperature_hist.append(this_temperature/wdt)
                # temperature_hist.append(this_temperature)
                temperature_medians.append(med_temperature)
                temperature_on = True
            except:
                temperature_on = False
        
        else:
            print("Skipping ",fit," because it is not within the z and/or NHI range.")

    #save
    if temperature_on is True:
        pickle.dump({'U_base':U_hist_bin,'U_hist':U_hist,'U_med':U_medians,
                     'XHI_base':X_hist_bin,'XHI_hist':X_hist,'XHI_med':X_medians,
                     'XMgII_base':Xmg2_hist_bin,'XMgII_hist':Xmg2_hist,'XMgII_med':Xmg2_medians,
                     'XSiII_base':Xsi2_hist_bin,'XSiII_hist':Xsi2_hist,'XSiII_med':Xsi2_medians,
                     'XSiIII_base':Xsi3_hist_bin,'XSiIII_hist':Xsi3_hist,'XSiIII_med':Xsi3_medians,
                     'XCII_base':XC2_hist_bin,'XCII_hist':XC2_hist,'XCII_med':XC2_medians,
                     'XCIII_base':XC3_hist_bin,'XCIII_hist':XC3_hist,'XCIII_med':XC3_medians,
                     'XCIV_base':XC4_hist_bin,'XCIV_hist':XC4_hist,'XCIV_med':XC4_medians,
                     ##CBW
                     'temperature_base':temperature_hist_bin,'temperature_hist':temperature_hist,
                     'name':allname,'z':allz,'logNHI':allnhi},open(outfile,"w"))
    else:
        pickle.dump({'U_base':U_hist_bin,'U_hist':U_hist,'U_med':U_medians,
                     'XHI_base':X_hist_bin,'XHI_hist':X_hist,'XHI_med':X_medians,
                     'XMgII_base':Xmg2_hist_bin,'XMgII_hist':Xmg2_hist,'XMgII_med':Xmg2_medians,
                     'XSiII_base':Xsi2_hist_bin,'XSiII_hist':Xsi2_hist,'XSiII_med':Xsi2_medians,
                     'XSiIII_base':Xsi3_hist_bin,'XSiIII_hist':Xsi3_hist,'XSiIII_med':Xsi3_medians,
                     'XCII_base':XC2_hist_bin,'XCII_hist':XC2_hist,'XCII_med':XC2_medians,
                     'XCIII_base':XC3_hist_bin,'XCIII_hist':XC3_hist,'XCIII_med':XC3_medians,
                     'XCIV_base':XC4_hist_bin,'XCIV_hist':XC4_hist,'XCIV_med':XC4_medians,
                     'name':allname,'z':allz,'logNHI':allnhi},open(outfile,"w"))


    return



def main():
    # grid='{}/Strong_LLSs/Cloudy/Cloudy_grids/grid_minimal_HM05_carbalpha.pkl'.format(COSCGMLegacyPath)
    grid='{}/Strong_LLSs/Cloudy/Cloudy_grids/grid_cgm_extensive_HM05_carbalpha.pkl'.format(COSCGMLegacyPath)
    
    interp=init_interpolator(grid)

    samples = ["L13","W16","W17","L18"]
    for sample in samples:
        getionization(sample,interp)


if __name__ == '__main__':
    main()

