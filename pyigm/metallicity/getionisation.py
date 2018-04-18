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
Compute the X_ion ionization correction along the chains

"""

def getions_x_ICF(result,ion,base,interp):

    #print('Get XSiII  for ', result['info']['name'])
    xhichain=interp['HI'](result['pdfs'])
    x_chain=interp[ion](result['pdfs'])
    diff=x_chain-xhichain
    hist,edge=np.histogram(diff,bins=base)
    hist=hist/(1.*len(diff))

    return (hist,np.median(diff))


"""
Compute the X_ion along the chains

"""

def getions_x(result,ion,base,interp):

    #print('Get XSiII  for ', result['info']['name'])
    xchain=interp[ion](result['pdfs'])
    hist,edge=np.histogram(xchain,bins=base)
    hist=hist/(1.*len(xchain))
    return (hist,np.median(xchain))



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

def init_interpolator(grid, ions=["HI", "SiII", "SiIII", "CII", "CIII", "CIV", "MgII"]):
    
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


def getionisation(grid_interp, infile_list, outfile, ions=["HI", "SiII", "SiIII", "CII", "CIII", "CIV", "MgII"]):
    
    ##Aliases of input variable name, so it's clearer below
    ionlist = ions

    #store arrays for distributions
    allname=[]
    allz=[]
    allnhi=[]

    #1D histo
    wdt=0.15
    
    U_hist=[]
    U_medians=[]
    U_hist_bin=np.arange(-7.,1.,wdt)
    
    temperature_hist=[]
    temperature_medians=[]
    temperature_hist_bin=np.arange(1.,5.,wdt)

    ############
    
    x_hist = {}
    x_medians = {}
    x_hist_bin = {}
    for ion in ionlist:
        x_hist[ion] = []
        x_medians[ion] = []
    
    x_hist_bin["HI"] = np.arange(-6.,0.4,wdt)
    x_hist_bin["SiII"] = np.arange(-3.,3.,wdt)
    x_hist_bin["SiIII"] = np.arange(-6.,0.4,wdt)
    x_hist_bin["CII"] = np.arange(-6.,0.4,wdt)
    x_hist_bin["CIII"] = np.arange(-6.,0.4,wdt)
    x_hist_bin["CIV"] = np.arange(-6.,0.4,wdt)
    x_hist_bin["MgII"] = np.arange(-6.,0.4,wdt)

    ############

    #compress results
    for fit in infile_list:

        #open the grid 
        try:
            ##Python2
            fil=open(fit)
            modl=pickle.load(fil)
        except:
            ##Python3
            fil=open(fit,'rb')
            modl=pickle.load(fil, encoding='latin1')
        
        fil.close()
        

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
            
            ############
            
            for ion in ionlist:
                this_X, med_X = getions_x(result,ion,x_hist_bin[ion],grid_interp)
                x_hist[ion].append(this_X/wdt)
                x_medians[ion].append(med_X)
            
            ############
            
            try:
                #grab the temperature histogram 
                this_temperature,med_temperature=gettemperature(result,temperature_hist_bin,grid_interp)
                temperature_hist.append(this_temperature/wdt)
                # temperature_hist.append(this_temperature)
                temperature_medians.append(med_temperature)
                temperature_on = True
            except:
                temperature_on = False
        
        else:
            print("Skipping ",fit," because it is not within the z and/or NHI range.")
    
    outdict = {}
    outdict['name'] = allname
    outdict['z'] = allz
    outdict['logNHI'] = allnhi
    
    outdict['U_base'] = U_hist_bin
    outdict['U_hist'] = U_hist
    outdict['U_med'] = U_medians
    for ion in ionlist:
        keyname_base = "X{}_base".format(ion)
        keyname_hist = "X{}_hist".format(ion)
        keyname_med = "X{}_med".format(ion)
        outdict[keyname_base] = x_hist_bin[ion]
        outdict[keyname_hist] = x_hist[ion]
        outdict[keyname_med] = x_medians[ion]
    
    
    if temperature_on is True:
        outdict['temperature_base'] = temperature_hist_bin
        outdict['temperature_hist'] = temperature_hist
    
    
    #save
    pickle.dump(outdict,open(outfile,"w"))
    
    
    return



