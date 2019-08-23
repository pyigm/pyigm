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




def getions_uparam(result,base,uvb="HM05logU"):
    """
    Compute the ionization parameter (logU) along the chains
    
    Inputs:
    --------------------------
    result: list of lists (of floats/walkers)
        Usually from the MCMC output (pickle.load('mcmc_out.pkl')),
        with length = number_of_walkers/trials (e.g., 100,000) and
        length[0] = number_of_axes (e.g., NHI, z, [X/H], n_H).
        It is can also be done artificially like:
            result['pdfs'] = [
                [17.2, 0.7, -1.0, -2.3],
                [17.1, 0.7, -1.1, -2.4],
                [17.2, 0.8, -1.0, -2.4],
                [17.2, 0.7, -1.0, -2.3],
            ]
    base: np.array
        Histogram edges to use for logU binning, np.arange(...)
    uvb: str
        Which UVB to integrate (see cloudyiongrid.py)
    
    Outputs:
    --------------------------
    histogram: np.array
        The logU histogram y values (x values are based on "base" input)
    median: float
        The median logU value of this system
    
    """
    
    #print('Get ionization parameter for ', result['info']['name'])
    
    #init UVB at given redshift
    uvb_spectrum=cld.Cldyion(uvb=uvb)
    uvb_spectrum.redshift=result['info']['z']
    
    ##This should be taken care of when calling cld.Cldyion(uvb=uvb)
    ##  but I'm going to do this again, just in case.
    uvb_spectrum.uvbtype=uvb
    #Get the freqeunce Log of Hz and UVB Log of erg/s/cm2/Hz
    uvb_spectrum.inituvb()
    
    #now integrate 4pi Jnu / h  in d log nu
    sort=np.argsort(uvb_spectrum.source['lgnu'])
    uvb_spectrum.source['lgnu']=uvb_spectrum.source['lgnu'][sort]
    uvb_spectrum.source['lgfpiJ']=uvb_spectrum.source['lgfpiJ'][sort]
    
    #define integal quantity (integrate in d log nu)
    lognu=np.log(10**uvb_spectrum.source['lgnu'])
    hplanck = 6.6260755e-27 # erg/s
    integrand=10**uvb_spectrum.source['lgfpiJ']/hplanck
    
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



def getions_uparam_unbinned(result,uvb="HM05logU"):
    """
    Compute the ionization parameter (logU) along the chains, unbinned
    
    Inputs:
    --------------------------
    result: list of lists (of floats/walkers)
        Usually from the MCMC output (pickle.load('mcmc_out.pkl')),
        with length = number_of_walkers/trials (e.g., 100,000) and
        length[0] = number_of_axes (e.g., NHI, z, [X/H], n_H).
        It is can also be done artificially like:
            result['pdfs'] = [
                [17.2, 0.7, -1.0, -2.3],
                [17.1, 0.7, -1.1, -2.4],
                [17.2, 0.8, -1.0, -2.4],
                [17.2, 0.7, -1.0, -2.3],
            ]
    base: np.array
        Histogram edges to use for logU binning, np.arange(...)
    uvb: str
        Which UVB to integrate (see cloudyiongrid.py)
    
    Outputs:
    --------------------------
    logU: np.array
        The logU values for every walker
    
    """
    
    ##Added by CBW
    
    #print('Get ionization parameter for ', result['info']['name'])
    
    #init UVB at given redshift
    uvb_spectrum=cld.Cldyion(uvb=uvb)
    uvb_spectrum.redshift=result['info']['z']
    
    ##This should be taken care of when calling cld.Cldyion(uvb=uvb)
    ##  but I'm going to do this again, just in case.
    uvb_spectrum.uvbtype=uvb
    #Get the freqeunce Log of Hz and UVB Log of erg/s/cm2/Hz
    uvb_spectrum.inituvb()
    
    #now integrate 4pi Jnu / h  in d log nu
    sort=np.argsort(uvb_spectrum.source['lgnu'])
    uvb_spectrum.source['lgnu']=uvb_spectrum.source['lgnu'][sort]
    uvb_spectrum.source['lgfpiJ']=uvb_spectrum.source['lgfpiJ'][sort]
    
    #define integal quantity (integrate in d log nu)
    lognu=np.log(10**uvb_spectrum.source['lgnu'])
    hplanck = 6.6260755e-27 # erg/s
    integrand=10**uvb_spectrum.source['lgfpiJ']/hplanck
    
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



def getions_x_ICF(result,ion,base,interp):
    """
    Compute the x_ion ionization correction (ICF) along the chains
    
    Inputs:
    --------------------------
    result: list of lists (of floats/walkers)
        Usually from the MCMC output (pickle.load('mcmc_out.pkl')),
        with length = number_of_walkers/trials (e.g., 100,000) and
        length[0] = number_of_axes (e.g., NHI, z, [X/H], n_H).
        It is can also be done artificially like:
            result['pdfs'] = [
                [17.2, 0.7, -1.0, -2.3],
                [17.1, 0.7, -1.1, -2.4],
                [17.2, 0.8, -1.0, -2.4],
                [17.2, 0.7, -1.0, -2.3],
            ]
    ion: str
        The ion you want to get the ICF for. This HAS to exist within init_interpolator()
    base: np.array
        Histogram edges to use for logU binning, np.arange(...)
    interp: RegularGridInterpolator
        The interpolated Cloudy grid; see init_interpolator() (below)
    
    Outputs:
    --------------------------
    histogram: np.array
        The histogram y values (x values are based on "base" input)
    median: float
        The median value of this system
    
    """
    

    #print('Get XSiII  for ', result['info']['name'])
    xhichain=interp['HI'](result['pdfs'])
    x_chain=interp[ion](result['pdfs'])
    diff=x_chain-xhichain
    hist,edge=np.histogram(diff,bins=base)
    hist=hist/(1.*len(diff))

    return (hist,np.median(diff))


def getions_x(result,ion,base,interp):

    """
    Compute the x_ion ionization along the chains
    
    Inputs:
    --------------------------
    result: list of lists (of floats/walkers)
        Usually from the MCMC output (pickle.load('mcmc_out.pkl')),
        with length = number_of_walkers/trials (e.g., 100,000) and
        length[0] = number_of_axes (e.g., NHI, z, [X/H], n_H).
        It is can also be done artificially like:
            result['pdfs'] = [
                [17.2, 0.7, -1.0, -2.3],
                [17.1, 0.7, -1.1, -2.4],
                [17.2, 0.8, -1.0, -2.4],
                [17.2, 0.7, -1.0, -2.3],
            ]
    ion: str
        The ion you want to get the ICF for. This HAS to exist within init_interpolator()
    base: np.array
        Histogram edges to use for logU binning, np.arange(...)
    interp: RegularGridInterpolator
        The interpolated Cloudy grid; see cloudyiongrid.py
    
    Outputs:
    --------------------------
    histogram: np.array
        The histogram y values (x values are based on "base" input)
    median: float
        The median value of this system
    
    """
    
    xchain=interp[ion](result['pdfs'])
    hist,edge=np.histogram(xchain,bins=base)
    hist=hist/(1.*len(xchain))
    return (hist,np.median(xchain))



def gettemperature(result,base,interp):
    """
    Compute the temperature along the chains
    
    Inputs:
    --------------------------
    result: list of lists (of floats/walkers)
        Usually from the MCMC output (pickle.load('mcmc_out.pkl')),
        with length = number_of_walkers/trials (e.g., 100,000) and
        length[0] = number_of_axes (e.g., NHI, z, [X/H], n_H).
        It is can also be done artificially like:
            result['pdfs'] = [
                [17.2, 0.7, -1.0, -2.3],
                [17.1, 0.7, -1.1, -2.4],
                [17.2, 0.8, -1.0, -2.4],
                [17.2, 0.7, -1.0, -2.3],
            ]
    base: np.array
        Histogram edges to use for logU binning, np.arange(...)
    interp: RegularGridInterpolator
        The interpolated Cloudy grid; see cloudyiongrid.py
    
    Outputs:
    --------------------------
    histogram: np.array
        The histogram y values (x values are based on "base" input)
    median: float
        The median value of this system
    
    """
    

    temperaturechain=interp['temperature'](result['pdfs'])
    hist,edge=np.histogram(temperaturechain,bins=base)
    hist=hist/(1.*len(temperaturechain))

    return (hist,np.median(temperaturechain))


"""
Initialise the interpolator over the grid

"""

def init_interpolator(grid, ions=["HI", "SiII", "SiIII", "CII", "CIII", "CIV", "MgII"]):
    """
    Initialise the interpolator over the grid
    
    Inputs:
    --------------------------
    grid: str
        The path to the location of the Cloudy grid (*.pkl), output of
        grid_compress.py (or something similar).
    
    Outputs:
    --------------------------
    interpolated_grid: RegularGridInterpolator
        The interpolated Cloudy grid
    
    """
    
    
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
    
    """
    Get the histograms of ionisation for a set of ions (e.g., x_SiII)
    
    Inputs:
    --------------------------
    grid_interp: RegularGridInterpolator
        The interpolated Cloudy grid, output of init_interpolator()
    infile_list: list of str
        A list of the *.pkl files, as outputted by the MCMC code
    outfile: str
        The file you want to save the histograms to
    ions: list of str
        The ions you want to get the x_ion histograms for.
    
    Outputs:
    --------------------------
    outdict: *.pkl pickle file
        File of binned histograms
    
    """
    
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
    
    ##In case the ion we need isn't pre-defined above
    for ion in ionlist:
        try:
            dummy = x_hist_bin[ion]
        except:
            x_hist_bin[ion] = np.arange(-6.,0.4,wdt)

    ############

    #compress results
    for fit in infile_list:

        #open the grid 
        try:
            ##Python2
            fil=open(fit)
            result=pickle.load(fil)
        except:
            ##Python3
            fil=open(fit,'rb')
            result=pickle.load(fil, encoding='latin1')
        
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



