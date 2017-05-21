"""
Reconstruct the ionsation properties of the different systems

"""

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

mblue='#1E90FF'
mred='#DC143C'

from astropy import constants as const
from scipy import integrate 
from scipy.interpolate import interp1d


#import UVB interpolator
import sys
# sys.path.append('../cloudy/py/')
import cloudyiongrid as cld


"""

Compute the ionisation paramater along the chains - assumes HM12

"""

def getions_uparam(result,base):

    #print 'Get ionisation paramater for ', result['info']['name']
  
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
    
    #now compute the ionisation paramater
    den=result['tags'].index('dens')
    Uparam=np.log10(phi)-result['pdfs'][:,den]-np.log10(const.c.to('cm/s').value)
   
    hist,edge=np.histogram(Uparam,bins=base)
    hist=hist/(1.*len(result['pdfs'][:,den]))

    return (hist, np.median(Uparam))


"""

Compute the XHI along the chains

"""

def getions_xhi(result,base,interp):

    #print 'Get XHI  for ', result['info']['name']

    #grab XHI [tested against cloudy]
    xhichian=interp['HI'](result['pdfs'])
    hist,edge=np.histogram(xhichian,bins=base)
    hist=hist/(1.*len(xhichian))

    return (hist,np.median(xhichian))


"""
Compute the X_SiII ionisation correction along the chains

"""

def getions_xsi2(result,base,interp):

    #print 'Get XSiII  for ', result['info']['name']
    xhichian=interp['HI'](result['pdfs'])
    xsi2chian=interp['SiII'](result['pdfs'])
    diff=xsi2chian-xhichian
    hist,edge=np.histogram(diff,bins=base)
    hist=hist/(1.*len(diff))

    return (hist,np.median(diff))

"""
Compute the X_SiIII along the chains

"""

def getions_xsi3(result,base,interp):

    #print 'Get XSiII  for ', result['info']['name']
    xsi3chian=interp['SiIII'](result['pdfs'])
    hist,edge=np.histogram(xsi3chian,bins=base)
    hist=hist/(1.*len(xsi3chian))
    return (hist,np.median(xsi3chian))



"""
Compute the X_CII along the chains

"""

def getions_xc2(result,base,interp):

    xC2chian=interp['CII'](result['pdfs'])
    hist,edge=np.histogram(xC2chian,bins=base)
    hist=hist/(1.*len(xC2chian))

    return (hist,np.median(xC2chian))

"""
Compute the X_CIII along the chains

"""

def getions_xc3(result,base,interp):

    xC3chian=interp['CIII'](result['pdfs'])
    hist,edge=np.histogram(xC3chian,bins=base)
    hist=hist/(1.*len(xC3chian))

    return (hist,np.median(xC3chian))


"""
Compute the X_CIV along the chains

"""

def getions_xc4(result,base,interp):

    xC4chian=interp['CIV'](result['pdfs'])
    hist,edge=np.histogram(xC4chian,bins=base)
    hist=hist/(1.*len(xC4chian))

    return (hist,np.median(xC4chian))


"""
Initialise the interpolator over the grid

"""

def init_interpolator(grid):

    
    #open the grid 
    fil=open(grid)
    modl=pickle.load(fil)
    fil.close()
    
    #define the dimension of the problem
    mod_axistag=modl[0]
    mod_axisval=[]
    ndim=len(mod_axistag)
    nmodels=1
    for tt in mod_axistag:
        nmodels=nmodels*(modl[1][tt]).size
        #append axis value in a list
        mod_axisval.append(modl[1][tt])
   
    #extract x_HI
    xhigrid=modl[2]['HI']
    #set up interpolator for this ions 
    xhi_interp=interpolate.RegularGridInterpolator(mod_axisval,xhigrid,method='linear',bounds_error=False,fill_value=-np.inf)
        
    #extract x_SiII
    xsi2grid=modl[2]['SiII']
    #set up interpolator for this ions 
    xsi2_interp=interpolate.RegularGridInterpolator(mod_axisval,xsi2grid,method='linear',bounds_error=False,fill_value=-np.inf)
 
    #extract x_SiIII
    xsi3grid=modl[2]['SiIII']
    #set up interpolator for this ions 
    xsi3_interp=interpolate.RegularGridInterpolator(mod_axisval,xsi3grid,method='linear',bounds_error=False,fill_value=-np.inf)

    #extract x_CII
    xc2grid=modl[2]['CII']
    #set up interpolator for this ions 
    xc2_interp=interpolate.RegularGridInterpolator(mod_axisval,xc2grid,method='linear',bounds_error=False,fill_value=-np.inf)
    
    #extract x_CIII
    xc3grid=modl[2]['CIII']
    #set up interpolator for this ions 
    xc3_interp=interpolate.RegularGridInterpolator(mod_axisval,xc3grid,method='linear',bounds_error=False,fill_value=-np.inf)

    #extract x_CIV
    xc4grid=modl[2]['CIV']
    #set up interpolator for this ions 
    xc4_interp=interpolate.RegularGridInterpolator(mod_axisval,xc4grid,method='linear',bounds_error=False,fill_value=-np.inf)

    #return
    interp={'HI':xhi_interp,'SiII':xsi2_interp,'SiIII':xsi3_interp,'CII':xc2_interp,'CIII':xc3_interp,'CIV':xc4_interp}
   
    return interp


"""

Main function

"""

if __name__ == "__main__":
    
    ##CBW adjust
    #define the model and grid
    modelsdir='minimal'
    grid='grid_minextsupz_carbalpha.pkl'

    #init the interpolator
    interp=init_interpolator(grid)

    #load the models
    figname=modelsdir
    listin=glob.glob(modelsdir+'/*emcee*hd5')
    print 'Found {} fits'.format(len(listin))
    
    #store arrays for distributions
    allname=[]

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

    #compress results
    for fit in listin:

        print 'Process ', fit

        fl=open(fit)
        result=pickle.load(fl)
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

    
        #check condition if needed
        if((result['info']['z'] >= 2.5) & (result['info']['z'] <= 3.5) & (result['info']['NHI'] < 19.0) ):
            if('LLS' not in figname):
                figname=figname+'LLSz'

            #store values of NHI
            allname.append(result['info']['name'])
    
            #grab the U histogram
            this_U,med_U=getions_uparam(result,U_hist_bin)
            U_hist.append(this_U/wdt)
            U_medians.append(med_U)
            
            #grab the X HI histogram 
            this_X,med_X=getions_xhi(result,X_hist_bin,interp)
            X_hist.append(this_X/wdt)
            X_medians.append(med_X)
            
            #grab the X Si2 ionisation correction histogram 
            this_Xsi2,med_Xsi2=getions_xsi2(result,Xsi2_hist_bin,interp)
            Xsi2_hist.append(this_Xsi2/wdt)
            Xsi2_medians.append(med_Xsi2)
            
            #grab the X Si3 fraction histogram 
            this_Xsi3,med_Xsi3=getions_xsi3(result,Xsi3_hist_bin,interp)
            Xsi3_hist.append(this_Xsi3/wdt)
            Xsi3_medians.append(med_Xsi3)

            #grab the XC2 histogram 
            this_XC2,med_XC2=getions_xc2(result,XC2_hist_bin,interp)
            XC2_hist.append(this_XC2/wdt)
            XC2_medians.append(med_XC2)
            
            #grab the XC3 histogram 
            this_XC3,med_XC3=getions_xc3(result,XC3_hist_bin,interp)
            XC3_hist.append(this_XC3/wdt)
            XC3_medians.append(med_XC3)

            #grab the XC4 histogram 
            this_XC4,med_XC4=getions_xc4(result,XC4_hist_bin,interp)
            XC4_hist.append(this_XC4/wdt)
            XC4_medians.append(med_XC4)

    #save
    pickle.dump({'U_base':U_hist_bin,'U_hist':U_hist,'U_med':U_medians,
                 'XHI_base':X_hist_bin,'XHI_hist':X_hist,'XHI_med':X_medians,
                 'XSiII_base':Xsi2_hist_bin,'XSiII_hist':Xsi2_hist,'XSiII_med':Xsi2_medians,
                 'XSiIII_base':Xsi3_hist_bin,'XSiIII_hist':Xsi3_hist,'XSiIII_med':Xsi3_medians,
                 'XCII_base':XC2_hist_bin,'XCII_hist':XC2_hist,'XCII_med':XC2_medians,
                 'XCIII_base':XC3_hist_bin,'XCIII_hist':XC3_hist,'XCIII_med':XC3_medians,
                 'XCIV_base':XC4_hist_bin,'XCIV_hist':XC4_hist,'XCIV_med':XC4_medians,
                 'name':allname},open(modelsdir+"/"+figname+'_ionisation.pkl',"w"))
