"""
Given a grid of cloudy models and a list of measured ions, run MCMC to find
best fit values for free parameters in the grid model

data        -> a list of tuples containing column densities and associated errors for all the
               observed ions

               [(ion, log column, log error, flag), ...]

               Upper limits are flagged by  -1
               Lower limits are flagged by  -2

               For limits, error contains the typical rms of the spectrum at the wave of the line

infodata    -> information about the data (NHI, z, etc.. used to define the priors) in a dictionary
               [logNHI,errlogNHI,hiflag,redshift,errorredshift,outputname]

               HI top-hat prior is flagged by -3

models      -> grid of cloudy models to be used

nwalkers    -> The number of walkers in the MCMC run. [500]
nsamp       -> The number of samples in the MCMC chains. [250]

threads     -> number of processors

outsave     -> where to place output stuff

optim       -> Initialise with optimisation rather than optimising the likelihood for starting ball
               local: find the max likelihood given the intial guess at grid centre (typically fails)
               global: search for global max in likelihood (maybe a good choice)
               false: no optimisation - start from random position in the grid

effnhi      -> If True, use the NHI in the grid to compute the prior rather than the input value in cloudy.
               This is useful for models in which the input NHI does not correspond to the output one
               (e.g. CIE with high temperatures or in presence of molecular phase)

Written by Michele Fumagalli in Durham, Summer 2015
michele.fumagalli@durham.ac.uk

Copyright (C) 2015 Michele Fumagalli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function, absolute_import, division, unicode_literals

#start with some general imports

import matplotlib as mpl

#The next line works on queue systems for DISPLAY issues
mpl.use('Agg')

import warnings
import pdb
import copy

#Here some general import
import matplotlib.pyplot as plt

try:
    import corner
except ImportError:
    raise ImportError("Install corner to use metallicity.mcmc")

import pickle

#Specify decimal places for high precision arithmetic
try:
    import mpmath as mmath
except ImportError:
    raise ImportError("Install mpmath to use metallicity.mcmc")
else:
    mmath.mp.dps = 200

try:
    import emcee
except ImportError:
    raise ImportError("Install emcee to use metallicity.mcmc")

import numpy as np
import os
#enable warning handling in numpy
#np.seterr(all='raise')

from scipy import interpolate
#enable warning handling in scipy
#sp.seterr(all='raise')
import scipy.optimize as op

from linetools import utils as ltu

class Emceebones(object):
    """
    This is a class that does all the fun stuff like bookeeping, plots, driver for emcee etc..
    """

    def __init__(self,data,infodata,model,nwalkers,nsamp,threads,outsave,optim,effnhi):
        """ First, do some bookeeping like finding ions and preparing the grid
        Parameters
        ----------
        data : list of tuples
          observations
        infodata :
        model :
        nwalkers : int
          Number of walkers
        nsamp : int
          Number of samples
        threads : int
          Number of threads (for multi-processing)
        outsave
        optim : str
        effnhi

        Returns
        -------

        """
        print("Initialise class mcmc...")

        #save some variables
        self.nwalkers=nwalkers
        self.nsamp=nsamp
        self.threads=threads
        self.info=infodata
        self.outsave=outsave
        self.paramguess=0
        self.final=0
        self.optim=optim      #flag to switch type of optimisation
        self.effnhi=effnhi    #flag to use effective NHI

        #load the grid of models to memory
        self.loadmodel(data,model)

        #At this point self.data, self.mod_ions, self.colm contain only what you need
        #Initialise the interpolator next

    def loadmodel(self, data, model):
        """ Utility function that loads the model

        Parameters
        ----------
        data :
          observations
        model : pickle file

        Returns
        -------

        """


        #load the model in a format that can be handled later on
        fil=open(model)
        modl=pickle.load(fil)
        fil.close()

        #unpack axis tag, axis value, grid column, grid ions
        self.mod_axistag=modl[0]
        self.mod_axisval=[]

        #define the dimension of the problem
        self.ndim=len(self.mod_axistag)
        self.nmodels=1
        for tt in self.mod_axistag:
            self.nmodels=self.nmodels*(modl[1][tt]).size
            #append axis value in a list
            self.mod_axisval.append(modl[1][tt])

        print("The problem as dimension {}".format(self.ndim))
        print("There are {} models".format(self.nmodels))


        #now queue up the model and data
        self.nions=0
        self.data=[]          #obs with corresponding ions in model
        self.mod_colm=[]      #columns for the ions in model matched to obs
        self.mod_colm_tag=[]  #ions in obs & model

        #loop over all the ions observed
        for obs in data:

            # Check for zero error (causes unexepcted failure)
            if obs[2] <= 0.:
                raise ValueError("Cannot have 0 error on the column density, even in a limit.  Fix {}".format(obs[0]))

            #check obs one by one for corresponding entries in model
            if modl[2].has_key(obs[0]):

                #stack the observable
                #[list of tuples for observations with (ion,column,error)]
                self.data.append(obs)
                #stack the columns for the models
                #[list of ndim arrays of columns for each ion]
                self.mod_colm.append(modl[3]['N_'+obs[0]])
                self.mod_colm_tag.append(obs[0])
                self.nions=self.nions+1

            else:

                print("{} not in the grid: skip this ion!".format(obs[0]))

        print("Handling {} ions".format(self.nions))

        #at last store the NHI grid if useful for effective NHI
        if(self.effnhi):
            print('Using effective NHI.. Load values')
            self.nhigrid=modl[3]['N_HI']

        return

    def plotinfo(self, sampler, use_pkl=False):
        """ Function that makes some useful plots of stuff

        Parameters
        ----------
        sampler

        Returns
        -------

        """
        print('Writing outputs...')

        #pickle the results to disk
        if use_pkl:
            wout=open(self.outsave+'/'+self.info['name']+'_emcee.pkl','w')
            pickle.dump(self.final,wout)
            wout.close()
        else:  # hd5
            import h5py
            import json
            with h5py.File(self.outsave+'/'+self.info['name']+'_emcee.hd5', 'w') as f:
                # Input
                in_group = f.create_group('inputs')
                for in_key in ['data', 'ions', 'guess']:
                    in_group[in_key] = self.final[in_key]
                for key in self.final['info']:
                    in_group.attrs[key] = self.final['info'][key]
                # Output
                out_group = f.create_group('outputs')
                mcmc_dict = dict(nwalkers=self.nwalkers, nsamp=self.nsamp,
                                 nburn=self.burn, nthread=self.threads)
                out_group.attrs['mcmc'] = unicode(json.dumps(mcmc_dict))
                for out_key in ['tags', 'results', 'pdfs', 'best_fit',
                                'effNHI', 'acceptance']:
                    out_group[out_key] = self.final[out_key]

        #Start by plotting the chains with initial guess and final values
        fig=plt.figure()
        xaxis=np.arange(0,self.nsamp,1)

        #Loop over each parameter and plot chain
        for ii in range(self.ndim):

            #add chains for this param
            thischain=sampler.chain[:,:,ii]
            ax = fig.add_subplot(self.ndim, 1, ii+1)
            for jj in range(self.nwalkers):
                ax.plot(xaxis,thischain[jj,:],color='grey')
            ax.set_ylabel(self.mod_axistag[ii])

            #overplot burnt in cut
            ax.axvline(x=self.burn,color='red',linewidth=3)
            #overplot median
            ax.axhline(y=np.median(self.final['pdfs'][:,ii]),color='blue',linewidth=3)
            #overplot starting point
            ax.axhline(y=self.paramguess[ii],color='black',linewidth=3)

        #Last plot
        ax.set_xlabel('Steps')
        fig.savefig(self.outsave+'/'+self.info['name']+'_chains.pdf')

        #now do a corner plot
        samples = sampler.chain[:,self.burn:, :].reshape((-1,self.ndim))
        cfig = corner.corner(samples, labels=self.mod_axistag, quantiles=[0.05,0.5,0.95],verbose=False)
        cfig.savefig(self.outsave+'/'+self.info['name']+'_corner.pdf')

        #now plot the residuals
        rfig=plt.figure()

        #plot values
        xaxis=np.arange(0,self.nions,1)
        axlab=[]

        plt.scatter(xaxis,self.final['best_fit'], color='blue',alpha=0.2,s=100)

        for ii in range(self.nions):

            axlab.append(self.data[ii][0])

            if(self.data[ii][3] == -1):
                #upper limit
                plt.scatter([ii],[self.data[ii][1]], color='red',s=100,marker='v')
            elif(self.data[ii][3] == -2):
                #lower limit
                plt.scatter([ii],[self.data[ii][1]], color='red',s=100,marker='^')
            else:
                plt.errorbar([ii],[self.data[ii][1]],xerr=0,yerr=[self.data[ii][2]],color='red',marker='.')


        plt.xticks(xaxis,axlab,rotation='vertical')
        plt.ylabel('Log Column')
        rfig.savefig(self.outsave+'/'+self.info['name']+'_residual.pdf')

        print('All done with system {}!'.format(self.info['name']))

    def setup_emc(self):
        """ As named
        Returns
        -------

        """
        emc=Emceeutils()

        #simple pass of utility variables to emcee
        emc.ndim=self.ndim
        emc.nmodels=self.nmodels
        emc.nions=self.nions
        emc.nwalkers=self.nwalkers
        emc.nsamp=self.nsamp
        emc.threads=self.threads
        emc.data=self.data
        emc.info=self.info
        emc.mod_colm=self.mod_colm
        emc.mod_colm_tag=self.mod_colm_tag
        emc.mod_axistag=self.mod_axistag
        emc.mod_axisval=self.mod_axisval
        emc.effnhi=self.effnhi
        #
        if(emc.effnhi):
            emc.nhigrid=self.nhigrid
        # Save internally
        self.emc = emc
        #
        return

    def __call__(self):

        print('Preparing emcee...')

        #now init the emcee utlility class
        self.setup_emc()
        emc = self.emc  # For backwards compatability

        #if in use, pass the HI grid to compute effective NHI
        #now initialise the interpolator
        emc.init_interpolators()

        #DEBUG sequence: check the interpolator (Works ok on 8 May 2015)
        #emc.test_interpolators()

        #DEBUG sequence: check the prior+likelihood
        ##('col', 'red', 'met', 'dens')
        #param=[17.7,3.226,-2.0,-2]
        #zzz=np.arange(17.01,18.5,0.01)
        #self.info['hiflag']=0
        #self.info['eNHI']=0.5
        #val=[]
        #for z in zzz:
        #    param[0]=z
        #    getv=emc(param)
        #    val.append(getv)
        #    print z, getv
        #
        #plt.plot(zzz,val)
        #plt.show()


        #DEBUG sequence: init fake data
        #param=[self.info['NHI'],self.info['z'],0.82,-3.12]
        #for cc in range(self.nions):
        #   print  "('{}',{},0.05,0),".format(self.data[cc][0],emc.interpol[cc](param)[0])
        #exit()

        #now explore the likelihood for initial guess to MCMC
        #Given the use of flat priors, use prior * likelihood to
        #settle the redshift and columns in the right spot
        print("Finding initial guess...")

        #set up a -like function for min finding
        nll = lambda *args: -1.*emc(*args)

        #Pick starting value in middle of grid as first step
        initguess=[]
        for pp in range(self.ndim):
            tag=self.mod_axistag[pp]
            if tag == 'col':
                #Handle special case column density
                initguess.append(self.info['NHI'])
            elif tag == 'red':
                #Handle special case redshift
                initguess.append(self.info['z'])
            else:
                #Handle all other paramaters
                if self.optim in ['guess', 'guess_NHI']:
                    if tag == 'met':
                        initguess.append(self.info['met'])
                    elif tag == 'dens':
                        initguess.append(self.info['dens'])
                    else:
                        raise ValueError("Not ready for this one")
                else:
                    initguess.append(np.median(self.mod_axisval[pp]))

        #now go through a bunch of cases in which the starting ball is
        #set according to local, global, or no optimisation
        if self.optim == 'local':
            raise NotImplementedError("The following is not well implemented")
            #optimise with local mimimum search
            print('Searching local minimum [poor choice!]...')
            mcmcguess = op.minimize(nll, initguess)
            print('Optimisation completed? {}'.format(mcmcguess["success"]))
            #now prepare the starting ball
            if mcmcguess["success"]:
                #if this worked ok, go ahead
                self.paramguess = mcmcguess["x"]
                #set up the starting ball at this position
                pos = [self.paramguess + 1e-2*np.random.randn(self.ndim) for i in range(self.nwalkers)]
                print('Now starting emcee ...')
            else:
                #if not, reset to random
                print('Optimisation failed: starting emcee from random positions')
                self.optim=False

        if self.optim == 'global':
            raise NotImplementedError("The following is not well implemented")
            #optimise with global mimimum search
            print('Searching global minimum [takes a bit]...')
            mcmcguess = op.basinhopping(nll,initguess)
            print('Optimisation completed? {}'.format(mcmcguess["message"][0]))
            #set up the starting ball at this position
            if 'successfully' in mcmcguess["message"][0]:
                #if this worked ok, go ahead
                self.paramguess = mcmcguess["x"]
                #set up the starting ball at this position pos[nwalkers,ndim]
                pos = [self.paramguess + 1e-2*np.random.randn(self.ndim) for i in range(self.nwalkers)]
                print('Initial guess for {} are {}'.format(self.mod_axistag,self.paramguess))
                print('Now starting emcee...')
            else:
                #if not, reset to random
                print('Optimisation failed: starting emcee from random position')
                self.optim=False

        if self.optim in ['guess', 'guess_NHI']:
            print("Skip optimisation... Initialise at input guess")
            #set up a ball centred at initguess
            pos = [initguess + 1e-2*np.random.randn(self.ndim) for i in range(self.nwalkers)]
            self.paramguess = initguess
            if self.optim == 'guess_NHI':
                # Don't center on NHI
                ballNHI=np.random.random(self.nwalkers)*(2*self.info['eNHI'])+(self.info['NHI']-self.info['eNHI'])
                for pp in range(self.ndim):
                    tag=self.mod_axistag[pp]
                    if tag == 'col':
                        #For column density, randomise within the error range
                        for i in range(self.nwalkers):
                            pos[i][pp] = ballNHI[i]
                    if tag == 'met':
                        for i in range(self.nwalkers):
                            pos[i][pp] -= (ballNHI[i]-self.info['NHI'])
                    else:
                        pass

        if self.optim == False:
            print("Skip optimisation... Initialise at random with the grid")

            #set up a ball centred at initguess
            pos = [initguess + 1e-2*np.random.randn(self.ndim) for i in range(self.nwalkers)]
            self.paramguess = initguess

            #now loop on each axis of the model grid to tweak
            #the starting ball initialisation as appropriate
            for pp in range(self.ndim):

                tag=self.mod_axistag[pp]
                if tag == 'col':
                    #For column density, randomise uniformly within the error bars (+/- 1 sigma)
                    ballst=np.random.random(self.nwalkers)*(2*self.info['eNHI'])+(self.info['NHI']-self.info['eNHI'])
                    #now assign value to starting balls
                    for i in range(self.nwalkers):
                        pos[i][pp]=ballst[i]

                elif tag == 'red':
                    #For redshift, randomise within the error bars (+/- 1 sigma)
                    ballst=np.random.random(self.nwalkers)*(2*self.info['errz'])+(self.info['z']-self.info['errz'])
                    #now assign value to starting balls
                    for i in range(self.nwalkers):
                        pos[i][pp]=ballst[i]

                else:
                    #For the other parameters, randomise inside the full extent of the grid
                    ballst=np.random.random(self.nwalkers)*\
                        (np.max(self.mod_axisval[pp])-np.min(self.mod_axisval[pp]))+np.min(self.mod_axisval[pp])
                    #now assign value to starting balls
                    for i in range(self.nwalkers):
                        pos[i][pp]=ballst[i]

        print('Running {} chains for {} steps on {} processors'.format(self.nwalkers,self.nsamp,self.threads))

        #init the sampler and run
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, emc, threads=self.threads)
        sampler.run_mcmc(pos, self.nsamp)

        #Check the acceptance fraction
        accept=np.mean(sampler.acceptance_fraction)
        print("Mean acceptance fraction: {0:.3f}".format(accept))
        print("[good range 0.25-0.5. Low is bad!]")

        #remove burnt in to generate "clean" PDFs
        self.burn = 45
        if self.burn >= self.nsamp:
            raise ValueError("Burn out exceeds number of samples!")

        #get nsamples * ndim pdfs
        samples = sampler.chain[:, self.burn:, :].reshape((-1,self.ndim))

        #now compute percentiles for each paramater
        percentile=[10,16,25,50,75,84,90]
        results = np.percentile(samples,percentile,axis=0)
        medians = results[3,:]

        #print outputs
        print('---------------------')
        print('Best fit values [median,10,90%]')
        for pp in range(self.ndim):
            print(self.mod_axistag[pp], results[3,pp], results[0,pp], results[6,pp])
        print('---------------------')

        #extract the "best fit" model
        print("Extract best fit model...")
        self.best_fit=[]
        for ion in range(self.nions):
            self.best_fit.append(emc.interpol[ion](medians)[0])

        #grab effective NHI if in use
        if(self.effnhi):
            best_fit_nhi=emc.hiinterp(medians)[0]
        else:
            best_fit_nhi=False

        #find the peaks in probability space
        self.final={'ions':self.mod_colm_tag,'percent':percentile,'results':results,'acceptance':accept,
                    'pdfs':samples,'guess':self.paramguess,'tags':self.mod_axistag,'best_fit':self.best_fit,
                    'data':self.data,'info':self.info,'effNHI':best_fit_nhi}

        #now plot some information and write output
        self.plotinfo(sampler)

        return

class Emceeutils():
    """
    This is a class that is called at each iteration to evaluate priors etc...
    """

    def __init__(self):

        print("Initialise Emceeutils...")

        self.ndim=0
        self.nmodels=0
        self.nions=0
        self.nwalkers=0
        self.nsamp=0
        self.threads=0
        self.data=0
        self.info=0
        self.mod_colm=0
        self.mod_colm_tag=0
        self.mod_axistag=0
        self.mod_axisval=0
        self.interpol=[]
        self.hiinterp=0
        self.effnhi=False
        self.nhigrid=0

        return

    def lnprior(self, param):
        """ This is the function that computes the prior. It assuems non informative priors
        on paramaters, plus Gaussians around measured NHI and redshift.
        If flag -3 is raised in HI column density, then assume a flat prior within that range

        Values off the grid are suppressed by -inf

        Parameters
        ----------
        param :

        """

        #Check if inside the grid
        inside = True
        ii=0
        for pp in param:
            condition=(np.min(self.mod_axisval[ii]) <= pp <= np.max(self.mod_axisval[ii]))
            inside = inside & condition
            ii=ii+1

        #now eval the prior if inside the grid
        if (inside):
            #If passed inside test, add the values of log for Gaussian priors on red and column density
            #Eval product of gaussians
            #Also handle the special case of a top-hat prior in HI if needed.

            #start with redshift
            zindex=self.mod_axistag.index('red')
            pri_red=-1*np.log(np.sqrt(2*np.pi)*self.info['errz'])-\
                (self.info['z']-param[zindex])**2/(2*self.info['errz']**2)

            #now compute prior with column density
            nhiindex=self.mod_axistag.index('col')

            #establish the mapping between NHI in the grid and input NHI
            if(self.effnhi):
                #find current NHI given grid
                currentnhi=self.hiinterp(param)
            else:
                #leave NHI as in input
                currentnhi=param[nhiindex]

            if(self.info['hiflag'] < -2):
                #handle the tophat function
                tophat=(self.info['NHI']-self.info['eNHI'] <= currentnhi <= self.info['NHI']+self.info['eNHI'])
                if(tophat):
                    pri_nhi=0
                else:
                    pri_nhi=-np.inf
            else:
                #gaussian priors
                pri_nhi=-1*np.log(np.sqrt(2*np.pi)*self.info['eNHI'])-(self.info['NHI']-currentnhi)**2/(2*self.info['eNHI']**2)

            #add together [log product]
            prior=pri_red+pri_nhi

        else:
            #if outside return -inf
            prior=-np.inf

        return prior


    def lnlhood(self, param):
        """ This is the function that evaluates the likelihood at each point in NDIM space

        Parameters
        ----------
        param :

        Returns
        -------

        """
        likeit=0

        #loop over all the ions
        for ii in range(self.nions):

            #parity check : make sure data and models are aligned
            if(self.data[ii][0] != self.mod_colm_tag[ii]):
                raise ValueError('Mismtach between observables and models. This is a big mistake!!!')

            #now call the interpolator for models given current ion
            mod_columns=self.interpol[ii](param)

            #check if upper limit
            if(self.data[ii][3] == -1):
                #integrate the upper limit of a Gaussian - cumulative distribution
                arg=((self.data[ii][1]-mod_columns)/(np.sqrt(2)*self.data[ii][2]))[0]
                thislike=mmath.log(0.5+0.5*mmath.erf(arg))
                likeit=likeit+float(thislike)
                #print self.data[ii][0], float(thislike), self.data[ii][1], mod_columns

            #check if lower limit
            elif(self.data[ii][3] == -2):

                #integrate the lower limit of a Gaussian - Q function
                arg=((self.data[ii][1]-mod_columns)/(np.sqrt(2)*self.data[ii][2]))[0]
                thislike=mmath.log(0.5-0.5*mmath.erf(arg))
                likeit=likeit+float(thislike)
                #print self.data[ii][0], float(thislike), self.data[ii][1], mod_columns

            #if value, just eval Gaussian
            else:

                #add the likelihood for this ion
                thislike=-1*np.log(np.sqrt(2*np.pi)*self.data[ii][2])-(self.data[ii][1]-mod_columns)**2/(2*self.data[ii][2]**2)
                likeit=likeit+thislike

        return likeit


    def init_interpolators(self):
        """ This seeds the linear interpolator. For each parameter in the model,
        perform linear interpolation along each dimension of the column for one ion
        to return the column density value at an input position

        This is using scipy.interpolate.RegularGridInterpolator
        so it is assumend that the points are on a regualr grid
        While this can be changed, computing a regular rather than a sparse
        grid speed up the computation quite significantly

        A list of interpolators - one per ion - is returned
        """

        #loop over all ions
        for ii in range(self.nions):

            print('Seeding the interpolator for {}'.format(self.data[ii][0]))

            #parity check: make sure the ordering of the data and column is the same
            if(self.data[ii][0] != self.mod_colm_tag[ii]):
                raise ValueError('Mismtach between observables and models. This is a big mistake!!!')

            #stack up the interpolator for current ions
            self.interpol.append(interpolate.RegularGridInterpolator(self.mod_axisval,self.mod_colm[ii],
                                                                     method='linear',bounds_error=False,fill_value=-np.inf))

        #if using effective NHI, init the HI interpolator
        if(self.effnhi):
            self.hiinterp=interpolate.RegularGridInterpolator(self.mod_axisval,self.nhigrid,
                                                              method='linear',bounds_error=False,fill_value=-np.inf)


            #START Quick debug sequence for interpolator [works with dust grid]
            #nhitest=np.arange(17.2,20.5,0.1)
            #denstest=np.arange(-4,0.0,0.1)
            #redtest=0.3
            #mettest=0.8
            #fstartest=1.4
            #intgrd=np.zeros((nhitest.size,denstest.size))
            #fltgrd=np.zeros((nhitest.size,denstest.size))
            #
            #for ttdd in range(denstest.size):
            #    for tthh in range(nhitest.size):
            #        intgrd[tthh,ttdd]=self.hiinterp([nhitest[tthh],redtest,mettest,denstest[ttdd],fstartest])
            #        fltgrd[tthh,ttdd]=nhitest[tthh]
            #
            #plt.imshow(intgrd)
            ##plt.imshow(intgrd-fltgrd)
            #plt.colorbar()
            #plt.show()
            #exit()
            #END of debug sequence

        return


    def test_interpolators(self):
        """
        This is utility function to test the interpolation procedures

        Differently from the rest of the code - this is not general - works
        only for a specific grid (minimal).

        """

        int_dens=np.arange(-4.0,0.2,0.1)
        int_colm=np.arange(17.0,20.6,0.1)
        int_red=np.arange(0.0,4.6,0.1)
        int_met=np.arange(-4.0,1.2,0.1)

        #loop over ions
        for ii in range(self.nions):

            print('Check interpolator')


            #plane 1
            original_1=self.mod_colm[ii][:,:,0,0]
            imgplt=plt.imshow(original_1,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_orig_a_{}.png'.format(ii))
            plt.clf()

            interp_1=np.zeros((int_colm.size,int_red.size))
            xx=0
            for cc in int_colm:
                yy=0
                for rr in int_red:
                    interp_1[xx,yy]=self.interpol[ii]([cc,rr,self.mod_axisval[2][0],self.mod_axisval[3][0]])
                    yy=yy+1
                xx=xx+1
            imgplt=plt.imshow(interp_1,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_inter_a_{}.png'.format(ii))
            plt.clf()


            #plane 2
            original_2=self.mod_colm[ii][0,:,:,0]
            imgplt=plt.imshow(original_2,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_orig_b_{}.png'.format(ii))
            plt.clf()


            interp_2=np.zeros((int_red.size,int_met.size))
            xx=0
            for cc in int_red:
                yy=0
                for rr in int_met:
                    interp_2[xx,yy]=self.interpol[ii]([self.mod_axisval[0][0],cc,rr,self.mod_axisval[3][0]])
                    yy=yy+1
                xx=xx+1
            imgplt=plt.imshow(interp_2,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_inter_b_{}.png'.format(ii))
            plt.clf()

            #plane 3
            original_3=self.mod_colm[ii][0,0,:,:]
            imgplt=plt.imshow(original_3,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_orig_c_{}.png'.format(ii))
            plt.clf()


            interp_3=np.zeros((int_met.size,int_dens.size))
            xx=0
            for cc in int_met:
                yy=0
                for rr in int_dens:
                    interp_3[xx,yy]=self.interpol[ii]([self.mod_axisval[0][0],self.mod_axisval[1][0],cc,rr])
                    yy=yy+1
                xx=xx+1
            imgplt=plt.imshow(interp_3,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_inter_c_{}.png'.format(ii))
            plt.clf()

            #plane 4
            original_4=self.mod_colm[ii][:,0,0,:]
            imgplt=plt.imshow(original_4,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_orig_d_{}.png'.format(ii))
            plt.clf()

            interp_4=np.zeros((int_colm.size,int_dens.size))
            xx=0
            for cc in int_colm:
                yy=0
                for rr in int_dens:
                    interp_4[xx,yy]=self.interpol[ii]([cc,self.mod_axisval[1][0],self.mod_axisval[2][0],rr])
                    yy=yy+1
                xx=xx+1
            imgplt=plt.imshow(interp_4,origin='lower')
            imgplt.set_clim(7.0,12.0)
            plt.savefig('debug_inter_d_{}.png'.format(ii))
            plt.clf()

        return

    def __call__(self, param):
        """ This is the method that hides the computation of the likelihood/prior behind a callable function

        Parameters
        ----------
        param

        Returns
        -------

        """
        #call the likelihood function [log]
        lh = self.lnlhood(param)

        #call the prior function [log]
        lp = self.lnprior(param)

        #return likelihood*prior [lor] or inf
        if not np.isfinite(lh) or not np.isfinite(lp):
            return -np.inf
        else:
            return lh + lp


def mcmc_ions(data,infodata,model,nwalkers=500,nsamp=250,threads=12,
              outsave='emceeout', optim=False, effnhi=True, testing=False):
    """ This is the main, which does not do much

    Parameters
    ----------
    data : list
      List of tuples containing the column density constraints
    infodata : dict
      dict containing observation info, e.g. NHI, z
    model : str
      Name of the grid
    nwalkers : int, optional
    nsamp : int, optional
    threads : int, optional
    outsave : str, optional
    optim : str, optional
      Optimize?  [Not recommended]
    effnhi : bool, optional
    testing : bool, optional

    Returns
    -------

    """

    print('Getting started...')

    # Do not run with less than 100 samples or 400 walkers
    if not testing:
        nsamp=np.max([nsamp,100])
        nwalkers=np.max([nwalkers,400])

    # initialise the mcmc for this problem
    # load the observations and the model
    mcmc=Emceebones(data,infodata,model,nwalkers,
                    nsamp,threads,outsave,optim,effnhi)

    # make space for output if needed
    if not os.path.isdir(outsave):
        os.makedirs(outsave)

    # Run emcee on this dataset
    results=mcmc()


    return









