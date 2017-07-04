#!/usr/bin/env python
"""
This is the code that loop over data and runs mcmc chains.

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

import pdb
import numpy as np
import sys
import os

from pyigm.metallicity.mcmc import mcmc_ions
from astropy.table import Table


def read_guesses_file(guesses_file, row_index_to_run):
    """
    For use with --wotta, this reads the "input guesses" file
    and sets "optim=guess" for the MCMC.
    """
    
    
    ##Define the columns in the input file
    keys=('ok_to_run_in',
          'sightline_in',
          'met_guess_in',
          'dens_guess_in',
          'carbalpha_guess_in',
          'carbalpha_use_in',
          'logUconstraint_use_in',
          'uvb_in',
          'notes_in')
    
    input_dict={}
    for key in keys:
        input_dict[key] = []
    
    ##Read the input file
    infile=open(guesses_file, 'r')
    for line in infile:
        if line and not line.startswith("#"):
            linespl = line.strip().split()
            ##
            for ls in xrange(len(linespl)):
                input_dict[keys[ls]].append(linespl[ls])
    
    infile.close()


    ##We only care about the row that is given as the command line argument
    thisrow = [input_dict[key][row_index_to_run] for key in keys]
    
    ##Return the array as separate values, not as an array
    return *thisrow



################################################
################################################
################################################



def run_mcmc_wotta(args):
    """
    For use with --wotta, this:
    (1) Read initial guesses file,
    (2) Set a few options that Wotta uses, and
    (3) Calls run_mcmc() to run the MCMC code on the desired line
    
    Parameters
    ---------
    args : see all of the arguments detailed in main()
    
    """
    
    ##Note: row_to_run should be indexed from 1 not 0!! We subtract 1 below!
    
    ##CBW adjust
    input_file_dir = "./input_files/"
    args.outsave='MCMC_FULL'
    args.optim='guess'
    
    if not args.nwalkers:
        args.nwalkers=400
    if not args.nsamp:
        args.nsamp=400
    if not args.logUmean:
        args.logUmean=-2.968
    if not args.logUsigma:
        args.logUsigma=0.481
    
    ##Re-index the row number
    row_index = int(args.row) - 1
    
    ##Read the input file and get the desired row
    ok_to_run, \
        sightline, \
        args.met, \
        args.dens, \
        args.carbalpha, \
        carbalpha_use, \
        args.logUconstraint_use, \
        args.UVB, \
        notes = read_guesses_file(args.guessesfile, row_index)


    args.fileinput=input_file_dir+"/mcmc."+sightline+".in"        

    ##Grid file
    ##Test whether to include the carbalpha parameter
    if str(carbalpha_use).lower() == 'false':
        ##Do NOT use the carbalpha parameter
        args.grid = os.getcwd()+'/../Cloudy_grids_'+uvb+'/grid_minimal.pkl'
    else:
        ##Use the carbalpha parameter
        args.grid = os.getcwd()+'/../Cloudy_grids_'+uvb+'/grid_minimal_carbalpha.pkl'



    ########################
    ## Read sightline's data
    ########################
    
    ##Get NHI (and the number of systems)
    data = Table.read(args.fileinput, format='ascii', comment="#")
    HI_rows = np.where(data['ion'] == 'HI')[0]
    if len(HI_rows) > 1:
        print("Too many 'HI'. You can only run one system.")
        sys.exit()
    
    
    args.sightline=data['name'][HI_rows[0]].strip()



    ########################
    ##Check to make sure this row was marked to run
    ########################
    
    # if ok_to_run == 'Y':
    #     run_mcmc(args)
    # 
    # else:
    #     print(sightline+" was NOT marked to run!")
    #     
    #     ##Can't be interactive when running on the CRC
    #     # run_anyway = raw_input("Do you want to run anyway (y/n)? ")
    #     # if run_anyway.lower() == 'y':
    #     #     run_mcmc(args)
    
    ##Always run it. It would stink to sit in the CRC queue and
    ##  then fail as soon as it starts because the user forgot
    ##  to change the stinking "Run me?" column.
    run_mcmc(args)
    
    return


################################################
################################################


def run_mcmc(args):

    """ 
    (1) Read observational column density data and
    (2) Run MCMC

    Parameters
    ---------
    args : see all of the arguments detailed in main()

    """

    #load the file
    sightline = args.sightline.strip()
    data = Table.read(args.fileinput, format='ascii', comment='#')
    
    #find the ions for the current system
    sets = np.where(data['name'] == sightline)[0]
    
    #now loop over ions for this system and initialise the tuples
    observ=[]
    for ii in range(len(sets)):
        if(data['ion'][sets[ii]] == 'HI'):
            #isolate hydrogen and info
            obsinfo={'NHI':data['logn'][sets[ii]],'eNHI':data['elogn'][sets[ii]],'hiflag':data['flag'][sets[ii]],
                     'z':data['zabs'][sets[ii]],'errz':data['ezabs'][sets[ii]],'name':sightline}
        else:
            #append ion value to list
            observ.append((data['ion'][sets[ii]],data['logn'][sets[ii]],data['elogn'][sets[ii]],data['flag'][sets[ii]]))
    

    ########################
    ## Run the MCMC
    ########################
    
    ##Should we optimize?
    print('Ready to run mcmc for {}'.format(sightline))
    if args.optim in ['guess', 'guess_NHI']:
        obsinfo['met'] = args.met
        obsinfo['dens'] = args.dens
        obsinfo['carbalpha'] = args.carbalpha
        optim = args.optim
    else:
        optim = False
    
    #pick optimised values for 12 processors - cosma (proc*NN walkers, proc*YY samples)
    mcmc=mcmc_ions(observ,obsinfo, args.grid, 
                   logUconstraint=args.logUconstraint, logUmean=args.logUmean, logUsigma=args.logUsigma
                   UVB=args.UVB,
                   nwalkers=(args.nwalkers),
                   nsamp=(args.nsamp), optim=optim, threads=args.nthread,
                   outsave=args.outsave, testing=args.testing)

    print('All done with this batch')
    
    #return


################################################
################################################


def main(args=None):
    import argparse
    #get the call
    parser = argparse.ArgumentParser(description='Running grid on shared memory system')
    parser.add_argument('-sightline', type=str, help='Name of the System to analyze')
    parser.add_argument('-fileinput')
    parser.add_argument('-outsave')
    parser.add_argument('-grid')
    parser.add_argument('-logUconstraint', type=str, help='Should we use logU constraint on density')
    parser.add_argument('-logUmean', type=str, help='If we use logUconstraint, what is the mean for the Gaussian?')
    parser.add_argument('-logUsigma', type=str, help='If we use logUconstraint, what is the sigma for the Gaussian?')
    parser.add_argument('-UVB', type=str, help='The UVB to use when we are using the logU constraint on density')
    parser.add_argument('-nthread', type=int, help='Number of threads')
    parser.add_argument('-nwalkers', type=int, help='Number of walkers')
    parser.add_argument('-nsamp', type=int, help='Number of samples')
    parser.add_argument('-optim', type=str, help='Optimization method')
    parser.add_argument('-dens', type=float, help='Guess at density (optim=guess)')
    parser.add_argument('-met', type=float, help='Guess at metallicity (optim=guess)')
    parser.add_argument('-carbalpha', type=float, help='Guess at carbalpha (optim=guess)')
    parser.add_argument("--testing", help="Set to test (over-rides minimum nwalkers)", action="store_true")
    parser.add_argument("--wotta", help="Read in files using Wotta's file format (there is a guesses file, and each sightline has separate input file). If specified, the guesses file contains: 'run now?' column; the sightline name; metallicity initial guess; density initial guess; carbon/alpha ratio (carbalpha) initial guess; whether or not to allow carbalpha to vary; whether to use the Wotta+16 logUconstraint as a prior on the density; the UVB to use; and any notes (these are not used). Then, all that needs to be specified here is: -guessesfile=__; -row=__ (in the guessesfile); -nthread=__; -nwalkers=__; and -nsamp=__.", action="store_true")
    parser.add_argument('-guessesfile')
    parser.add_argument("-row", type=str, help="Row (sightline) in the guesses file to run, one-indexed (not zero-indexed)")
    pargs = parser.parse_args()


    # Run
    if pargs.wotta == True:
        run_mcmc_wotta(pargs)
    else:
        run_mcmc(pargs)


# Example
# pyigm_mtlmcmc J012156.03+144823.8_z2.662 alldata.txt savehere grid_minimal.pkl -nthread=1 -nwalkers=80 -nsamp=80

##Wotta Example
# pyigm_mtlmcmc --wotta -guessesfile="MCMC_initial_guesses-run_me.dat" -row=2 -nthread=12 -nwalkers=400 -nsamp=400
