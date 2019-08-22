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

import numpy as np
import sys
import os

from pyigm.metallicity.mcmc import mcmc_ions
from astropy.table import Table


def read_guesses_file(guesses_file, row_index_to_run):
    """
    For use with --wotta, this reads the "input guesses" file
    and sets "optim=guess" for the MCMC
    ...unless met_guess_in == "False" or
      dens_guess_in == "False" (see run_mcmc_wotta())
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
          'comment_out',
          'notes_in')
    
    fmts=(str,
        str,
        float,
        float,
        float,
        str,
        str,
        str,
        str,
        str,
        )
    
    input_dict={}
    for key in keys:
        input_dict[key] = []
    
    ##Read the input file
    infile=open(guesses_file, 'r')
    for line in infile:
        if line and not line.startswith("#"):
            ##For each line, use "|" as the delimiter. Then we'll
            ##  want to remove any whitespace surrounding the leftover
            ##  text (but not WITHIN it).
            ##  E.g., "  my comment   "  -->  "my comment"
            linespl = line.strip().split()
            # linespl = [i.strip() for i in line.strip().split('|')]
            ##
            for ls in range(len(linespl)):
                try:
                    input_dict[keys[ls]].append(fmts[ls](linespl[ls]))
                except:
                    input_dict[keys[ls]].append(linespl[ls])
    
    infile.close()


    ##We only care about the row that is given as the command line argument
    thisrow = [input_dict[key][row_index_to_run] for key in keys]
    
    ##Return this row and also everything
    ##  That way, if we want to use this function elsewhere, we can
    return thisrow, input_dict



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
    if not args.carbalphamean:
        args.carbalphamean=0.0
    if not args.carbalphasigma:
        args.carbalphasigma=0.5

    ##Re-index the row number
    row_index = int(args.row) - 1
    
    ##Read the input file and get the desired row
    [ok_to_run, \
        sightline, \
        args.met, \
        args.dens, \
        args.carbalpha, \
        carbalpha_use, \
        logUconstraint_use, \
        args.UVB, \
        comment_out, \
        notes], all_guesses = read_guesses_file(args.guessesfile, row_index)
    
    if (str(args.met).lower() == "false") or (str(args.dens).lower() == "false"):
        ##The important one
        args.optim=False
        ##Set these just in case it makes a difference
        args.met=False
        args.dens=False


    ##Take care of logU guess
    ##If the guesses_file says either "True" or "False", then:
    ##  (1) just keep the default logUmean and logUsigma
    ##  (2) set the use_logUconstraint? option to be equal to the value in the file
    ##Else, convert the guess to logUmean and logUsigma values to OVERRIDE the default
    ##  and set use_logUconstraint? option to True
    if logUconstraint_use.lower() == 'true' or logUconstraint_use.lower() == 'false':
        args.logUconstraint = str(logUconstraint_use)
    else:
        try:
            lum, lus = logUconstraint_use.split(',')
            args.logUconstraint = 'True'
            args.logUmean=float(lum)
            args.logUsigma=float(lus)
        except:
            print("")
            print("ERROR!! There was a problem determining if you wanted to use the")
            print("    logUconstraint based on the guesses file. Assuming you do NOT.")
            print("")
            args.logUconstraint = 'False'



    args.fileinput=input_file_dir+"/mcmc."+sightline+".in"        

    ##Grid file
    ##Test whether to include the carbalpha parameter, and use the correct UVB
    ##...But set a default first
    if not args.grid:
        args.grid = 'grid_minimal'
    
    if str(carbalpha_use).lower() == 'false':
        ##Do NOT use the carbalpha parameter
        ca_text = ""
        args.carbalphaconstraint = 'False'
    else:
        ##Use the carbalpha parameter
        ca_text = "_carbalpha"
        try:
            cam, cas = carbalpha_use.split(',')
            args.carbalphamean = float(cam)
            args.carbalphasigma = float(cas)
            args.carbalphaconstraint = 'True'
        except:
            args.carbalphaconstraint = 'Flat'
            # args.carbalphamean = 'Flat'
            # args.carbalphasigma = 'Flat'

    args.grid = "{}_{}{}.pkl".format(args.grid, args.UVB, ca_text)



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
    
    ##Add our own info to the saved information
    obsinfo = {}
    obsinfo['UVB'] = args.UVB
    obsinfo['carbalphaconstraint'] = args.carbalphaconstraint
    obsinfo['carbalphamean'] = args.carbalphamean
    obsinfo['carbalphasigma'] = args.carbalphasigma
    obsinfo['logUconstraint'] = args.logUconstraint
    obsinfo['logUmean'] = args.logUmean
    obsinfo['logUsigma'] = args.logUsigma



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
    run_mcmc(args, obsinfo)
    
    return


################################################
################################################


def run_mcmc(args, obsinfo=None):

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
    if not obsinfo:
        obsinfo = {}
    for ii in range(len(sets)):
        if(data['ion'][sets[ii]] == 'HI'):
            #isolate hydrogen and info
            obsinfo['NHI']=data['logn'][sets[ii]]
            obsinfo['eNHI'] = data['elogn'][sets[ii]]
            obsinfo['hiflag'] = data['flag'][sets[ii]]
            obsinfo['z'] = data['zabs'][sets[ii]]
            obsinfo['errz'] = data['ezabs'][sets[ii]]
            obsinfo['name'] = sightline
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
                   logUconstraint=args.logUconstraint, logUmean=args.logUmean, logUsigma=args.logUsigma,
                   carbalphaconstraint=args.carbalphaconstraint, carbalphamean=args.carbalphamean, carbalphasigma=args.carbalphasigma,
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
    parser.add_argument('-grid', type=str, help='Full path+basename to Cloudy grid. E.g., "/afs/crc.nd.edu/group/CGMND/Cloudy_grids/grid_cgm_extensive" (this is the default). Do NOT include UVB, carbalpha, or .pkl extension; you specify these with the other options, on a per-absorber basis.')
    parser.add_argument('-logUconstraint', type=str, help='Should we use logU constraint on density? Can be: "True", "False", or comma-separated values for logUmean and logUsigma (e.g., "-3.1,0.2"), which assumes "True". Note: This final (CSV) format overrides -logUmean and -logUsigma options!')
    parser.add_argument('-logUmean', type=str, help='If we use logUconstraint, what is the mean for the Gaussian? Note: This may be overridden by -logUconstraint option!')
    parser.add_argument('-logUsigma', type=str, help='If we use logUconstraint, what is the sigma for the Gaussian? Note: This may be overridden by -logUconstraint option!')
    parser.add_argument('-UVB', type=str, help='The UVB to use when we are using the logU constraint on density')
    parser.add_argument('-nthread', type=int, help='Number of threads')
    parser.add_argument('-nwalkers', type=int, help='Number of walkers')
    parser.add_argument('-nsamp', type=int, help='Number of samples')
    parser.add_argument('-optim', type=str, help='Optimization method')
    parser.add_argument('-dens', type=float, help='Guess at density (optim=guess); if "False", then optim=False')
    parser.add_argument('-met', type=float, help='Guess at metallicity (optim=guess); if "False", then optim=False')
    parser.add_argument('-carbalpha', type=float, help='Guess at carbalpha; if "True", uses carbalpha grid. If comma-separated values for carbalphamean and carbalphasigma (e.g., "0.1,0.4"), which assumes "True". If "False", does not use carbalpha grid')
    parser.add_argument('-carbalphamean', type=float, help='Guess at carbalpha; if "True", uses carbalpha grid. If comma-separated values for carbalphamean and carbalphasigma (e.g., "0.1,0.4"), which assumes "True". If "False", does not use carbalpha grid')
    parser.add_argument('-carbalphasigma', type=float, help='Guess at carbalpha; if "True", uses carbalpha grid. If comma-separated values for carbalphamean and carbalphasigma (e.g., "0.1,0.4"), which assumes "True". If "False", does not use carbalpha grid')
    parser.add_argument("--testing", help="Set to test (over-rides minimum nwalkers)", action="store_true")
    parser.add_argument("--wotta", help="If used, reads in files using Wotta's file format (there is a guesses file, and each sightline has separate input file). If specified, the guesses file contains: dummy column at the front (not used here); the sightline name; metallicity initial guess; density initial guess; carbon/alpha ratio (carbalpha) initial guess (required, even if you don't want to use carbalpha); whether or not to allow carbalpha to vary; whether to use the Wotta+16 logUconstraint as a prior on the density; the UVB to use; ions to comment out (not used here);  and any additional notes (not used here). Then, all that needs to be specified here is: -guessesfile=__; -row=__ (in the guessesfile, usually automatically done by the supercomputer submission script); -nthread=__ (also done by the submission script); -nwalkers=__; and -nsamp=__.", action="store_true")
    parser.add_argument('-guessesfile')
    parser.add_argument("-row", type=str, help="Row (sightline) in the guesses file to run, one-indexed (not zero-indexed)")
    pargs = parser.parse_args()


    # Run
    if pargs.wotta == True:
        run_mcmc_wotta(pargs)
    else:
        run_mcmc(pargs)



if __name__ == '__main__':
    args=sys.argv
    main(args)


# Example
# pyigm_mtlmcmc J012156.03+144823.8_z2.662 alldata.txt savehere grid_minimal.pkl -nthread=1 -nwalkers=80 -nsamp=80

##Wotta Example
# pyigm_mtlmcmc --wotta -guessesfile="MCMC_initial_guesses-run_me.dat" -row=2 -nthread=12 -nwalkers=400 -nsamp=400
