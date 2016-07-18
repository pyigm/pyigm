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
 
def run_mcmc(args):

    from pyigm.metallicity.mcmc import mcmc_ions
    import numpy as np
    from astropy.table import Table

    #load the file
    sy = args.sy
    data = Table.read(args.fileinput, format='ascii')

    #find the ions for the current system
    sets = np.where(data['name'] == sy.strip())
    sets = sets[0]
               
    #now loop over ions for this system and initialise the tuples
    observ=[]
    for ii in range(len(sets)):
        if(data['ion'][sets[ii]] == 'HI'):
            #isolate hydrogen and info
            obsinfo={'NHI':data['logn'][sets[ii]],'eNHI':data['elogn'][sets[ii]],'hiflag':data['flag'][sets[ii]],
                     'z':data['zabs'][sets[ii]],'errz':data['ezabs'][sets[ii]],'name':sy.strip()}
        else:
            #append ion value to list
            observ.append((data['ion'][sets[ii]],data['logn'][sets[ii]],data['elogn'][sets[ii]],data['flag'][sets[ii]]))

    #run the mcmc
    print 'Ready to run mcmc for {}'.format(sy)
    if args.optim in ['guess', 'guess_NHI']:
        obsinfo['met'] = args.met
        obsinfo['dens'] = args.dens
        optim = args.optim
    else:
        optim = False

    #pick optimised values for 12 processors - cosma (proc*NN walkers, proc*YY samples)
    mcmc=mcmc_ions(observ,obsinfo, args.grid, nwalkers=(args.nwalkers),
                   nsamp=(args.nsamp), optim=optim, threads=args.nthread,
                   outsave=args.outsave, testing=args.testing)

    print 'All done with this batch'
    
    #return

def main(args=None):
    import argparse
    #get the call
    parser = argparse.ArgumentParser(description='Running grid on shared memory system')
    parser.add_argument('sy', type=str, help='Name of the System to analyze')
    parser.add_argument('fileinput')
    parser.add_argument('outsave')
    parser.add_argument('grid')
    parser.add_argument('-nthread', type=int, help='Number of threads')
    parser.add_argument('-nwalkers', type=int, help='Number of walkers')
    parser.add_argument('-nsamp', type=int, help='Number of samples')
    parser.add_argument('-optim', type=str, help='Optimization method')
    parser.add_argument('-dens', type=float, help='Guess at density (optim=guess)')
    parser.add_argument('-met', type=float, help='Guess at metallicity (optim=guess)')
    parser.add_argument("--testing", help="Set to test (over-rides nwalkers)", action="store_true")
    pargs = parser.parse_args()

    # Run
    run_mcmc(pargs)


# Example
# pyigm_mtlmcmc J012156.03+144823.8_z2.662 alldata.txt savehere grid_minimal.pkl -nthread=1 -nwalkers=80 -nsamp=80
