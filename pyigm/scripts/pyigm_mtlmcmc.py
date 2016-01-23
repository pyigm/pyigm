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


 
def run_mcmc(argv):

    import argparse
    from pyigm.metallicity.mcmc import mcmc_ions
    import numpy as np

    #get the call
    parser = argparse.ArgumentParser(description='Running grid on shared memory system')
    parser.add_argument('listmod')
    parser.add_argument('fileinput')
    parser.add_argument('outsave')
    parser.add_argument('grid')
    parser.add_argument('proc', type=int)
    args = parser.parse_args(argv)
 
    #read the input file
    fl = open(args.listmod)

    #read in the data file
    print 'Loading the data...'
    #load the file
    data=np.loadtxt(args.fileinput, dtype={'names':('name','zabs','ezabs','ion','logn','elogn','flag','sample'),
                                           'formats':('S50','f8','f8','S5','f5','f5','f5','S5')})
    
    #loop on systems to fit in this batch 
    for sy in fl:
     
        #find th eions for the current system
        sets=np.where(data['name'] == sy.strip())
        sets=sets[0]
               
        #now loop over ions for this system and initialise the tuples
        observ=[]
        for ii in range(len(sets)):
            if(data['ion'][sets[ii]] == 'HI'): 
                               
                #isolate hydrogen and info
                #print ii, sets[ii], data['ion'][sets[ii]], sy
                obsinfo={'NHI':data['logn'][sets[ii]],'eNHI':data['elogn'][sets[ii]],'hiflag':data['flag'][sets[ii]],
                         'z':data['zabs'][sets[ii]],'errz':data['ezabs'][sets[ii]],'name':sy.strip()}
            else:
                
                #append ion value to list 
                observ.append((data['ion'][sets[ii]],data['logn'][sets[ii]],data['elogn'][sets[ii]],data['flag'][sets[ii]]))
                
        #run the mcmc
        print 'Ready to run mcmc for {}'.format(sy)
        
        #pick optimised values for 12 processors - cosma (proc*NN walkers, proc*YY samples)
        mcmc=mcmc_ions.mcmc_ions(observ,obsinfo,args.grid,nwalkers=(args.proc*80),nsamp=(args.proc*40),
                                 optim=False,threads=args.proc,outsave=args.outsave)
        #mcmc=mcmc_ions.mcmc_ions(observ,obsinfo,args.grid,nwalkers=(10),nsamp=(50),
        #                         optim=False,threads=args.proc,outsave=args.outsave)

    print 'All done with this batch'
    
    return

if __name__ == "__main__":
    import sys
    run_mcmc(sys.argv[1:])
    
  











