"""
Module for analyzing metallicity with MCMC

"""


import numpy as np
import sys
import os

from astropy.table import Table
from pyigm.metallicity.mcmc import mcmc_ions


################################################
################################################


def run_full_mcmc(infile='MgIIsurvey_MTL.ascii', nthread=4, met_guess='-1.0', dens_guess='-2.5', carbalpha_guess='0.0', carbalpha_use=True, logUconstraint=False,UVB='HM05'):
    
    """ Run MCMC

    Parameters
    ---------
    nthread : int, optional
      Number of threads per system

    """
    
    ########################    
    ##CBW adjust parameters
    ########################
    
    ##NOTE: In mcmc.py, we have the following lines:
    ##  > if not testing:
    ##  >     nwalkers=np.max([nwalkers,400])
    ##  >     nsamp=np.max([nsamp,100])
    ##which means the options here may be overridden later
    ##(depending on your choices here).
    
    nwalkers = int(400) ##make sure it's an integer
    nsteps = int(400) ##make sure it's an integer
    testing = False
    optim = 'guess' ##or False
    outsave = 'MCMC_FULL'
    
    ##Grid file
    ##Test whether to include the carbalpha parameter
    if str(carbalpha_use).lower() == 'false':
        ##Do NOT use the carbalpha parameter
        grid_fil = os.getcwd()+'/../Cloudy_grids_'+UVB+'/grid_minimal.pkl'
    else:
        ##Use the carbalpha parameter
        grid_fil = os.getcwd()+'/../Cloudy_grids_'+UVB+'/grid_minimal_carbalpha.pkl'
    
    
    ########################
    ## Read data
    ########################
    
    ##Get NHI (and the number of systems)
    mcmc_fil = infile
    data = Table.read(mcmc_fil, format='ascii', comment="#") ##CBW edited to add: comment="#"
    HI_rows = np.where(data['ion'] == 'HI')[0]
    if len(HI_rows) > 1:
        print("Too many 'HI'. You can only run one system.")
        sys.exit()
    
    
    ##Load the file
    data = Table.read(mcmc_fil, format='ascii', comment="#")
    
    
    ############
    ##System name
    ############
    sy=data['name'][HI_rows[0]].strip()
    
    #find the ions for the current system
    sets = np.where(data['name'] == sy)[0]
    
    #now loop over ions for this system and initialise the tuples
    observ=[]
    for ii in range(len(sets)):
        if(data['ion'][sets[ii]] == 'HI'):
            #isolate hydrogen and info
            obsinfo={'NHI':data['logn'][sets[ii]],'eNHI':data['elogn'][sets[ii]],'hiflag':data['flag'][sets[ii]],
                     'z':data['zabs'][sets[ii]],'errz':data['ezabs'][sets[ii]],'name':sy}
        else:
            #append ion value to list
            observ.append((data['ion'][sets[ii]],data['logn'][sets[ii]],data['elogn'][sets[ii]],data['flag'][sets[ii]]))
    
    
    ########################
    ## Run the MCMC
    ########################
    
    ##Should we optimize?
    print 'Ready to run mcmc for {}'.format(sy)
    if optim in ['guess', 'guess_NHI']:
        obsinfo['met'] = float(met_guess) ##make sure it's a float
        obsinfo['dens'] = float(dens_guess) ##make sure it's a float
        obsinfo['carbalpha'] = float(carbalpha_guess) ##make sure it's a float
    else:
        optim = False
    
    logUconstraint = str(logUconstraint) ##convert to string
    nthread = int(nthread) ##make sure it's an integer
    
    ##Run
    mcmc=mcmc_ions(observ, obsinfo, grid_fil, 
                   logUconstraint=logUconstraint, UVB=UVB,
                   nwalkers=(nwalkers),
                   nsamp=(nsteps), optim=optim, threads=nthread,
                   outsave=outsave, testing=testing)
    
    print('All done!!')



################################################
################################################


def read_input_file(guesses_file, row_index_to_run):

    ##Define the columns in the input file
    ok_to_run_in = []
    sightline_in = []
    met_guess_in = []
    dens_guess_in = []
    carbalpha_guess_in = []
    carbalpha_use_in = []
    logUconstraint_use_in = []
    uvb_in = []
    notes_in = []

    ##Read the input file
    infile=open(guesses_file, 'r')
    for line in infile:
        if line and not line.startswith("#"):
            linespl = line.strip().split()
            ##
            ok_to_run_in.append(linespl[0])
            sightline_in.append(linespl[1])
            met_guess_in.append(linespl[2])
            dens_guess_in.append(linespl[3])
            carbalpha_guess_in.append(linespl[4])
            carbalpha_use_in.append(linespl[5])
            logUconstraint_use_in.append(linespl[6])
            uvb_in.append(linespl[7])
            notes_in.append(linespl[8])


    ##We only care about the row that is given as the command line argument
    ok_to_run = ok_to_run_in[row_index_to_run]
    sightline = sightline_in[row_index_to_run]
    met_guess = met_guess_in[row_index_to_run]
    dens_guess = dens_guess_in[row_index_to_run]
    carbalpha_guess = carbalpha_guess_in[row_index_to_run]
    carbalpha_use = carbalpha_use_in[row_index_to_run]
    logUconstraint_use = logUconstraint_use_in[row_index_to_run]
    uvb = uvb_in[row_index_to_run]
    notes = notes_in[row_index_to_run]



    return ok_to_run, \
        sightline, \
        met_guess, \
        dens_guess, \
        carbalpha_guess, \
        carbalpha_use, \
        logUconstraint_use, \
        uvb, \
        notes



################################################
################################################
################################################



def main(row_to_run):
    
    ##row_to_run should be indexed from 1 not 0!! We subtract 1 below!
    
    ##CBW adjust
    guesses_file = "MCMC_initial_guesses-run_me.dat"
    input_file_dir = "./input_files/"
    nthread = 12
    
    ##Re-index the row number
    row_index_to_run = int(row_to_run) - 1
    
    ##Read the input file and get the desired row
    ok_to_run, \
        sightline, \
        met_guess, \
        dens_guess, \
        carbalpha_guess, \
        carbalpha_use, \
        logUconstraint_use, \
        uvb, \
        notes = read_input_file(guesses_file, row_index_to_run)


    ##Check to make sure this row was marked to run
    if ok_to_run == 'Y':
        mcmc_infile=input_file_dir+"/mcmc."+sightline+".in"
        
        run_full_mcmc(
            infile=mcmc_infile,
            nthread=nthread,
            met_guess=met_guess,
            dens_guess=dens_guess,
            carbalpha_guess=carbalpha_guess,
            carbalpha_use=carbalpha_use,
            logUconstraint=logUconstraint_use,
            UVB=uvb)
    
    else:
        print(sightline+" was NOT marked to run!")
        
        ##Can't be interactive when running on the CRC
        # run_anyway = raw_input("Do you want to run anyway (y/n)? ")
        # if run_anyway.lower() == 'y':
        #     run_full_mcmc(
        #         infile=mcmc_infile,
        #         nthread=nthread,
        #         met_guess=met_guess,
        #         dens_guess=dens_guess,
        #         carbalpha_guess=carbalpha_guess,
        #         carbalpha_use=carbalpha_use,
        #         logUconstraint=logUconstraint_use,
        #         UVB=uvb)
    
    return



################################################
################################################


# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1:
        print("You need to specify which row (sightline) in the guesses_file to run! (Indexed from 1 not 0 to make your life easier.)")
        print("Exiting...")
        sys.exit()
    
    else:
        row_to_run = int((sys.argv)[1])
        main(row_to_run)
