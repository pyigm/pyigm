#!/bin/bash
#$ -M my@email.com     ##Email address to send job status.
#$ -m ae               ##Email me when my job: (a)borts, (b)egins, and/or (e)nds
#$ -pe smp 12          ##Number of processors to use (factors of "12" for "long" queue, factors of "8" for "debug" queue). The number will be saved/usable as "${NSLOTS}".
#$ -q long             ##Which queue/machine to use (can be "long" or "debug")
#$ -r n                ##If the job is killed, should it try to auto-restart? (y/n)
#$ -t 1-17             ##Define an array job. CANNOT START AT "0". The numbers will be saved/usable as "${SGE_TASK_ID}".


##Note: ${SGE_TASK_ID} comes from whatever the
##  "-t" (from above) is at in the job array.
##It starts at 1 (it CANNOT start at 0)

##Prints script name to STDOUT
basename "$0"
##Prints date+time to STDOUT
##  (to get an idea for how long this took to run)
date
##Flush open file buffers every 300 sec (5 min) so
##  you can see any STDOUT BEFORE the job is completed
fsync -d 300 $SGE_STDOUT_PATH &


########################
##Actual code to run
########################

##Change to the desired directory
cd ~/MCMC_Project/

##Define the python virtual environment location
##  (can be relative to the above path)
pyvenv_loc="./pythonvirtualenv"

##Activate the python virtual environment
source "${pyvenv_loc}"/bin/activate


##Find the "pyigm_mtlmcmc.py" within the virtual environment.
##  You could type this by hand if you'd rather.
# pyigm_mtlmcmc_loc=$(find "${pyvenv_loc}" -name "pyigm_mtlmcmc.py")

##There should only be ONE found... hopefully. But just
##  in case not, we can do this, which also allows for
##  spaces in file path names (ick!). We take the [0] element
##  of the array below.
##  https://blog.famzah.net/2016/10/20/bash-process-null-terminated-results-piped-from-external-commands/
##  https://stackoverflow.com/questions/1116992/capturing-output-of-find-print0-into-a-bash-array
pyigm_mtlmcmc_loc=()
while IFS='' read -r -u"$FD" -d $'\0' file; do
    pyigm_mtlmcmc_loc+=("$file")
done {FD}< <(find -L "${pyvenv_loc}" -name "pyigm_mtlmcmc.py" -print0)


##Use a job array to loop over the row number
##  (by getting the job array instance number)
##  and automatically get nthread
##  (number of processors to use; it will get from the header).
##
##Note: if you want to use fewer than 400 walkers
##  (nwalkers) or 400 steps (nsamp), you need to
##  include the "--testing" option. See all options below.
python "${pyigm_mtlmcmc_loc[0]}" \
    --wotta \
    -grid="/path/to/grid_cgm" \
    -guessesfile="MCMC_initial_guesses-run_me.dat" \
    -row=${SGE_TASK_ID} \
    -nthread=${NSLOTS} \
    -nwalkers=400 \
    -nsamp=400

deactivate


########################
########################


##Prints date+time to STDOUT
##  (to get an idea for how long this took to run)
date



########################
##All of pyigm_mtlmcmc.py options:
########################
##
## -sightline=__       type=str, Name of the System to analyze (auto-detected if using --wotta)
## -fileinput=__       type=str, The observed column density information (auto-detected if using --wotta)
## -outsave=__         type=str, The directory in which to save output
## -grid=__            type=str, Full path+basename to Cloudy grid.
##                               E.g., "/afs/crc.nd.edu/group/CGMND/Cloudy_grids/grid_cgm_extensive" (this is the default).
##                               Do NOT include UVB, carbalpha, or .pkl extension; you specify these with
##                               the other options, on a per-absorber basis.
## -logUconstraint=__  type=str, Should we use logU constraint on density? (auto-detected if using --wotta) Can be: "True", "False",
##                               or comma-separated values for logUmean and logUsigma (e.g., "-3.1,0.2"), which assumes "True".
##                               NOTE: This final (CSV) format overrides -logUmean and -logUsigma options!
## -logUmean=__        type=str, If we use logUconstraint, what is the mean for the Gaussian?
##                               NOTE: This can be overridden by -logUconstraint option!
## -logUsigma=__       type=str, If we use logUconstraint, what is the sigma for the Gaussian?
##                               NOTE: This can be overridden by -logUconstraint option!
## -UVB=__             type=str, The UVB to use when we are using the logU constraint on density (auto-detected if using --wotta)
## -nthread=__         type=int, Number of threads
## -nwalkers=__        type=int, Number of walkers
## -nsamp=__           type=int, Number of samples/steps
## -optim=__           type=str, Optimization method
## -dens=__            type=float, Guess at density (optim=guess) (auto-detected if using --wotta)
## -met=__             type=float, Guess at metallicity (optim=guess) (auto-detected if using --wotta)
## -carbalpha=__       type=float, Guess at carbalpha (optim=guess) (auto-detected if using --wotta)
## --testing           If used, will over-ride minimum nwalkers and minimum nsamp
## --wotta             If used, reads in files using Wotta's file format (there is a guesses file,
##                     and each sightline has separate input file). If specified, the guesses file contains: dummy column at the front (not used here);
##                     the sightline name; metallicity initial guess; density initial guess; carbon/alpha ratio (carbalpha) initial guess
##                     (required, even if you don't want to use carbalpha); whether or not to allow carbalpha to vary;
##                     whether to use the Wotta+16 logUconstraint as a prior on the density; the UVB to use; ions to comment out (not used here); 
##                     and any additional notes (not used here).
##                     Then, all that needs to be specified here is: -guessesfile=__; -row=__ (in the guessesfile, usually
##                     automatically done by the supercomputer submission script); -nthread=__ (also done by the submission script);
##                     -nwalkers=__; and -nsamp=__.
## -guessesfile=__     type=str, The --wotta options input file holding the initial guesses
## -row=__             type=str, Row (sightline) in the guesses file to run, one-indexed (not zero-indexed)
##


