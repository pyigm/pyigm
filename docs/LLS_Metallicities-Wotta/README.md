## Create MCMC input files

Create the file that holds the observed column densities, etc. (See examples in `input_files/`.) Make sure each name and redshift is unique so they don't overwrite each other. (I usually name them as `mcmc.SIGHTLINE_zREDSHIFT.in` so that both the sightline AND redshift info will be printed on the output plots.)

Any line beginning with "#" will be ignored (this is a good way to comment out ions, while still keeping the info).

Flags are: `0`=detection; `-1`=upper limit; `-2`=lower limit

Copy these to a directory on the CRC called `input_files/` in your working directory.



## Run MCMC

Note: This was designed to run a job array (runs of multiple different sightlines simultaneously) on Notre Dame's supercomputing cluster. A couple things might change for you.

In short, the call sequence is the following:

    queue_script.sh
        (gets options/initial guess input from MCMC_initial_guesses-run_me.dat)
        (selects appropriate line from file)
        (gets ion column density input from ./input_files/mcmc.SIGHTLINE.in)
    PyIGM/scripts/pyigm_mtlmcmc.py
    PyIGM/metallicity/mcmc.py


First, select the correct Cloudy grid in `queue_script.sh`. You shouldn't have to adjust anything, since you select the carbalpha and UVB options in the "initial guesses" file (see below). But if you need a different Cloudy grid (i.e., different parameter space, etc.), then edit `queue_script.sh` to send the correct grid to `PyIGM/scripts/pyigm_mtlmcmc.py`.

Then, edit `MCMC_initial_guesses-run_me.dat` to only include the systems you want to run. The name should correspond exactly with the name of the input file in the `input_files/` directory (without the `mcmc.` and `.in` --- following the above example, it would be of the format `SIGHTLINE_zREDSHIFT`). This "initial guesses" file provides the code with an initial guess on metallicity, an initial guess on (total hydrogen number) density, an initial guess on carbon/alpha (I normally just choose 0.0), whether you want to allow carbon/alpha to vary ("carbalpha" in the code), whether you want to use the constraint on log U ("logUconstraint" in the code), and which UVB you want to use ("HM05" or "HM12").

`PyIGM/scripts/pyigm_mtlmcmc.py` is called on a specific line (indexed from '1', so normal counting) of the `MCMC_initial_guesses-run_me.dat` file. When called, it reads in `MCMC_initial_guesses-run_me.dat`, it then reads in the `input_files/mcmc.SIGHTLINE.in`, and executes the MCMC. For example, the script `queue_script.sh` sets up a computing cluster job array that individually runs lines 1-17 (see the header) of the `MCMC_initial_guesses-run_me.dat` file. Here is a snippet from `queue_script.sh`:

    BASH$ cd ~/Lehner13-MCMC/
    BASH$ source ./pythonvirtualenv/bin/activate
    BASH$ python pyigm_mtlmcmc.py \
              --wotta \
              -guessesfile="MCMC_initial_guesses-run_me.dat" \
              -row=${SGE_TASK_ID} \
              -nthread=${NSLOTS} \
              -nwalkers=400 \
              -nsamp=400 \
    BASH$ deactivate

Any outputs from a submitted job script will be in `./MCMC/queue_script.sh.o.JOBNUMBER.SGE_TASK_ID`

If you just want to run a quick job (e.g., without allowing carbon/alpha to vary, which only takes a couple of minutes), then you can call this straight from the command line and wait for it to run. Just know that any output will be to STDOUT and will not be saved in a file (for Notre Dame's computing cluster, the default is to save to a file).

    BASH$ cd ~/Lehner13-MCMC/
    BASH$ source ./pythonvirtualenv/bin/activate
    BASH$ python pyigm_mtlmcmc.py \
              --wotta \
              -guessesfile="MCMC_initial_guesses-run_me.dat" \
              -row=1 \
              -nthread=12 \
              -nwalkers=400 \
              -nsamp=400 \
    BASH$ deactivate



## Extract MCMC output

Please see `README-understanding_MCMC_output` for how to extract MCMC output information from the `*.pkl` format files.



