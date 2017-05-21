#!/bin/bash            ##Which shell to use to execute this file
#$ -M cwotta@nd.edu    ##Email address to send job status.
#$ -m ae               ##Email me when my job: (a)borts, (b)egins, and/or (e)nds
#$ -pe smp 12          ##Number of processors to use (factors of "12" for "long" queue, factors of "8" for "debug" queue)
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


##Actual code to run
cd ~/Lehner13-MCMC-CRC/
source ./pythonvirtualenv/bin/activate
python mcmc_met.py "${SGE_TASK_ID}"
deactivate


##Prints date+time to STDOUT
##  (to get an idea for how long this took to run)
date
