#!/bin/bash
####---------------------------------------------------------------------------
####   Set SLURM / Job parameters
####---------------------------------------------------------------------------
### Set the account from which HPC resources are used. For now, this is
### fixed for our project.
#SBATCH -A special00005
### Set a name identifying the job
#SBATCH -J simpleFoam 
### Do not send email notifcations, check job status with 'squeue`
#SBATCH --mail-type=NONE
### Set the number of processes. This is equivalent to the number of
### cores used usually
#SBATCH -n 4
### Set the amount of memory per core in megabyte
#SBATCH --mem-per-cpu=700
### Set a run time limit for the job. This limit must not exceed the limit
### of the chosen partition, e.g. on test30m the runtime limit is
### 30 minutes.
### The format is hh:mm:ss
#SBATCH -t 00:15:00

####---------------------------------------------------------------------------
####    The actual command to submit as job 
####---------------------------------------------------------------------------
srun simpleFoam -parallel 

