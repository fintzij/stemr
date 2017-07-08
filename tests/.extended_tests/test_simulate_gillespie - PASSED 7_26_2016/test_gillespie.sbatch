#!/bin/bash
#
# test_gillespie.sbatch
#
#SBATCH --job-name myjob      # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition largemem  # Slurm partition to use
#SBATCH --ntasks 1            # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 4-00:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=4096    # Memory limit for each tasks (in MB)
#SBATCH --mail-type=END       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=fintzij@uw.edu # Email to which notifications will be sent

# set outfile
#set outfile = "~/"test_simulate_gillespie.$sim_num.$$.Rout

# execute commands
R CMD BATCH --quiet --no-save --no-restore "--args $sim_num" test_simulate_gillespie.R test_simulate_gillespie.$sim_num.$$.Rout
