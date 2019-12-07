#!/bin/bash
#SBATCH --job-name=ex19omp    ### Job Name
#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)
#SBATCH --time=0-00:10:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=10    ### Number of cores per task

# If Intel C++ is used, need this
module load intel-psxe/2015-update1

# Setting the number of threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./ex19omp 50
