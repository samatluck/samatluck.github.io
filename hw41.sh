#!/bin/bash
#SBATCH --qos=normal
#SBATCH --job-name=PETSc1  
#SBATCH --time=00:10:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 
#SBATCH --gres=mic:0  
#SBATCH --output=hw4.%j.out
module load intel-psxe
module load petsc/3.5.4
./benchmarkJacobi 10
