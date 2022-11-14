#!/bin/bash

#SBATCH --job-name="H2MSR"
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8GB
#SBATCH --account=innovation

module load 2022r2
module load openmpi
module load julia

srun julia MAIN_v2.jl > run.log