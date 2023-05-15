#!/bin/bash

#SBATCH --job-name="EU-ETS-H2"
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8GB
#SBATCH --account=research-tpm-ess
#SBATCH --array=27-50

module load 2022r2
module load julia

srun julia --threads=17 MAIN.jl --start_scen $SLURM_ARRAY_TASK_ID --stop_scen $SLURM_ARRAY_TASK_ID --start_sens 1 --stop_sens 100 > run_$SLURM_ARRAY_TASK_ID.log