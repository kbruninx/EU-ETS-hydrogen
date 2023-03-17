#!/bin/bash

#SBATCH --job-name="EU-ETS-H2"
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8GB
#SBATCH --account=lp_elect_gen_modeling
#SBATCH --array=1-4

module load 2022r2
module load julia

srun julia --threads=13 MAIN.jl --start_scen $SLURM_ARRAY_TASK_ID --stop_scen $SLURM_ARRAY_TASK_ID --start_sens 17 --stop_sens 100 > run_$SLURM_ARRAY_TASK_ID.log