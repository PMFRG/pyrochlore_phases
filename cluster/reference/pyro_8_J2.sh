#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu 1G
#SBATCH -J "pyro_8"
#SBATCH --mail-type END
#SBATCH --mail-user noah01@physik.fu-berlin.de

#SBATCH --output=./logging/slurm_%j.out
#SBATCH --error=./logging/slurm_%j.err
#SBATCH --array=1-4

scontrol show job $SLURM_JOBID

module purge
module load lang
module load JuliaHPC

export JULIA_DEPOT_PATH=/scratch/hpc-prf-pm2frg/pm2frg18/.julia

julia -t ${SLURM_CPUS_PER_TASK} pyro_8_J2.jl ${SLURM_ARRAY_TASK_ID}

