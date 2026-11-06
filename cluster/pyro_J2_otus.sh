#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mail-type END
#SBATCH --mail-user noah01@physik.fu-berlin.de

#SBATCH --output=./logging/slurm_%j.out
#SBATCH --error=./logging/slurm_%j.err

# We need:
# - name
# - memory
# - time

set -eu

scontrol show job $SLURM_JOBID

module purge
module load lang
module load JuliaHPC

export JULIA_DEPOT_PATH=/scratch/hpc-prf-pm2frg/pm2frg18/.julia

julia -t ${SLURM_CPUS_PER_TASK} pyro_J2.jl "$@" # This needs to contain Jx, Jy, ...

exit 



