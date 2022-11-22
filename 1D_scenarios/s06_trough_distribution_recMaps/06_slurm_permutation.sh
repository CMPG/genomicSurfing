#!/bin/bash
#SBATCH --mem-per-cpu=2G
#SBATCH --output=./log_files/permutation_output_%j_%a.txt
#SBATCH --error=./log_files/permutation_error_%j_%a.txt
#SBATCH --time=01:00:00


lvl=${1}
mod_prefix=${2}
samp=${3}
winsize=${4}
rep=${SLURM_ARRAY_TASK_ID}
mod_sufix=${5}
mapType=${6}
minRecRegionLength=${7}


module load GCC
module load R


script=06_null_dist_troughs.R

chmod 777 ./${script}

printf "Rscript ${script} ${lvl} ${mod_prefix} ${samp} ${winsize} ${rep} ${mod_sufix} ${mapType} ${minRecRegionLength} & \n"
Rscript ${script} ${lvl} ${mod_prefix} ${samp} ${winsize} ${rep} ${mod_sufix} ${mapType} ${minRecRegionLength} &

wait

printf "\n end of TASK ID ${SLURM_ARRAY_TASK_ID} \n"