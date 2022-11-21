#!/bin/bash
#SBATCH --mem-per-cpu=5G
#SBATCH --output=./log_files/s07a_output_%j_%a.txt
#SBATCH --error=./log_files/s07a_error_%j_%a.txt
#SBATCH --time=24:00:00


module load GCC
module load R

main_folder=${1}
prefix=${2}
model=${3}
sufix=${4}
file_location=${5}
from=${6}
to=${7}

script=07_convert_time_to_diversity_loss.R

printf "\n Rscript ${script} ${main_folder} ${prefix} ${model} ${sufix} ${file_location} ${from} ${to}& \n"

Rscript ${script} ${main_folder} ${prefix} ${model} ${sufix} ${file_location} ${from} ${to} &

wait

printf "\n end of TASK ID ${SLURM_ARRAY_TASK_ID} \n"