#!/bin/bash
#SBATCH --output=./log_allele_count/output_%j_%a.txt
#SBATCH --error=./log_allele_count/error_%j_%a.txt
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=5G

sim_type=${1}
opt=${2}
model=${3}
grid_or_1D=${4}

rep=${SLURM_ARRAY_TASK_ID}

script=02_get_allele_counts_for_pi_calcul.sh

chmod 777 ${script}

# srun --ntasks=1 ${script} ${fnd} ${mig} ${grw} &

printf "\n srun --ntasks=1 ${script} ${sim_type} ${opt} ${model} ${grid_or_1D} & \n"
srun --ntasks=1 ${script} ${sim_type} ${opt} ${model} ${grid_or_1D} &

wait

printf "\n end of sim TASK ID ${SLURM_ARRAY_TASK_ID} \n"