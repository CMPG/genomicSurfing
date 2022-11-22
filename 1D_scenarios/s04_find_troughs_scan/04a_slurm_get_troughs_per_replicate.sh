#!/bin/bash
#SBATCH --mem-per-cpu=5G
#SBATCH --output=./log_files/s04a_output_%j_%a.txt
#SBATCH --error=./log_files/s04a_error_%j_%a.txt
#SBATCH --time=24:00:00


if [ $USER == "schlichta" ]
then
    printf "running on CMPG Matrix"
elif [ $USER == "flavia" ]
then
	printf "running on inspiron f5566"
else

	module load GCC
	module load R

fi

level=${1}
model=${2}
samps=${3}
wsize=${4}
rep=${SLURM_ARRAY_TASK_ID}

script=04a_get_troughs_per_replicate.R

printf "\n Rscript ${script} ${level} ${model} ${samps} T ${wsize} ${rep} & \n"
Rscript ${script} ${level} ${model} ${samps} T ${wsize} ${rep} &

wait

printf "\n end of TASK ID ${SLURM_ARRAY_TASK_ID} \n"