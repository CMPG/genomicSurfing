#!/bin/bash
#SBATCH --mem-per-cpu=3G
#SBATCH --output=./log_get_pi/output_%j_%a.txt
#SBATCH --error=./log_get_pi/error_%j_%a.txt
#SBATCH --time=02:00:00

# specific to cluster
module load GCC
module load R


script=03_Get_Pi_Profile.R
founders=${1}
mig=${2}
growth=${3}
inds=${4}
act_demes=${5}
old_demes=${6}
window_size=${7}
smp_file=${8}
rep=${SLURM_ARRAY_TASK_ID}

cd sim_files/${rep}
printf "\n $1 $2 $3 $4 $5 $6 $7 $8 $rep \n"

sfolder="genomic_profile_${window_size}"
mkdir -p ${sfolder}

ckout=genomic_profile_active_demes_f${founders}_m${mig}_g${growth}_r${rep}_resamp_${inds}inds
if [ -f $ckout.txt.gz ]
	then
	rm ./${sfolder}/$ckout.txt.gz
fi

cd $sfolder

echo $PWD

# ls -l ../../

printf "\n Rscript ../../../${script} ${founders} ${mig} ${growth} ${inds} ${act_demes} ${old_demes} ${window_size} ${rep} ${smp_file} & \n" 
Rscript ../../../${script} ${founders} ${mig} ${growth} ${inds} ${act_demes} ${old_demes} ${window_size} ${rep} ${smp_file} &

wait

gzip -f ${ckout}.txt

printf "\n end of TASK ID ${SLURM_ARRAY_TASK_ID} \n"