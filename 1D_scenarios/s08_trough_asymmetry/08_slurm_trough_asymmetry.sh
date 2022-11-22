#!/bin/bash
#SBATCH --mem-per-cpu=2G
#SBATCH --output=./log_files/asymmetry_output_%j_%a.txt
#SBATCH --error=./log_files/asymmetry_error_%j_%a.txt
#SBATCH --time=01:00:00


lvl=${1}
mod_prefix=${2}
samp=${3}
winsize=${4}
rep=${SLURM_ARRAY_TASK_ID}
mod_sufix=${5}
mapType=${6}
minRecRegionLength=${7}

ini=$(date +"%T")
echo "Current time : $ini"

module load GCC
module load R

script=08_asymmetry_in_troughs.R

chmod 777 ./${script}

# Rscript
# lvl = as.numeric(args[1])
# mod_prefix=as.character(args[2])
# samp= as.numeric(args[3])
# winsize=as.character(args[4])
# rep=as.numeric(args[5])
# mod_sufix=as.character(args[6])
# mapType=as.character(args[7])

printf "Rscript ${script} ${lvl} ${mod_prefix} ${samp} ${winsize} ${rep} ${mod_sufix} ${mapType} ${minRecRegionLength} & \n"
Rscript ${script} ${lvl} ${mod_prefix} ${samp} ${winsize} ${rep} ${mod_sufix} ${mapType} ${minRecRegionLength} &

wait

end=$(date +"%T")
echo "Current time : $end"

StartDate=$(date -u -d "$ini" +"%s")
FinalDate=$(date -u -d "$end" +"%s")
printf "Duration: "
date -u -d "0 $FinalDate sec - $StartDate sec" +"%H:%M:%S"

printf "\n end of TASK ID ${SLURM_ARRAY_TASK_ID} \n"