#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1
#SBATCH --output=./log_files/del_4c_output_%j_%a.txt
#SBATCH --error=./log_files/del_4c_error_%j_%a.txt
#SBATCH --time=01:00:00


founders=${1}
mig=${2}
growth=${3}
inds=${4}
level=${5}
janela=${6}

printf "input params: ${founders} ${mig} ${growth} ${inds} ${janela} ${level} \n"
# summary_data_lvl_10_f10_m10_g5_10inds_win_200k_194.txt
# complete_data_lvl_10_f10_m10_g5_10inds_wind_200k_194.

folder=genomic_profile_${janela}

out_1=summary_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_win_${janela}
out_2=complete_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_win_${janela}

printf "folder: ${folder}\n"
printf "out_1: ${out_1}\n"
printf "out_2: ${out_2} \n"

if [ -f ./${folder}/${out_1}.txt.gz ]
then
    rm ./${folder}/summary_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_win_${janela}_*.txt.gz

fi

if [ -f ./${folder}/${out_2}.txt.gz ]
then
    rm ./${folder}/complete_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_wind_${janela}_*.txt.gz
fi