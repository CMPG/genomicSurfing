#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1
#SBATCH --output=./log_files/conCAT_4b_output_%j_%a.txt
#SBATCH --error=./log_files/conCAT_4b_error_%j_%a.txt
#SBATCH --time=01:00:00


founders=${1}
mig=${2}
growth=${3}
inds=${4}
level=${5}
janela=${6}
from=${7}
to=${8}

printf "input params: ${founders} ${mig} ${growth} ${inds} ${janela} ${level} \n"
# summary_data_lvl_10_f10_m10_g5_10inds_win_200k_194.txt
# complete_data_lvl_10_f10_m10_g5_10inds_wind_200k_194.

folder=genomic_profile_${janela}

out_1=summary_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_win_${janela}
out_2=complete_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_win_${janela}

printf "folder: ${folder}\n"
printf "out_1: ${out_1}\n"
printf "out_2: ${out_2} \n"
printf "current dir: ${PWD} \n"

if [ -f ${out_1}.txt.gz ]
	then
	rm ${out_1}.txt.gz
    rm ${out_2}.txt.gz
fi

# for rep in $(seq 1 200); do
for rep in $(seq ${from} ${to}); do
    #       summary_data_lvl_10_f10_m10_g5_10inds_win_200k_199
    infile1=summary_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_win_${janela}_${rep}
    infile2=complete_data_lvl_${level}_f${founders}_m${mig}_g${growth}_${inds}inds_wind_${janela}_${rep}
    
    printf "infile1: ${infile1} \n"
    printf "infile2: ${infile2} \n"

    if [[ -f ./${folder}/${infile1}.txt ]]
	then
        gzip -f ./${folder}/${infile1}.txt
        printf "zipping infile1 \n"
    fi

    if [[ -f ./${folder}/${infile2}.txt ]]
	then
        gzip -f ./${folder}/${infile2}.txt
        printf "zipping infile2 \n"
    fi
    
    zcat ./${folder}/${infile1}.txt.gz > ./${folder}/tmp.txt
    tail -n +2 ./${folder}/tmp.txt > ./${folder}/tmp2.txt

    cat ./${folder}/tmp2.txt >> ./${folder}/$out_1.txt

    zcat ./${folder}/${infile2}.txt.gz > ./${folder}/tmp3.txt
    tail -n +2 ./${folder}/tmp3.txt > ./${folder}/tmp4.txt

    cat ./${folder}/tmp4.txt >> ./${folder}/$out_2.txt
done

if [ -f ./${folder}/tmp2.txt ]
then
    rm ./${folder}/tmp2.txt
    rm ./${folder}/tmp.txt
fi

if [ -f ./${folder}/tmp3.txt ]
then
    rm ./${folder}/tmp3.txt
    rm ./${folder}/tmp4.txt
fi

gzip -f ./${folder}/${out_1}.txt
gzip -f ./${folder}/${out_2}.txt