#!/bin/bash
#SBATCH --output=./log_sim/output_%j_%a.txt
#SBATCH --error=./log_sim/error_%j_%a.txt
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --tmp=5G


prefix=${1}
model=${prefix}_${2}
nsamps=${3}



fwd_path="${HOME}/$USER/fwd_RangeExpansion"

# cluster specific: needs to have gcc library loaded
module load GCC


models=(${model})


bash_file="01_run_model_vcf.sh"
slim_file="01_1D_range_exp_vcf.slim"

for m in "${models[@]}"; do
    mkdir -p ./$m
    
    files=(${bash_file} ${slim_file} "params_${m}.sh")
    
    for f in "${files[@]}"; do
        cp $f $fwd_path/$m/$f
    done
done


chmod 777 $fwd_path/*.sh
chmod 777 $fwd_path/*/*.sh

srun --ntasks=1 $fwd_path/$model/${bash_file} $model ${SLURM_ARRAY_TASK_ID} ${nsamps} &

wait


printf "\n end of sim TASK ID ${SLURM_ARRAY_TASK_ID} \n"