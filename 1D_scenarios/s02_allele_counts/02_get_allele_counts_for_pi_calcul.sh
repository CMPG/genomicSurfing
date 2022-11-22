#!/bin/bash
sim_type=${1}
opt=${2}

# just deciding which path/fileName to use for the sampling file
if [ ${opt} == 2 ]; then
    model=${3}
    grid=${4}
    fpath=2d_c1_gr${grid}_l100_d100_${model}
    samp_edge_file="samp_edge_demes_gr${4}.txt"
    echo "folder path as separate values of model and grid"
else if [ ${opt} == 1 ]; then
    fpath=${3}
    samp_edge_file="samp_edge_demes_${4}.txt"
    echo "folder path as a single value"
fi
fi

echo ${fpath}

rep=${SLURM_ARRAY_TASK_ID}

# where VCF tools is located
vcftools="${HOME}/bin/vcftools"

# Model folder path
# mod_folder=${HOME}/fwd_RangeExpansion/1d_c5_e5_l100_d100_f${found}_m${mig}_g${growth}
mod_folder=${HOME}/${sim_type}/${fpath}

# Work directory for current replicate:
wdir=${mod_folder}/sim_files/${rep}/vcf_files

printf "${wdir}\n"

# active demes first...
# read line by line of sampling file (see end of loop)
while read -r line; do
    # "out_gen_25475_p956.vcf.gz"
    # grab info of which generation and which popID file to use
    gen=$(echo $line | awk -F ' ' '{print $1}')
    pop=$(echo $line | awk -F ' ' '{print $2}')
    
    # decide name of infile, depending if 1D or 2D
    if [ ${opt} == 2 ]; then
        infile=$(echo "out_gen_"${gen}"_p"${pop}".vcf.gz")
    else if [ ${opt} == 1 ]; then
        infile="out_p"${pop}"_gen_"${gen}".vcf.gz"

    fi;fi;

    # name VCFtools outfiles
    freq_file=counts2_p${pop}_gen_${gen}
    out_file=count_size_g${gen}_p${pop}_r${rep}

    printf "\n${gen}\n ${pop}\n ${infile}\n ${freq_file}\n ${out_file} \n"
    
    # outputs allele counts into ${freq_file}.frq.count
    printf "\n ${vcftools} --gzvcf ${wdir}/${infile} --counts2 --stdout > ${wdir}/${freq_file}.frq.count \n"
    ${vcftools} --gzvcf ${wdir}/"${infile}" --counts2 --stdout > ${wdir}/${freq_file}.frq.count

    # remove header?
    sed -i '1d' ${wdir}/$freq_file.frq.count

    # if outfile from previous run exists, delete it
    if [ -f ${wdir}/$out_file.txt ]; then
        rm ${wdir}/${out_file}.txt
    fi

    # keeping only columns of interest
    cat  ${wdir}/${freq_file}.frq.count | cut -f4,2,6 | paste > ${wdir}/tmp.txt
    # adding columns with pop and gen information at the end of the table
    # final table --> POS, size, allele_count, pop, gen
    sed -i "s/$/\t${pop}/" ${wdir}/tmp.txt
    sed -i "s/$/\t${gen}/" ${wdir}/tmp.txt

    # save temp file to final name (using cat)
    cat ${wdir}/tmp.txt > ${wdir}/${out_file}.txt

    # gzip outfiles; rm tmp files & freq files;
    gzip -f ${wdir}/${out_file}.txt
    rm ${wdir}/${freq_file}.frq.count
    rm ${wdir}/tmp.txt

done < ./${samp_edge_file}
# done < ./sampling_active_edge_demes.txt