# How to run Simulations and analyse data

<!-- TOC -->

- [How to run Simulations and analyse data](#how-to-run-simulations-and-analyse-data)
- [0.0 Preparation](#00-preparation)
  - [0.1 Order of actions](#01-order-of-actions)
  - [0.2 Parameter files](#02-parameter-files)
- [1. Simulations](#1-simulations)
  - [1.1 Burn-in](#11-burn-in)
  - [1.2 1D simulations](#12-1d-simulations)
- [2. Data Analysis](#2-data-analysis)
  - [2.1 Files organization](#21-files-organization)
    - [2.1.1 Files used for each step:](#211-files-used-for-each-step)
    - [2.1.2 Files description](#212-files-description)
  - [2.2 order of things](#22-order-of-things)
  - [2.2 Allele counts: per generation \& per replicate](#22-allele-counts-per-generation--per-replicate)
  - [2.3 Genome scan: window the genome and calculate pi](#23-genome-scan-window-the-genome-and-calculate-pi)
  - [2.4 Map troughs based on a threshold](#24-map-troughs-based-on-a-threshold)
    - [2.4.1 per replicate](#241-per-replicate)
    - [2.4.2 combine files](#242-combine-files)
  - [2.5 Convert diversity data to proportion lost](#25-convert-diversity-data-to-proportion-lost)

<!-- /TOC -->


# 0.0 Preparation
## 0.1 Order of actions

1. **create new folder**: `${HOME}/${USER}/fwd_RangeExpansion`
2. **create parameter files** (for example: `params_1d_c5_e5_l100_d100_f30_m10_g5.sh`) and **put it in a subfolder**: `${HOME}/${USER}/fwd_RangeExpansion/param_files`
3. **run simulations**

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 0.2 Parameter files

`chrL`: genome length in bp
`mu`: mutation rate
`rho`: recombination rate
`core`: number of demes at the core
`maxN`: carrying capacity (max num. inds)
`mig`: migration rate in decimal (SLiM way of defining migration)
`tgrw`: time to reach caryying capacity
`r1`: 
`ert1`: growth factor at each generation
`burn`: last generation of burn in tree file
`aftgen`: generation to start expansion
`end`: last expansion generation
`Fndrs`: number of founders

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

# 1. Simulations

## 1.1 Burn-in


[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 1.2 1D simulations


```sh
# replicates ID numbers
from=1; to=200;
model_prefix="1d_c5_e5_l100_d100"
model_sufix="f20_m10_g5 10"
number_of_samples=10

cd ${HOME}/${USER}/fwd_RangeExpansion
mkdir -p log_sim
sbatch --array=${from}-${to} 01_slurm_run_model_vcf.sh ${model_sufix} ${model_prefix} ${number_of_samples}

```

This code will:
- make new folder named: `$prefix_$sufix`
- copy files to model folder:
  - `params_$prefix_$sufix.sh`
  - `01_run_model_vcf.sh`
  - `01_1D_range_exp_vcf.slim`
- create replicate number subfolders (based on array command from SLURM)
- copy slim file to replicate folder (`01_1D_range_exp_vcf.slim`)
- creates SLiM log file: `vcf_end_slim_**.output` (`**`: replicate number)
- run sims in SLiM:
  - output files: `out_p*_gen_2500*_n*.vcf.gz`. `*` in order: (1) population ID; (2) generation; (3) sample size;
  - samping frequency is defined in SLiM script, depending on the need. Normally based on variables `tg` (slim script); `tgrw` (param file), which define time to reach carrying capacity (from founders to Nmax).

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

# 2. Data Analysis

blabla
[back to top &uarr;](#how-to-run-simulations-and-analyse-data)


## 2.1 Files organization

### 2.1.1 Files used for each step:

**several steps**
- samp_edge_demes_1D_10.txt
- 10k_win_coord.csv
- FUNCTIONS_1D_2D_PLOTS.R

**step 2**
- 02_get_allele_counts_for_pi_calcul.sh
- 02_slurm_get_allele_counts_for_pi_calcul.sh


**step 3**
- 03_SLURM_Get_Pi_Profile.sh
- 03b_SLURM_conCAT_replicate_files.sh
- 03_FUNCTIONS_Get_Pi_Profile.R
- 03_Get_Pi_Profile.R

**step 4**
- 04b_slurm_ConCat_trough_files.sh
- 04a_slurm_get_troughs_per_replicate.sh
- 04a_get_troughs_per_replicate.R
- 04_FUNCTIONS_get_troughs_per_replicate.R

**step 5**
- 07_convert_time_to_diversity_loss.R
- 07_slurm_convert_time_to_diversity_loss.sh

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

### 2.1.2 Files description

- `samp_edge_demes_1D_**.txt` &rarr; `**`: number of samples; no header & space delimited;
<br>
   | generation | pop id |
   | ---------- | ------ |
   | 25002      | 6      |
   | 25007      | 7      |
   | 25012      | 8      |
   | 25017      | 9      |


- `**k_win_coord.csv` &rarr; `**`: window size in kb, e.g. 10; with header and comma separated
<br>
   | ini   | end   | wsize | mid_win | id  |
   | ----- | ----- | ----- | ------- | --- |
   | 18751 | 28751 | 10000 | 23751   | 4   |
   | 26251 | 36251 | 10000 | 31251   | 5   |
   | 33751 | 43751 | 10000 | 38751   | 6   |
   | 41251 | 51251 | 10000 | 46251   | 7   |

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 2.2 order of things

1. `cd` to model folder
2. make a sub folder: `sim_files`
3. move all replicate folders into `sim_files`
4. copy all analysis files to model folder
5. run scripts in order (see below)

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 2.2 Allele counts: per generation & per replicate

```sh

dir=${HOME}/${sim_type}/${model}
cd ${dir}
printf "cd ${dir}\n"

mkdir -p log_allele_count

printf "# $PWD \n"

sim_type="fwd_RangeExpansion";# defines the parent folder for a single spacial dist (1D or 2D) 
opt_s02=3;# defines location of sampling file
mprefix="1d_c5_e5_l100_d100"
msufix="f20_m10_g5"
model=${mprefix}_${msufix}
smp_file="10"
script_02=02_slurm_get_allele_counts_for_pi_calcul.sh

printf "sbatch --array=${list_reps} ${script_02} ${sim_type} ${opt_s02} ${model} ${smp_file} \n"

sbatch --array=${list_reps} ${script_02} ${sim_type} ${opt_s02} ${model} ${smp_file}

```
**example command line**: `sbatch --array=1-200 $script fwd_RangeExpansion 3 $model 10`

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 2.3 Genome scan: window the genome and calculate pi

```sh

dir=${HOME}/${sim_type}/${model}
cd ${dir}
printf "cd ${dir}\n"

mkdir -p log_get_pi


mprefix="1d_c5_e5_l100_d100"
msufix="f20_m10_g5"
model=${mprefix}_${msufix}

founders=20
migration=10
growth=5
samples=10
windows="10k"

active=T
old_demes=F
smp_file="10"


script_s03=03_SLURM_Get_Pi_Profile.sh

printf "${founders} ${migration} ${growth} ${samples} ${active} ${old_demes} ${windows} ${smp_file} \n"

sbatch --array=${list_reps} ${script_s03} ${founders} ${migration} ${growth} ${samples} ${active} ${old_demes} ${windows} ${smp_file}
```

**example command line**: `sbatch --array=1-200 $script 20 10 5 10 T F 10k 10`

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 2.4 Map troughs based on a threshold

### 2.4.1 per replicate

```sh

mkdir -p log_files
level=0.1
samples=10



mprefix="1d_c5_e5_l100_d100"
msufix="f20_m10_g5"
model=${mprefix}_${msufix}



script_s04a=04a_slurm_get_troughs_per_replicate.sh

printf "${model} \n"
printf "${mprefix} ${windows} \n"

mkdir -p genomic_profile_${windows}

sbatch --array=${list_reps_s03} ${script_s04a} ${level} ${mprefix} ${samples} ${windows}
```
**example command line**: `sbatch --array=1-200 $script 0.1 $mprefix 10 10k`

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

### 2.4.2 combine files

```sh

# same vars as s04a
samples=10
founders=20
migration=10
growth=5

# vars **different** from s04a
level=10
from=1
to=200

script_s04b=04b_slurm_ConCat_trough_files.sh

printf "sbatch ${script_s04b} ${founders} ${migration} ${growth} ${samples} ${level} ${windows} ${from} ${to} \n"

sbatch ${script_s04b} ${founders} ${migration} ${growth} ${samples} ${level} ${windows} ${from} ${to}

```

**example command line**: `sbatch --array=1 $script 20 10 5 10 10 10k 1 200`

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 2.5 Convert diversity data to proportion lost


```bash

target_folder="fwd_RangeExpansion";# ------|
main_folder="nope";#                       |
prefix="1d_c5_l100_d100";#                 | --> just used to define file location of data 
model="f20_m10_g5";#                       |
sufix="nope";#                             |
file_location=2;# -------------------------|
from=1
to=200
nsamp=10

slurm_script="07_slurm_convert_time_to_diversity_loss.sh"

cd /storage/homefs/fs19b061/${target_folder}/${prefix}_${model}_${sufix}

printf "${slurm_script} ${main_folder} ${prefix} ${model} ${file_location}"
sbatch --array=1 ${slurm_script} ${main_folder} ${prefix} ${model} ${sufix} ${file_location} ${from} ${to} ${nsamp}
```

**example command line**: `sbatch array=1 $script nope 1d_c5_l100_d100 f20_m10_g5 nope 2 1 200 10`

**work directory for R script** `${HOME}/${USER}/fwd_RangeExpansion/1d_c5_e5_l100_d100_f20_m10_g5/sim_files`

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)