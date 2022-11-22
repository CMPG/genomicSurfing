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
    - [2.1.3 order of things](#213-order-of-things)
  - [2.2 Allele counts: per generation \& per replicate](#22-allele-counts-per-generation--per-replicate)
  - [2.3 Genome scan: window the genome and calculate pi](#23-genome-scan-window-the-genome-and-calculate-pi)
  - [2.4 Map troughs based on a threshold](#24-map-troughs-based-on-a-threshold)
    - [2.4.1 per replicate](#241-per-replicate)
    - [2.4.2 combine files](#242-combine-files)
  - [2.5 Convert diversity data to proportion lost](#25-convert-diversity-data-to-proportion-lost)
  - [2.6 trough distribuitions depending on recombination landscape (permutation)](#26-trough-distribuitions-depending-on-recombination-landscape-permutation)
  - [2.7 trough asymmetry](#27-trough-asymmetry)
- [3. Reference to the Figures](#3-reference-to-the-figures)
  - [3.1 Genome Scan and trough detection](#31-genome-scan-and-trough-detection)

<!-- /TOC -->


# 0.0 Preparation
## 0.1 Order of actions

1. **create new folder**: `${HOME}/${USER}/fwd_RangeExpansion`
2. **create parameter files** (for example: `params_1d_c5_e5_l100_d100_f30_m10_g5.sh`) and **put it in a subfolder**: `${HOME}/${USER}/fwd_RangeExpansion/param_files`
3. **run simulations**

[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

## 0.2 Parameter files

* `chrL`: genome length in bp
* `mu`: mutation rate
* `rho`: recombination rate
* `core`: number of demes at the core
* `maxN`: carrying capacity (max num. inds)
* `mig`: migration rate in decimal (SLiM way of defining migration)
* `tgrw`: time to reach caryying capacity
* `r1`: 
* `ert1`: growth factor at each generation
* `burn`: last generation of burn in tree file
* `aftgen`: generation to start expansion
* `end`: last expansion generation
* `Fndrs`: number of founders

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

These steps below have to happen after steps 1-4, but not in any particular order:

**step 5**
- 07_convert_time_to_diversity_loss.R
- 07_slurm_convert_time_to_diversity_loss.sh

**step 6**
- 06_slurm_permutation.sh
- 06_null_dist_troughs.R
- 06_FUNCTIONS_null_dist_troughs.R
- recMap_100k_m**.txt &rarr; recombination map file 

**step 7** (can only be run after step 6)
- 08_asymmetry_in_troughs.R
- 08_slurm_trough_asymmetry.sh
- 08_FUNCTIONS_asymmetry_in_troughs.R





[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

### 2.1.2 Files description

- `samp_edge_demes_1D_**.txt` &rarr; `**`: number of samples; no header & space delimited;



   | generation | pop id |
   | ---------- | ------ |
   | 25002      | 6      |
   | 25007      | 7      |
   | 25012      | 8      |
   | 25017      | 9      |



- `**k_win_coord.csv` &rarr; `**`: window size in kb, e.g. 10; with header and comma separated


   | ini   | end   | wsize | mid_win | id  |
   | ----- | ----- | ----- | ------- | --- |
   | 18751 | 28751 | 10000 | 23751   | 4   |
   | 26251 | 36251 | 10000 | 31251   | 5   |
   | 33751 | 43751 | 10000 | 38751   | 6   |
   | 41251 | 51251 | 10000 | 46251   | 7   |



- `recMap_100k_m**.txt` &rarr; `**`: identifier of type of map


   | iniPOS  | endPOS  | rec.rate    | chunk id |
   | ------- | ------- | ----------- | -------- |
   | 0       | 299999  | 0.000000001 | 1        |
   | 300000  | 899999  | 0.00000001  | 2        |
   | 900000  | 999999  | 0.0000001   | 3        |
   | 1000000 | 1299999 | 0.000000001 | 4        |
   | 1300000 | 1899999 | 0.00000001  | 5        |
   | 1900000 | 1999999 | 0.0000001   | 6        |




[back to top &uarr;](#how-to-run-simulations-and-analyse-data)

### 2.1.3 order of things

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

## 2.6 trough distribuitions depending on recombination landscape (permutation)


```bash
lvl=10
mod_prefix="1d_c5_l100_d100"
samp=10
winsize="10k"
mod_sufix="f20_m10_g5"
main_dir="${HOME}/${USER}/recMap_fwd_RangeExpansion"

mod_tail_LIST=("recMap_100k_m2" "recMap_100k_m3")

for mod_tail in "${mod_tail_LIST[@]}"; do
  echo ${mod_tail}

  cd ${main_dir}/${mod_prefix}_${mod_sufix}_${mod_tail}
  
  echo ${PWD}

  mkdir -p results_permutation_10k
  
  sbatch --array=1-200%50 06_slurm_permutation.sh ${lvl} ${mod_prefix} ${samp} ${winsize} ${mod_sufix} ${mapType}
done

```


Resulting file:

`proba_rec_rate_in_troughs_lvl`**troughThresholdAsPercentage**`_`**completeModelName**`_`**sampleSize**`inds_win_`**windowsize**`k_r`**rep**`.txt`

| gen   | cpop | rep | obs.mean.rec         | probability       |
| ----- | ---- | --- | -------------------- | ----------------- |
| 25003 | 4    | 132 | 1.65147443650357e-08 | 0.636666666666667 |
| 25008 | 9    | 132 | 1.63510992993476e-08 | 0.896666666666667 |
| 25013 | 14   | 132 | 1.59076823935558e-08 | 1                 |
| 25018 | 19   | 132 | 1.61938738121737e-08 | 0.846666666666667 |
| 25023 | 24   | 132 | 1.59505868076624e-08 | 0.963333333333333 |
| 25028 | 29   | 132 | 1.62395872348557e-08 | 0.643333333333333 |



` permut_mean_rec_rate_in_troughs_by_case_lvl`**troughThresholdAsPercentage**`_`**completeModelName**`_`**sampleSize**`inds_win_`**windowSize**`k_r`**repNumber**`.txt`

| gen   | cpop | rep | perm.mean.rec        | obs.mean.rec         | is.perm.rate.higher | diff                 |
| ----- | ---- | --- | -------------------- | -------------------- | ------------------- | -------------------- |
| 25003 | 4    | 132 | 1.67645079714128e-08 | 1.65147443650357e-08 | TRUE                | 2.49763606377155e-10 |
| 25003 | 4    | 132 | 1.65731280923584e-08 | 1.65147443650357e-08 | TRUE                | 5.83837273227164e-11 |
| 25003 | 4    | 132 | 1.66954370533259e-08 | 1.65147443650357e-08 | TRUE                | 1.80692688290289e-10 |
| 25003 | 4    | 132 | 1.68642550852116e-08 | 1.65147443650357e-08 | TRUE                | 3.49510720175962e-10 |



## 2.7 trough asymmetry

```bash
lvl=10
mod_prefix="1d_c5_l100_d100"
samp=10
winsize="10k"
main_dir="${HOME}/${USER}/recMap_fwd_RangeExpansion"
mod_tail="recMap_100k"
minRecRegionLength=200
# mapType="m4"
# mapType="m6"
mapType="m7"
mod_sufix_LIST=("f20_m10_g5")

for mod_sufix in "${mod_sufix_LIST[@]}"; do
  echo ${mod_sufix}

  cd ${main_dir}/${mod_prefix}_${mod_sufix}_${mod_tail}_${mapType}
  
  echo ${PWD}

  mkdir -p results_asymmetry_${winsize}
 
  sbatch --array=1-200 08_slurm_trough_asymmetry.sh ${lvl} ${mod_prefix} ${samp} ${winsize} ${mod_sufix} ${mapType} ${minRecRegionLength}
done
```

# 3. Reference to the Figures

## 3.1 Genome Scan and trough detection



<!-- Figure 1 [[./Figures/MBE_SchlichtaF_Figure1_RGB.tiff]]

:::image source="./Figures/MBE_SchlichtaF_Figure1_RGB.tiff":::


![Figure 1][test]


[test]: Figures/MBE_SchlichtaF_Figure1_RGB.png "testesttest"


![Figure 1](genomicSurfing/Figures/MBE_SchlichtaF_Figure1_RGB_Figure1.tif)


![MBE_SchlichtaF_Figure1_RGB_Figure1](/assets/MBE_SchlichtaF_Figure1_RGB_Figure1.tif) -->

<!-- 
![alt text][logo]

[logo]: https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 2" -->