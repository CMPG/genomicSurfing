# Thu Oct 14 16:58:05 2021 ------------------------------
## ---------------------------
##
## Script name: 08_convert_time_to_diversity_loss.R
## modified from convert_time_to_diversity_loss.R
## Purpose of script: read same ol data file and instead of 
## using generations in the x-axis, use diversity loss
##
## Date Created: 2021-10-14
## ---------------------------


# -- -- -- -- -- -- -- -- --
# cluster cmd arguments ####
# -- -- -- -- -- -- -- -- --
ini_div=0.000125
args=commandArgs(trailingOnly=TRUE)

main_folder=as.character(args[1]) #"1D_fwd_RE_data"
prefix=as.character(args[2]) #"1d_c5_e5_l100_d100"
model=as.character(args[3]) #"f10_m10_g5"
sufix=as.character(args[4]) #""
file_location=as.numeric(args[5])# 1: Bern PC; 2: Ubelix; 3=IBU
from=as.numeric(args[6])
to=as.numeric(args[7])

if(ini_div != ""){
  ini_div=as.numeric(ini_div)
} else {
  ini_div=0.000125
}

sprintf("%s \n %s \n %s \n %s \n %s \n %s \n %s \n ", main_folder, prefix, model, sufix, file_location, from, to)

# -- -- -- -- -- 
# RStudio use ####
# -- -- -- -- --
# main_folder="1D_fwd_RE_data"
# prefix="1d_c5_e5_l100_d100"
# model="f20_m10_g5"
# file_location=1

# main_folder="nope"
# prefix="2d_c1_gr30_l100_d300"
# model="f10_m10_g5"
# sufix="mx100_tot_mig"
# file_location=2


# prefix="1d_c5_e5_l100_d100"
# model="f20_m10_g5"
# sufix=""
# main_folder="nope"
# file_location=2


# 2d_c1_gr30_l100_d300_f10_m10_g5_mx100_tot_mig


if(file_location == 1){
  #bern computer
  files_dir=sprintf("C:/Users/fs19b061/OneDrive - Universitaet Bern/%s/%s_%s%s/sim_files", main_folder, prefix, model, sufix)
} else if (file_location == 2){

  full_name=paste0(prefix, "_", model, "_", sufix)
  print(full_name)

  # ubelix
  if (unlist(strsplit(prefix, split = "_"))[1] == "2d"){
    nineth=unlist(strsplit(full_name, split = "_"))[9]
    print(nineth)

    if (!is.na(nineth)){
      if (nineth == "mx100"){
        if (prefix == "2d_c1_gr25_l100_d300"){
          files_dir=sprintf("/storage/homefs/fs19b061/2D_fwd_RangeExpansion/%s_%s_mx100_tot_mig/sim_files", prefix, model)
        } else {
            files_dir=sprintf("/storage/homefs/fs19b061/2D_fwd_RangeExpansion/main_models/%s_%s_mx100_tot_mig/sim_files", prefix, model)
        }
      }
      else if ((nineth == "recMap")){
        files_dir=sprintf("/storage/homefs/fs19b061/2D_recMap_fwd_RangeExpansion/%s_%s_recMap_100k/sim_files", prefix, model)
      }
    } else {
      stop(... =c("wrong path, bitch!"))
    }
  } else if (unlist(strsplit(prefix, split = "_"))[1] == "1d"){
    
    oitavo=unlist(strsplit(full_name, split = "_"))[8]
    print(oitavo)

    if (!is.na(oitavo)){
      if(oitavo == "recMap"){
        files_dir=sprintf("/storage/homefs/fs19b061/recMap_fwd_RangeExpansion/%s_%s_recMap_100k/sim_files", prefix, model)
          } else if (oitavo == "g5"){
            files_dir=sprintf("/storage/homefs/fs19b061/fwd_RangeExpansion/%s_%s/sim_files", prefix, model)
        }
      } else {
          stop(... =c("wrong path, bitch!"))}
  } # fim 1D
  
} else if (file_location == 3){
  # IBU cluster?
}

print(sprintf("pasta: %s", files_dir))


# -- -- -- -- -- -- -- -- -- -
# Libraries & Working dir ####
# -- -- -- -- -- -- -- -- -- -

library(data.table)
source("FUNCTIONS_1D_2D_PLOTS.R")


# -- -- -- -- --
# Load file ####
# -- -- -- -- --

name_col=data.frame(gen="gen",  rep="rep",  n="N", mps="mean_pi_samp", sd="sd", mda="median", qt2="qt2", qt3="qt3",  qt4="qt4",  qt1="qt1", qt5="qt5", se="se", ci="ci")

write.table(name_col, sprintf("%s/diversity_data_samp_%s_wsize_10k_10inds_all_reps.txt", files_dir, model), row.names = F, col.names = F)

write.table(name_col, sprintf("%s/diversity_data_raw_%s_wsize_10k_10inds_all_reps.txt", files_dir, model), row.names = F, col.names = F)

for(rec in seq(from, to, 1)){

  file_name=sprintf("%s/%s/genomic_profile_10k/genomic_profile_active_demes_%s_r%s_resamp_10inds.txt.gz", files_dir, rec, model, rec)
  
  if (file.exists(file_name)){
    dat=fread(file_name)
    print(rec)
    dat$ploss_raw=dat$mean_pi_raw/ini_div
    dat$ploss_samp=dat$mean_pi_samp/ini_div
    dt_samp=summarySE(data = dat, measurevar ="mean_pi_samp", groupvars = c("gen", "rep"))
    dt_raw=summarySE(data = dat, measurevar ="mean_pi_raw", groupvars = c("gen", "rep"))

    write.table(dt_samp, sprintf("%s/diversity_data_samp_%s_%s_wsize_10k_10inds_all_reps.txt",  files_dir, model, sufix), append = T, row.names = F, col.names = F)
    write.table(dt_raw, sprintf("%s/diversity_data_raw_%s_%swsize_10k_10inds_all_reps.txt", files_dir,  model, sufix), append = T, row.names = F, col.names = F)
  } else {

    print(sprintf("file: %s", file_name))
    print(sprintf("Inputfile for replicate %s DOES NOT EXISTS", rec))
    next
  }
}


