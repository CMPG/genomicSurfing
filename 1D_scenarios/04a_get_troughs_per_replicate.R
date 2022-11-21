# Tue Jun 22 11:36:21 2021 ------------------------------
## ---------------------------
##
## Script name: 04a_get_troughs_per_replicate.R
## previously: 04a_get_troughs.R
## previously: 16b of the same name
## Purpose of script:
##
##
## Date Created: 2021-06-22
## ---------------------------

# -- -- -- -- -- -- -- -- -- -- - 
# load libraries & functions ####
# -- -- -- -- -- -- -- -- -- -- -
library(data.table)
library(plyr)

source("./04_FUNCTIONS_get_troughs_per_replicate.R")


ini_div=0.000125

# -- -- -- -- -- 
# cmd inputs ####
# -- -- -- -- --

args=commandArgs(trailingOnly=TRUE)

level_trough = as.numeric(args[1])
model=as.character(args[2])
samps= as.numeric(args[3])
cluster=as.logical(args[4])
window_size=as.character(args[5])
rep=as.numeric(args[6])

if(!(exists("cluster"))){
  cluster=F
}

# -- -- -- -- -- -- -
# R STUDIO INPUT ####
# -- -- -- -- -- -- -
# "C:/Users/fs19b061/RangeExpansion/04_2D_model/04_get_troughs_per_replicate"

# ini_div=0.000125
# level_trough = 0.1
# model="f10_m5_g5"
# samps=10
# cluster=F
# window_size="200k"
# rep=10

# -- -- -- -- -- 
# LOAD DATA ####
# -- -- -- -- --

win_basepair = get_window_size(window_size)

if (Sys.info()[1] == "Windows"){
  df = fread(sprintf("C:/Users/fs19b061/OneDrive - Universitaet Bern/1_CMPG_matrix/010_1D_expansion/d36_data_s04_per_replicate/%s/genomic_profile_%s/genomic_profile_active_demes_%s_r%s_resamp_%sinds.txt.gz", rep, window_size, model, rep, samps))
  win_coord=fread(sprintf("./%s_win_coord.csv", window_size))
  head(win_coord)

}else if (Sys.info()[1] == "Linux"){
  df = fread(sprintf("sim_files/%s/genomic_profile_%s/genomic_profile_active_demes_%s_r%s_resamp_%sinds.txt.gz", rep, window_size, model, rep, samps))
  win_coord=fread(sprintf("%s_win_coord.csv", window_size))
  head(win_coord)
}


names(df)=c("gen", "cpop", "rep", "mean_pi_raw", "mean_pi_samp", "w_id")
head(df)
str(df)

df=merge(x=df, y=win_coord, by.x = "w_id", by.y = "id")


replicates=unique(df$rep)
used_gens=unique(df$gen)

f_dat=c()
a=Sys.time()
for (i in used_gens){
  lala=get_summary_of_troughs(df, i, level_trough, rep, win_basepair)
  if (!is.null(nrow(lala))){
    f_dat=rbind(f_dat, lala)
  } else {
    next
  }
}
  b=Sys.time()
  print(paste0("Done: replicate", rep, " - ", format(b-a, digits=4)))
  

write.table(f_dat, sprintf("genomic_profile_%s/complete_data_lvl_%s_%s_%sinds_wind_%s_%s.txt", window_size, level_trough*100, model, samps, window_size, rep), col.names = T, row.names = F)

tot_dat = ddply(f_dat, .(cpop, rep, gen), summarise, tot_win=sum(n_windows), mean_size=mean(size), tot_tr=mean(tot_tr), chr_tr=sum(size), prop_tr=sum(size)/1e8)

write.table(tot_dat, sprintf("genomic_profile_%s/summary_data_lvl_%s_%s_%sinds_win_%s_%s.txt", window_size, level_trough*100, model, samps, window_size, rep), col.names = T, row.names = F)
