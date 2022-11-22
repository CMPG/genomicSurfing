# Wed Dec 01 20:21:34 2021 ------------------------------
## ---------------------------
##
## Script name: 08_asymmetry_in_troughs.R
##
## Purpose of script: modified version for use in CLUSTER.
##  measure asymmetry in troughs 
##
## Date Created: 2021-12-01
## ---------------------------


#---------------------------------------------------------------------------------------------------
# 
# (c) Laurent Excoffier , September 2021
# 
# Estimating the null distribution of recombintion rates in chromosomal segments
#
# 09.09.2021 Changed the coordinate system to speed up the computations, but it seems to given different
#            recombination values. Originally made segments of 2kb and now switching to segments of 200 bp, 
#            which gives a much better correlation with original implementation but about 50  times faster
# 12.11.2021 Modification to examine the asymmetry of throughs that cross a low/medium rec rate
# 
#---------------------------------------------------------------------------------------------------

# -- -- -- -- -- -- -- ---
# LIBRARY & FUNCTIONS ####
# -- -- -- -- -- -- -- ---

source("06_FUNCTIONS_null_dist_troughs.R")
source("08_FUNCTIONS_asymmetry_in_troughs.R")

# -- -- -- -- ---
# IN CLUSTER ####
# -- -- -- -- ---

args=commandArgs(trailingOnly=TRUE)

lvl = as.numeric(args[1])
mod_prefix=as.character(args[2])
samp= as.numeric(args[3])
winsize=as.character(args[4])
rep=as.numeric(args[5])
mod_sufix=as.character(args[6])
mapType=as.character(args[7])
minRecRegionLength=as.numeric(args[8])


# -- -- -- -- -- 
# IN RSTUDIO ####
# -- -- -- -- --

# lvl = 10
# mod_prefix="1d_c5_l100_d100"
# samp=10
# winsize="10k"
# rep=200
# mod_sufix="f10_m10_g5"
# mapType="m1"
# minRecRegionLength=200
# 
# source("../06_tr_distribution/06_FUNCTIONS_null_dist_troughs.R")
# source("08_FUNCTIONS_asymmetry_in_troughs.R")

# -- -- -- -- -- -- -- -- -- -- -
# double checking parameter ####
# -- -- -- -- -- -- -- -- -- -- -


if (!exists("minRecRegionLength")){
  minRecRegionLength=200
}


# -- -- -- -- -- - 
#  File Names ####
# -- -- -- -- -- -

stime=Sys.time()

sprintf("lvl=%s; mod_prefix=%s; samp=%s; winsize=%s; rep=%s; mod_sufix=%s; mapType=%s; minRecRegionLength=%s", lvl, mod_prefix, samp, winsize, rep, mod_sufix, mapType, minRecRegionLength)


in_folder=sprintf("genomic_profile_%s", winsize)
inputFile=sprintf("complete_data_lvl_%s_%s_%sinds_wind_%s_%s.txt.gz", lvl, mod_sufix, samp, winsize, rep)

print(getwd())
sprintf("in file: %s/%s", in_folder, inputFile)


# -- -- -- -- -- --
# Reading data ####
# -- -- -- -- -- --

odat=read.table(paste0(in_folder, "/", inputFile), header = T)
data=data.frame(gen=odat$gen, cpop=odat$cpop, rep=odat$rep, ini=odat$ini, end=odat$end)


# ---------------------------------------------------------------------------------------------------
# Define the granularity of the analysis #Larger values wil improve speed but rec. rates 
# will be overestimated (for both simulated and observed troughs). A vlaue of 200 seems fine...
# minRecRegionLength=1

#Name of file where results are stored
out.file.name=paste("results_asymmetry_", winsize ,"/prop_troughs_in_rec_rate_category_", minRecRegionLength,"bp_MinRecRegLen_r", rep, "_", lvl,"lvl_", mod_sufix, "_", samp, "inds_", winsize, "wsize",".txt", sep="")


print(out.file.name)

# ----------------------------------------
#Get unique combinations of first three columns
cases=unique(data[,1:3])
cases[1,]
cases[2,]
trough.list=array(list(), nrow(cases))
for (i in 2:nrow(cases)) {
  low=as.numeric(row.names(cases[i-1,]))
  high=as.numeric(row.names(cases[i,]))-1
  trough.list[[i-1]]=data[low:high, 4:5]
}
#Add last case
trough.list[[nrow(cases)]]=data[(high+1):nrow(data), 4:5]

head(trough.list[[nrow(cases)]])
head(trough.list[[nrow(cases)-1]])
tail(trough.list[[nrow(cases)-1]])



#--- ok til here ...-----------------------

chromUnit = get_chr_map(mapType, minRecRegionLength)
meanObsRec=mean(chromUnit)

#Define a chromosome as a given number of repetitions of the rec structure (here 100Mb = 1000 units of 100Kb)
# chrom=rep(chromUnit, 1000)
chrom=rep(chromUnit, as.numeric(get_num_chr_units(mapType)))


sprintf("Number of ChrUnits of 100kb with a minRecRegionLength of %s: %s", minRecRegionLength, format(length(chrom), scientific = T))

# write.table(chrom, "with_map_function_chrom_object.txt", col.names = F, row.names = F)
# Define recombination rates  ----------------------
lowRec=1e-9
medRec=1e-8
highRec=1e-7

#---------------------------------------------------------------------------------------------------
# Main loop to generate new chromosome with the same number of throughs and sam total trough length 
# as observed, and computation of mean rec rates in these random troughs
#---------------------------------------------------------------------------------------------------
# t=1
for (t in 1:dim(trough.list))  {
  # for (t in 21:21)  {
  # for (t in 1:10)  {
  start_time <- Sys.time()
  
  
  # Assign current trhough for given cell and time -------------------------------------------------
  #Adopting the new coordinate system as a function of minimum rec region length
  troughs=trough.list[[t]]/minRecRegionLength
  numTroughs=nrow(troughs)
  cat(paste("Case #", t,  "/", dim(trough.list), " : ", numTroughs, " troughs\n", sep=""))
  
  
  prop.low=computeProportionOfAllThroughsInRecCategory(troughs, chrom, lowRec)
  prop.med=computeProportionOfAllThroughsInRecCategory(troughs, chrom, medRec)
  prop.high=computeProportionOfAllThroughsInRecCategory(troughs, chrom, highRec)
  
  
  out.file=file(out.file.name, open = "a")
  cat(paste(c(cases[t,], prop.low, prop.med, prop.high), sep="\t"), "\n", file=out.file)
  close(out.file)
  # cat(paste("Case ", t, "/",dim(trough.list), sep=""), "\n")
  cat(cat(paste(c(cases[t,], prop.low, prop.med, prop.high), sep="\t"), "\n"))
  
  end_time <- Sys.time()
  cat(paste("Computation time:", end_time-start_time, "\n\n"))
} 


