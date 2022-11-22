## ---------------------------
##
## Script name: 06_null_dist_troughs.R
##
## Purpose of script: see below; edited version (to run in cluster) by me
## Original code from LE.
##
## Date Created: 2021-09-15
## ---------------------------

# -- -- -- -- -- -- -- ---
# LIBRARY & FUNCTIONS ####
# -- -- -- -- -- -- -- ---

source("06_FUNCTIONS_null_dist_troughs.R")



#---------------------------------------------------------------------------------------------------
# 
# (c) Laurent Excoffier , September 2021
# 
# Estimating the null distribution of recombintion rates in chromosomal segments
#
# 09.09.2021 Changed the coordinate system to speed up the computations, but it seems to given different
#            recombination values. Originally made segments of 2kb and now switching to segments of 200 bp, 
#            which gives a much better correlation with original implementation but about 50  times faster
# 
#---------------------------------------------------------------------------------------------------


# Reading input file  ------------------------------------------------------------------------------

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

if (!exists("minRecRegionLength")){
  minRecRegionLength=200
}


# -- -- -- -- -- 
# IN RSTUDIO ####
# -- -- -- -- --

lvl = 10
mod_prefix="2d_c1_gr15_l100_d300"
samp=10
winsize="10k"
rep=132
# rep=NA
mod_sufix="f10_m10_g5"
mapType="m1"
minRecRegionLength=200
# source("../06_FUNCTIONS_null_dist_troughs.R")

# > lvl
#  10
# > mod_prefix
# "2d_c1_gr15_l100_d300"
# > samp
#  10
# > winsize
#  "10k"
# > rep
# 132
# > mod_sufix
# "f10_m10_g5"
# > mapType
# "m1"
# > minRecRegionLength
# 200


# -- -- -- -- -- - 
#  File Names ####
# -- -- -- -- -- -

stime=Sys.time()

# sprintf("lvl=%s; mod_prefix=%s; samp=%s; winsize=%s; rep=%s; mod_sufix=%s;", lvl, mod_prefix, samp, winsize, rep, mod_sufix)
sprintf("lvl=%s; mod_prefix=%s; samp=%s; winsize=%s; rep=%s; mod_sufix=%s; mapType=%s; minRecRegionLength=%s", lvl, mod_prefix, samp, winsize, rep, mod_sufix, mapType, minRecRegionLength)

in_folder=sprintf("genomic_profile_%s", winsize)

if (!is.na(rep)){
  inputFile=sprintf("complete_data_lvl_%s_%s_%sinds_wind_%s_%s.txt.gz", lvl, mod_sufix, samp, winsize, rep)
} else {
  inputFile=sprintf("complete_data_lvl_%s_%s_%sinds_win_%s.txt.gz", lvl, mod_sufix, samp, winsize, rep)
}


print(getwd())
# setwd("..")
sprintf("in file: %s/%s", in_folder, inputFile)


#Name of file where results are stored
out.file.name=sprintf("results_permutation_%s/proba_rec_rate_in_troughs_lvl%s_%s_%s_%sinds_win_%s_r%s.txt", winsize, lvl, mod_prefix, mod_sufix, samp, winsize, rep)

dist.file.name=sprintf("results_permutation_%s/permut_mean_rec_rate_in_troughs_by_case_lvl%s_%s_%s_%sinds_win_%s_r%s.txt", winsize, lvl, mod_prefix, mod_sufix, samp, winsize, rep)

sprintf("outfile: %s", out.file.name)

sprintf("distribution file: %s", dist.file.name)


# -- -- -- -- -- --
# Reading data ####
# -- -- -- -- -- --

odat=read.table(paste0(in_folder, "/", inputFile), header = T)
# odat=read.table("./genomic_profile_10k/complete_data_lvl_10_f10_m10_g5_10inds_win_10k.txt.gz", header = T)
                # "./genomic_profile_10k/complete_data_lvl_10_f10_m10_g5_10inds_wind_10k.txt.gz"
# format data to correspond to only the info used...
# GOAL:
# "gen" "cpop" "rep" "ini" "end"
# 25002 6 1 206251 216251
# 25002 6 1 243751 253751

# Original format:
# |"t_id"|"cpop"|"rep"|"gen"|"n_windows"|"mp_raw"|"mp_samp"|"ini"|"end"|"mid_tr"|"size"|"tot_tr"|

data=data.frame(gen=odat$gen, cpop=odat$cpop, rep=odat$rep, ini=odat$ini, end=odat$end)


# -- -- -- -- -- -- 
#  Study cases ####
# -- -- -- -- -- --

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



# mapType="m7"
# minRecRegionLength=100000


chromUnit = get_chr_map(mapType, minRecRegionLength)
meanObsRec=mean(chromUnit)

#Define a chromosome as a given number of repetitions of the rec structure (here 100Mb = 1000 units of 100Kb)
chrom=rep(chromUnit, as.numeric(get_num_chr_units(mapType)))


#Define number of collections of troughs randomly distributed over the genome  ---------------------
numRepeats=300

#---------------------------------------------------------------------------------------------------
# Main loop to generate new chromosome with the same number of throughs and same total trough length 
# as observed, and computation of mean rec rates in these random troughs
#---------------------------------------------------------------------------------------------------


out.file=file(out.file.name)
cat(paste(c(colnames(cases[1,]), "obs.mean.rec", "p-val")), "\n", file=out.file)
close(out.file)

dist.file=file(dist.file.name)
cat(paste(c(colnames(cases[1,]), "perm.mean.rec", "obs.mean.rec", "is.perm.rate.higher", "diff")), "\n", file=dist.file)
close(dist.file)

for (t in 1:dim(trough.list))  {
    # for debugging
   # for (t in 1:2)  {
   start_time <- Sys.time()
   
   
   # Assign current trhough for given cell and time -------------------------------------------------
   #Adopting the new coordinate system as a function of minimum rec region length
   troughs=trough.list[[t]]/minRecRegionLength
   numTroughs=nrow(troughs)
   
   # Compute obs through lengths -------------------------------------------------------------------
   compLength=function(trough) trough[2]-trough[1]
   troughLengths=apply(troughs, 1, compLength)
   totTroughLength=sum(troughLengths)
   freeChromLength=length(chrom)-totTroughLength
   
   cat(paste("Case",t, "with", numTroughs, "throughs of total length",totTroughLength/1000*minRecRegionLength,"kb\n"))
   
   numTroughFreeSegments=numTroughs+1
   # numTroughFreeSegments
   #Compute obs mean rec rate ----------------------------------------------------------------------
   obs.intervals=Comp.obs.intervals(troughs, length(chrom))
   mean(obs.intervals)
   sd(obs.intervals)
   obs.mean.rec=computeMeanRecRate(troughLengths, obs.intervals)
   
   # Assign random positions to troughs according to multinomial probability -----------------------
   intervals=rmultinom(numRepeats, size=freeChromLength*minRecRegionLength,
                   prob=rep(1.0/numTroughFreeSegments,numTroughFreeSegments))/minRecRegionLength
   mean(intervals)
   # sum(intervals)
   
   
   # tail(obs.intervals)
   # tail(troughLengths)
   
   #Compute mean rec of randomly located troughs...
   gt.rec=0
   for (i in 1:numRepeats) {
      
      #Permute trough order
      troughOrder=sample(1:numTroughs, numTroughs)
      
      curMeanRec=computeMeanRecRate(troughLengths[troughOrder], intervals[,i])
      # print(curMeanRec)

      dist.file=file(dist.file.name, open = "a")
      cat(paste(c(cases[t,], curMeanRec, obs.mean.rec, curMeanRec>obs.mean.rec, curMeanRec-obs.mean.rec), sep=""), "\n", file=dist.file)
      close(dist.file)


      if (curMeanRec>obs.mean.rec) gt.rec=gt.rec+1
      if (!i%%10) cat("\r",paste("step ", i, "/",numRepeats, sep=""))
   }
   gt.rec
   
   #Compute proba to have a larger mean rec rate in random troughs than the observed one -----------
   pval=gt.rec/numRepeats
   
   out.file=file(out.file.name, open = "a")
   cat(paste(c(cases[t,], obs.mean.rec, pval), sep=""), "\n", file=out.file)
   close(out.file)
   cat(paste("\nCase ", t, "/",dim(trough.list), sep=""), "\n")
   cat(paste("The probability that random troughs of same lengths randomly distributed over the genome have a larger rec. rate than that observed is", pval, "\n"))
   

   end_time <- Sys.time()
   cat(paste("Iteration time:", end_time-start_time, "\n"))
} 



etime=Sys.time()

print(sprintf("total time for this run: %s", format(etime-stime, units="secs")))


#Reading and ploting results  ----------------------------------------------------------------------

# res.file=out.file.name
# p.values=read.table(res.file, header=T)
# head(p.values)
# which(p.values[,5]<0.05)
# which(p.values[,5]>0.95)
# hist(p.values[,5], xlab="Proba. rand mean rec rate > obs mean rec rate", main = "Probability of trough mean rec rate lower than expected")


# # Comparison with previous analyses for testing only -------------------
# ori.res=read.table("proba_gt_rand_rec_rate_in_troughs-1bp-win.txt", header=T)
# hist(ori.res[,5])
# plot(p.values[,5], ori.res[,5])
# 
# old.res.2=read.table("proba_gt_rand_rec_rate_in_troughs-100bp-win-rep1.txt", header=T)
# hist(old.res.2[,5])
# plot(p.values[,5], old.res.2[,5])
# 
# old.res.2=read.table("proba_gt_rand_rec_rate_in_troughs-200bp-win-rep1.txt", header=T)
# hist(old.res.2[,5])
# plot(p.values[,5], old.res.2[,5])
