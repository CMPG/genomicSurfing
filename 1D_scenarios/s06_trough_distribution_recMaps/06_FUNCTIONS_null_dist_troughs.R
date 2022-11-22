## ---------------------------
## ORIGINAL PATH: C:\Users\fs19b061\RangeExpansion\05_recomb_map\06_tr_distribution
## Script name: 06_FUNCTIONS_null_dist_troughs.R
##
## Purpose of script: Collecting functions from original script:
## "heteroRecRateNull-fast.r" (LE version) -- > 
##                                          "06_null_dist_troughs.R" (Cluster version)
## Date Created: 2021-09-15
## ---------------------------


# -- -- -- -- -- 
# (TITLE HERE) ####
# -- -- -- -- --

get_num_chr_units = function(mapType){
   if ((mapType == "m1")|(mapType == "m2")|(mapType == "m3")|(mapType == "m4")){
      num_chr_units = 1000
   } else if ((mapType == "m6") | (mapType == "m7")){
      num_chr_units = 100
   }
   return(num_chr_units)
}



# If recMap is m1 (first model)

get_chr_map = function(mapType, minRecRegionLength){

  if (mapType == "" | mapType == "m1" | is.na(mapType)){
    # Define recombination rates  ----------------------
    lowRec=1e-9
    medRec=1e-8
    highRec=1e-7
    
    # Define chromosome structure   ----------------------
    lRregion=30e3
    mRregion=10e3
    hRregion=2e3
    
    lRregion=lRregion/minRecRegionLength
    mRregion=mRregion/minRecRegionLength
    hRregion=hRregion/minRecRegionLength
    
    chromUnit=c(rep(lowRec,lRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion)
    )
  } else if (mapType == "m2"){
    lowRec=1e-9
    medRec=1e-8
    
    lRregion=30e3
    mRregion=70e3
    
    
    lRregion=lRregion/minRecRegionLength
    mRregion=mRregion/minRecRegionLength
    
    chromUnit=c(rep(lowRec,lRregion),
                rep(medRec, mRregion)
                )
    
  } else if (mapType == "m3"){
    
    lowRec=1e-9
    medRec=1e-8
    
    lRregion=50e3
    mRregion=50e3
    
    lRregion=lRregion/minRecRegionLength
    mRregion=mRregion/minRecRegionLength
    
    chromUnit=c(rep(lowRec,lRregion),
                rep(medRec, mRregion)
                )
  } else if (mapType == "m4"){
    # Define recombination rates  ----------------------
    lowRec=1e-9
    medRec=1e-8
    highRec=1e-7
    
    # Define chromosome structure   ----------------------
    lRregion=30e3
    mRregion=60e3
    hRregion=10e3
    
    lRregion=lRregion/minRecRegionLength
    mRregion=mRregion/minRecRegionLength
    hRregion=hRregion/minRecRegionLength
    
    chromUnit=c(rep(lowRec,lRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion)
                )

  } else if (mapType == "m6"){
         # Define recombination rates  ----------------------
    lowRec=1e-9
    medRec=1e-8
    highRec=1e-7
    
    # Define chromosome structure   ----------------------
    lRregion=300e3
    mRregion=100e3
    hRregion=20e3
    
    lRregion=lRregion/minRecRegionLength
    mRregion=mRregion/minRecRegionLength
    hRregion=hRregion/minRecRegionLength
    
    chromUnit=c(rep(lowRec,lRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion),
                rep(medRec, mRregion)
    )

  } else if (mapType == "m7"){
         # Define recombination rates  ----------------------
    lowRec=1e-9
    medRec=1e-8
    highRec=1e-7
    
    # Define chromosome structure   ----------------------
    lRregion=300e3
    mRregion=600e3
    hRregion=100e3
    
    lRregion=lRregion/minRecRegionLength
    mRregion=mRregion/minRecRegionLength
    hRregion=hRregion/minRecRegionLength
    
    chromUnit=c(rep(lowRec,lRregion),
                rep(medRec, mRregion),
                rep(highRec, hRregion)
               )
  }
  
  return(chromUnit)
}

# -- -- -- -- -- 
# FUNCTIONS ####
# -- -- -- -- --
# 
# trough.pos=troughs[1:3,]
# totChromLength=length(chrom)
Comp.obs.intervals=function(trough.pos, totChromLength) {
   # print(head(trough.pos))
   num.troughs=nrow(trough.pos)
   inter=array(0, dim=c(num.troughs+1, 2))
   int.beg=1/minRecRegionLength
   int.end=trough.pos[1,1] #- 1/minRecRegionLength
   if (int.end>int.beg) {
      inter[1,1]= int.beg
      inter[1,2]= int.end
   }
   for (i in 2:num.troughs) {
      inter[i,1]=trough.pos[i-1, 2] #+ 1/minRecRegionLength
      inter[i,2]=trough.pos[i, 1] #- 1/minRecRegionLength
   }
      
   inter[num.troughs+1,1]=trough.pos[num.troughs, 2] #+ 1/minRecRegionLength
   inter[num.troughs+1,2]=totChromLength
   inter[,2]-inter[,1]
}


findTroughEndPoints=function(troughLengths, intervals) {
   numTroughs=length(troughLengths)
   begs=rep(0,numTroughs)
   ends=begs
   begs[1]=intervals[1] #+ 1/minRecRegionLength
   ends[1]=begs[1]+troughLengths[1]
   for (i in 2:numTroughs) {
      begs[i]=ends[i-1]+intervals[i] #+ 1/minRecRegionLength
      ends[i]=begs[i]+troughLengths[i]
   }
   list(tbegs=begs, tends=ends)
}

# Function computing the mean rec rate in randomly occuring segments of same length as those observed
# troughLengths = troughLengths[troughOrder]
# intervals=intervals[,1]

computeMeanRecRate=function(troughLengths, intervals) {
   trough.pos=findTroughEndPoints(troughLengths, intervals)
   numTroughs=length(troughLengths)
   
   sumRec=0
   for (i in 1:numTroughs) {
      sumRec=sumRec+sum(chrom[trough.pos$tbegs[i]:trough.pos$tends[i]])   
   }
   meanRec=sumRec/sum(troughLengths)
   meanRec
}

numTroughsWithGivenRecRegions=function(troughLengths, intervals) {
   trough.pos=findTroughEndPoints(troughLengths, intervals)
   numTroughs=length(troughLengths)
   num.ltr=0
   num.mtr=0
   num.htr=0
   for (i in 1:numTroughs) {
      if (length(which(chrom[trough.pos$tbegs[i]:trough.pos$tends[i]]==lowRec))) num.ltr=num.ltr+1
      if (length(which(chrom[trough.pos$tbegs[i]:trough.pos$tends[i]]==medRec))) num.mtr=num.mtr+1
      if (length(which(chrom[trough.pos$tbegs[i]:trough.pos$tends[i]]==highRec))) num.htr=num.htr+1  
   }
   list(num.ltr, num.mtr, num.htr)
}

