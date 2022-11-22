# Wed Dec 01 20:31:46 2021 ------------------------------
## ---------------------------
##
## Script name: 08_FUNCTIONS_asymmetry_in_troughs.R
##
## Purpose of script:
##
##
## Date Created: 2021-12-01
## ---------------------------


#---------------------------------------------------------------------------------------------------

# t=1
# troughs=trough.list[[t]]/minRecRegionLength
# numTroughs=nrow(troughs)


computeProportionOfAllThroughsInRecCategory=function(currenThroughs, curChrom, category) {
  
  
  curThrough=currenThroughs[1,]
  chromosome=curChrom
  
  computeLengthOfThroughInRecCategory=function(curThrough, chromosome, category) {
    bounds=round(curThrough)
    numPosInCategory=sum((chromosome[bounds[1]:bounds[2]]==category), na.rm=T)
    numPosInCategory
  }
  
  computeTotThroughLength=function(curThrough, chromosome, category) {
    bounds=round(curThrough)
    # print(bounds)
    troughLength=bounds[2]-bounds[1]+1
    troughLength
  }
  
  lengthTroughsinCategory=sum(apply(currenThroughs, 1, computeLengthOfThroughInRecCategory, curChrom, category))
  totTroughsLength=sum(apply(currenThroughs, 1, computeTotThroughLength))
  lengthTroughsinCategory/totTroughsLength
}

