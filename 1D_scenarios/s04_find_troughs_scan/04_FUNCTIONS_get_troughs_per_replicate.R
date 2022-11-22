# Tue Jun 22 11:45:13 2021 ------------------------------
## ---------------------------
##
## Script name: 04_FUNCTIONS_get_troughs_per_replicate.R
## previously: 04_FUNCTIONS_get_troughs.R
## (modified b15 of same name, cloud folder: 
## "010_1D_expansion\d30_sampling_effect)"
##
## Purpose of script:
## Date Created: 2021-06-22
## ---------------------------

# -- -- -- -- -- 
# FUNCTIONS ####
# -- -- -- -- --

get_window_size = function(x){
  if (x=="large"){
    return(2e6)
  } else if (x == "med" | x == "200k"){
    return(2e5)
  } else if (x == "small"){
    return(2e4)
  } else if (x == "100k"){
    return(1e5)
  } else if (x == "50k"){
    return(5e4)
  } else if (x == "10k"){
    return(1e4)
  }
}

  
# dummy vals for debug:
# cur_gen=used_gens[50]
# lvl=0.1
# cur_rep=replicates[1]
# win_basepair=2e6

# goes over genome scan windows, sequentially.
# those bellow threshold, will be classified as troughs
# windows that are next to each other and are bellow treshold, will be considered as a single trough
# windows coordinate: ini | end  | wsize | mid_win     | id
# genomic profile:    gen | cpop | rep   | mean_pi_raw | mean_pi_samp | w_id
# df is the merge of winCoord + genomic profile, by (id, w_id)
# df = gen, cpop, rep, mean_pi_raw, mean_pi_samp, + w_id + ini, end, wsize, mid_win

calc_troughs = function(df, cur_gen, lvl, cur_rep, win_basepair){
  # select windows from a single generation, single replicate and that are below the threshold
  dat=df[((df$gen == cur_gen) & (df$rep == cur_rep) & (df$mean_pi_samp <= ini_div*lvl)),]
  
  # if there are windows below the threshold...
  if (nrow(dat)>0){
    i=1
    j=1
    
    while (i <= nrow(dat)) {
      # record the index started iteration in this while loop
      a=i

      # while the distance between the mid_win position of subsequent rows is <= than a window size
      # and while index is <= nrow of dat
      while ((dat$mid_win[i+1] - dat$mid_win[i]) <= (win_basepair*0.75) & (i+1 <= nrow(dat))) {
        i=i+1
      }

      # i is now the index of the last window consecutive that is below the treshold
      # save index to variable b
      b=i

      # for all windows between indexes a:b
      for (k in a:b){
        # create new column in data frame &
        # assign it trough ID number
        # all windows between indexes a:b will have the same ID number
        dat$t_id[k]=j
      }
      
      # move on to next window
      j=j+1
      i=i+1
    }
    
    return(dat)
  } else {
    # if there are no troughs below the treshold,
    # just returns NA
    return(NA)
    }
  
}

# df = |w_id |gen |cpop |rep |mean_pi_raw |mean_pi_samp |ini |end |wsize |mid_win |
get_summary_of_troughs = function (df, cur_gen, lvl, cur_rep, win_basepair){

  # calculates troughs, returns a data frame
  # tr_dat --> |w_id |gen |cpop |rep |mean_pi_raw |mean_pi_samp |ini |end |wsize |mid_win | t_id|
  tr_dat=calc_troughs(df, cur_gen, lvl, cur_rep, win_basepair)
  if (is.null(nrow(tr_dat))){
    return(NA)
  } else {
    # summarise data, first grouping by: t_id, cpop, rep, gen
    # number of windows: how many rows for each group
    # PI of trough: mean(pi of windows)
    # trough coordinates: ini of 1st win, end of last window
    # mid_win of trough: mean(mid_win of windows)
    tmp= ddply(tr_dat, .(t_id, cpop, rep, gen), summarise, n_windows=length(t_id), mp_raw=mean(mean_pi_raw), mp_samp=mean(mean_pi_samp), ini=min(ini), end=max(end), mid_tr=mean(mid_win))
    
    # trough size: end - ini
    tmp$size=tmp$end-tmp$ini

    # number of individual troughs
    tmp$tot_tr=nrow(tmp)
    
    return(tmp)
  }
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     median = median(xx[[col]], na.rm=na.rm),
                     qt2 = quantile(xx[[col]], na.rm=na.rm)[2],
                     qt3 = quantile(xx[[col]], na.rm=na.rm)[3],
                     qt4 = quantile(xx[[col]], na.rm=na.rm)[4],
                     qt1 = quantile(xx[[col]], probs = 0.025, na.rm=na.rm),
                     qt5 = quantile(xx[[col]], probs = 0.975, na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  datac <- rename(datac, c("qt2.25%" = "qt2"))
  datac <- rename(datac, c("qt3.50%" = "qt3"))
  datac <- rename(datac, c("qt4.75%" = "qt4"))
  datac <- rename(datac, c("qt5.97.5%" = "qt5"))
  datac <- rename(datac, c("qt1.2.5%" = "qt1"))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
