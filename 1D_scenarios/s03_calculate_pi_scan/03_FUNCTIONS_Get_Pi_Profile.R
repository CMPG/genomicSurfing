# Wed Mar 31 13:13:21 2021 ------------------------------
## ---------------------------
##
## Script name: b14_FUNCTIONS_Get_Pi_Profile.R
##
## Purpose of script:
##
##
## Date Created: 2021-03-31
## ---------------------------

# Calculates Expected Heterozygosity
#Inputs: Frequency (p) and number of diploid individuals
get_pi_SAMPLING = function(freq, diplo){
  # sample size = 2n, because diploids
  samp=2*diplo
  # correction for samplin: (n/n-1)*Pi
  # pi: 1 - ((p^2) + (1-p)^2) --> Pi = 1 - (p^2 + q^2)
  (samp/(samp-1))*(1-((freq^2)+((1-freq)^2)))
}

# Exp. Het, without sampling correction
get_pi_RAW = function(freq){
  # pi: 1 - ((p^2) + (1-p)^2) --> pi = 1 - (p^2 + q^2)
 (1-((freq^2)+((1-freq)^2)))
}


# ini=win_coord$ini[1]
# fim=win_coord$end[1]
# wsize=win_coord$wsize[1]
# mid_win=win_coord$mid_win[1]
# id=win_coord$id[1]

# dt is the data used, contains allele freqs
# dt header: POS, size, allele_count, pop, gen, a_freq

subset_calc_pi_per_window = function(ini, fim, wsize, mid_win, id, dt){
  
  # gets all polymorphic sites inside the window coordinates into this temporary table
  tmp=dt[((dt$POS >= ini)&(dt$POS<fim)),]
  
  #calculates expected heterozygosity with function "get_pi", using derived allele frequency and pop size.
  # uses apply(df,1 ...) .--> goes through each row of the table
  # actually size divided by 2, because this version of output counts the number of alleles directly, not individuals...
  # and the function counts as if size == num.diplo.indviduals
  exp_het_samp=apply(tmp, 1, function(x) get_pi_SAMPLING(as.numeric(x['a_freq']), (as.numeric(x['size'])/2)))
  
  exp_het_raw=apply(tmp, 1, function(x) get_pi_RAW(as.numeric(x['a_freq'])))
  
  # passing resulting vector to the temp table
  tmp$pi_samp=exp_het_samp
  tmp$pi_raw=exp_het_raw
  
  
  
  #gets per bp values: sum all Pis and divides them by the window size
  mean_pi_samp=sum(tmp$pi_samp)/(wsize)
  mean_pi_raw=sum(tmp$pi_raw)/(wsize)
  
  return(data.frame(mean_pi_raw, mean_pi_samp, w_id=id))
  
}


#It just loops over all the mut files for all generations being analysed
# to be used in an apply "loop" of the "index" data frame
# cur_gen = fml$gen[2]
# cur_edge = fml$pop[2]

loop_to_get_Pi = function(cur_gen, cur_edge){
  a=Sys.time()
  library(rlist)
  gen=cur_gen
  edge=cur_edge
    
  print(paste0("GEN: ", gen, ", DEME: ", edge))
  
  dt=fread(sprintf("../vcf_files/count_size_g%s_p%s_r%s.txt.gz", gen, edge, rep))

  # POS, size, allele_count, pop, gen
  names(dt) = c("POS", "size", "allele_count", "pop", "gen")

  # allele freq: num. alleles / tot.number of alleles
  dt$a_freq = (dt$allele_count/(dt$size))
  # dt looks like this now:
  # POS, size, allele_count, pop, gen, a_freq


  # apply function: will subset genome in windows, according to the windows coordinate file
  # with dt, it will subset the loci within each window and then calculate pi
  calc_pis = apply(win_coord, 1, function (x) subset_calc_pi_per_window(ini = as.numeric(x['ini']), fim = as.numeric(x['end']), wsize = as.numeric(x['wsize']), mid_win = as.numeric(x['mid_win']), id = as.numeric(x['id']), dt))
  # each iteration returns: c(mean_pi_raw, mean_pi_samp, w_id)

  # transform returned list into table
  df_pis=list.rbind(calc_pis)

  # convert it to data.frame and add some info (gen, pop, rep...)
  res=data.frame(gen=gen, cpop=edge, rep=rep, mean_pi_raw=df_pis$mean_pi_raw, mean_pi_samp=df_pis$mean_pi_samp, w_id=df_pis$w_id)
  
  b=Sys.time()
  print(paste0("Done: g", gen, ", p", edge, " - ", format(b-a, digits=4)))
  return(res)

}


# format(Sys.time(), "%d_%b_%Hh_%Mm_%Ss")
