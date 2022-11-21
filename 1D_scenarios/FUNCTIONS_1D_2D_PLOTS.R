# -- -- -- -- -- 
# COLORS ####
# -- -- -- -- --


# level_trough = level_trough
# model_list = model_list
# measure_variable = mvar
# group_variable = gvar
# color_list = color_list
# win_size = ws
# f_trans = dividir
# resc_ymax = yscale
# plot_polig=F
# m_nbrs=c(1,1,1,2)
# m_height=1
# m_width=4
# main_lwd=4


plot_different_models_any_number=function(level_trough, model_list, win_size, measure_variable, group_variable, color_list, f_trans=1, resc_ymax=1, m_nbrs=NULL, m_height=NULL, m_width=NULL, plot_polig, main_lwd, polyg_alpha=NULL){
  
  if (is.null(m_nbrs)| is.null(m_height)|is.null(m_width)){
    m_height=1
    m_width=6
    m_nbrs=c(rep(1,5),2)
    l_matrix=matrix(m_nbrs, m_height, m_width)
  } else {
    l_matrix=matrix(m_nbrs, m_height, m_width)
  }
  
  if(is.null(polyg_alpha)){
    polyg_alpha=c(rep(0.01, length(model_list)))
  }
  
  # if (is.null(plot_polig)){
  #   plot_polig == T
  # } else { plot_polig = plot_polig}
  
  
  n_models=length(model_list)
  for (m in 1:n_models){
    file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_data_lvl_%s_*_*inds_win_%s.txt", level_trough*100, win_size)))
    assign(
      paste0("d", m),
      fread(sprintf("./r_%s/%s", model_list[m], file_name))
    )
  }
  
  
  for(m in 1:n_models){
    assign(
      paste0("color", m), unlist(color_list[m])
    )
  }
  
  for (m in 1:n_models){
    tmp_data=get(paste0("d", m))
    
    tryCatch({
      tmp_summary = summarySE(tmp_data, measurevar = measure_variable, groupvars = group_variable)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    assign(
      paste0("s", m),
      tmp_summary
    )
    
    
  }

  
  if (measure_variable == "tot_tr"){
    mv_label="Total number of troughs"
  } else if (measure_variable == "mean_size"){
    mv_label= "Trough mean size"
  } else if (measure_variable == "prop_tr"){
    mv_label="Proportion of the Chromosome in Troughs"
  } else {
    mv_label=measure_variable
  }
  
  if(f_trans != 1){
    mv_label = paste0(mv_label, " (scaled 1:", formatC(f_trans, format = "g"), ")")
  }
  
  
  if(group_variable == "gen"){
    gv_label="Time (generations)"
  } else if (group_variable== "cpop"){
    gv_label="Deme ID"
  } else {
    gv_label=group_variable
  }
  
  
  poss_y=c()
  min_y=c()
  for (m in 1:n_models){
    tmp_var=get(paste0("s", m))
    poss_y=c(poss_y, max(tmp_var$qt5))
    min_y=c(min_y, min(tmp_var$qt1))
  }
  
  poss_x=c()
  for (m in 1:n_models){
    tmp_var=get(paste0("s", m))
    poss_x=c(poss_x, max(tmp_var[group_variable]), min(tmp_var[group_variable]))
  }
  
  # if(min(poss_x)<25000){
  #   min_x=min(poss_x)
  #   max_x=max(poss_x)
  # } else {
  #   
  # }
  
  # max(poss_x)
  yli=max(poss_y)*resc_ymax
  min_yli=min(min_y)
  
    par(oma=c(0,2,0,0))
    # layout(matrix(c(1,1,1,1,1, 2), 1, 6))
    layout(l_matrix)
    mv_title=unlist(strsplit(mv_label, split = "\\("))[1]
    gv_title=unlist(strsplit(gv_label, split = " \\("))[1]
    plot(x=s1[,group_variable], y=s1[,measure_variable]/f_trans, type="l", col=color1, ylim = c(min_yli/f_trans, yli/f_trans), xlim=c(min(poss_x), max(poss_x)), main=sprintf("%s by %s;", mv_title, gv_title), xlab = gv_label, ylab = mv_label, cex.main=2, cex.lab=1.5, cex.axis=1.5, lwd=main_lwd)

    if (plot_polig == T){
      polygon(c(s1[,group_variable], rev(s1[,group_variable])), c(s1[,'qt5']/f_trans, rev(s1[,'qt1']/f_trans)), col=scales::alpha(color1, unlist(polyg_alpha[1])), border = scales::alpha(color1, 0.35))
    } 
    
    if (plot_polig == T){
      for (m in 2:n_models){
        sdat=get(paste0("s", m))
        acor=get(paste0("color", m))
        polygon(c(sdat[,group_variable], rev(sdat[,group_variable])), c(sdat[,'qt5']/f_trans, rev(sdat[,'qt1']/f_trans)), col=scales::alpha(acor, unlist(polyg_alpha[m])), border = scales::alpha(acor, 0.35))
        
      }
    }
      
     for (m in 2:n_models){
        sdat=get(paste0("s", m))
        acor=get(paste0("color", m))
        
        lines(x=sdat[,group_variable], y=sdat[,measure_variable]/f_trans, col=acor, lwd=main_lwd)
      }
      
      
    

    plot.new()
    par(oma=c(0,0,0,0))
    legend("right", legend = c(model_list), fill=c(unlist(color_list)), cex = 2, xpd = NA)
  
  
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


# Wed Jun 16 18:55:02 2021 ------------------------------


# win_size=ws;
# measure_variable="mean_size";
# group_variable="gen";
# f_trans=1; resc_ymax=1; main_lwd=4
# m=4
# model_list=model_list_rec
# color_list=color_list_rec
# polyg_alpha=polyg_alpha_rec
# lty_list=lty_list_rec

plot_different_models_NO_LAYOUT=function(level_trough, model_list, win_size, measure_variable, group_variable, color_list, f_trans=1, resc_ymax=1, plot_polig, main_lwd, polyg_alpha=NULL, cmain, clab, caxis, mgp_list, lty_list, set_y_max=NA, set_y_min=NA){
  
  if(is.null(polyg_alpha)){
    polyg_alpha=c(rep(0.01, length(model_list)))
  }

  
  n_models=length(model_list)
  for (m in 1:n_models){
    file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_*_lvl_%s_*_*_win_%s*", level_trough*100, win_size)))
    print(sprintf("%s: %s", m, model_list[m]))
    print(file_name)
    assign(
      paste0("d", m),
      fread(sprintf("./r_%s/%s", model_list[m], file_name))
    )
  }
  
  # n_models=length(model_list)
  # for (m in 1:n_models){
  # file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_data_lvl_%s_*_*inds_win_%s*", level_trough*100, win_size)))
  # sprintf("%s: %s", m, model_list[m])
  # cat(paste0(m, ": ", file_name), sep = "\n")
  # }
  
  for (m in 1:n_models){
    tmp = get(paste0("d", m))
  
    if (names(tmp)[1] != "cpop"){
      names(tmp) = c("cpop", "rep", "gen", "tot_win", "mean_size", "tot_tr", "chr_tr", "prop_tr")
      assign(paste0("d", m), tmp)
    }
  }
  
  for(m in 1:n_models){
    assign(
      paste0("color", m), unlist(color_list[m])
    )
  }
  
  for (m in 1:n_models){
    tmp_data=get(paste0("d", m))
    
    tryCatch({
      tmp_summary = summarySE(tmp_data, measurevar = measure_variable, groupvars = group_variable)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    assign(
      paste0("s", m),
      tmp_summary
    )
    
    
  }
  
  
  if (measure_variable == "tot_tr"){
    mv_label="Total number of troughs"
  } else if (measure_variable == "mean_size"){
    mv_label= "Trough mean size"
  } else if (measure_variable == "prop_tr"){
    mv_label="Proportion of the Chromosome in Troughs"
  } else {
    mv_label=measure_variable
  }
  
  if(f_trans != 1 & measure_variable != "mean_size"){
    mv_label = paste0(mv_label, " (scaled 1:", formatC(f_trans, format = "g"), ")")
  }
  
  if ((f_trans == 1e6) & (measure_variable=="mean_size")){
    mv_label=paste0(mv_label, " (Mb)")
  }
  
  
  if(group_variable == "gen"){
    gv_label="Time (generations)"
  } else if (group_variable== "cpop"){
    gv_label="Deme ID"
  } else {
    gv_label=group_variable
  }
  
  
  
  if (is.na(set_y_max)){
    poss_y=c()
    min_y=c()
    for (m in 1:n_models){
      tmp_var=get(paste0("s", m))
      poss_y=c(poss_y, max(tmp_var$qt5))
      min_y=c(min_y, min(tmp_var$qt1))
    }
    
    yli=max(poss_y)*resc_ymax
    min_yli=min(min_y)
    
  } else {
    yli=set_y_max
    min_yli=set_y_min
  }
  
  poss_x=c()
  for (m in 1:n_models){
    tmp_var=get(paste0("s", m))
    poss_x=c(poss_x, max(tmp_var[group_variable]), min(tmp_var[group_variable]))
  }
  


  
  # par(oma=c(0,2,0,0))
  # # layout(matrix(c(1,1,1,1,1, 2), 1, 6))
  # layout(l_matrix)
  mv_title=unlist(strsplit(mv_label, split = "\\("))[1]
  gv_title=unlist(strsplit(gv_label, split = " \\("))[1]
  plot(x=s1[,group_variable], y=s1[,measure_variable]/f_trans, type="l", col=color1, ylim = c(min_yli/f_trans, yli/f_trans), xlim=c(min(poss_x), max(poss_x)), main=sprintf("%s by %s;", mv_title, gv_title), xlab = gv_label, ylab = mv_label, cex.main=cmain, cex.lab=clab, cex.axis=caxis, lwd=main_lwd, mgp=c(mgp_list[1], mgp_list[2], mgp_list[3]), lty=lty_list[[1]])
  
  if (plot_polig == T){
    polygon(c(s1[,group_variable], rev(s1[,group_variable])), c(s1[,'qt5']/f_trans, rev(s1[,'qt1']/f_trans)), col=scales::alpha(color1, unlist(polyg_alpha[1])), border = scales::alpha(color1, 0.35), lty=lty_list[[1]])
  } 
  
  if (plot_polig == T){
    for (m in 2:n_models){
      sdat=get(paste0("s", m))
      acor=get(paste0("color", m))
      polygon(c(sdat[,group_variable], rev(sdat[,group_variable])), c(sdat[,'qt5']/f_trans, rev(sdat[,'qt1']/f_trans)), col=scales::alpha(acor, unlist(polyg_alpha[m])), border = scales::alpha(acor, 0.35), lty=lty_list[[m]])
      
    }
  }
  
  for (m in 2:n_models){
    sdat=get(paste0("s", m))
    acor=get(paste0("color", m))
    
    lines(x=sdat[,group_variable], y=sdat[,measure_variable]/f_trans, col=acor, lwd=main_lwd, lty=lty_list[[m]])
  }
  
  
  
  
  # plot.new()
  # par(oma=c(0,0,0,0))
  # legend("right", legend = c(model_list), fill=c(unlist(color_list)), cex = 2, xpd = NA)
  # 
  
}


legend_different_models_NO_LAYOUT=function(model_list, color_list, clegend, lty_list, lwd_list){
  legend("right", legend = c(model_list), fill=c(unlist(color_list)), cex = clegend, xpd = NA, lty = lty_list, lwd=lwd_list)
}

legend_LINES_models_NO_LAYOUT=function(model_list, color_list, clegend, lty_list=NA, lwd_list=NA, lpos=NA){
  if (is.na(lpos)){
    lpos="right"
  }
  
  if (is.na(lty_list)){
    lty_list=1
  }
  
  if (is.na(lwd_list)){
    lwd_list=2
  }
  
  legend(lpos, legend = c(model_list), fill=c(unlist(color_list)), col=c(unlist(color_list)), cex = clegend, xpd = NA, lty = lty_list, lwd=lwd_list)
}


# Tue Aug 31 00:19:14 2021 ------------------------------

# summ_type="simple"
# win_size=ws;
# measure_variable="mean_size";
# group_variable="gen";
# f_trans=1; resc_ymax=1; main_lwd=4

plot_recMap_models_NO_LAYOUT=function(level_trough, model_list, win_size, measure_variable, group_variable, color_list, f_trans=1, resc_ymax=1, plot_polig, main_lwd, polyg_alpha=NULL, cmain, clab, caxis, mgp_list, lty_list, summ_type){
  
  if(is.null(polyg_alpha)){
    polyg_alpha=c(rep(0.01, length(model_list)))
  }
  
  
  n_models=length(model_list)
  for (m in 1:n_models){
    #"r_1d_f10_m10_g5_recMap/summary_simple_recMap_lvl_10_f10_m10_g5_win_10k.txt"
    file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_%s_recMap_lvl_%s_*_win_%s*", summ_type, level_trough*100, win_size)))
    print(sprintf("%s: %s", m, model_list[m]))
    print(file_name)
    assign(
      paste0("d", m),
      fread(sprintf("./r_%s/%s", model_list[m], file_name))
    )
  }
  
  # n_models=length(model_list)
  # for (m in 1:n_models){
  # file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_data_lvl_%s_*_*inds_win_%s*", level_trough*100, win_size)))
  # sprintf("%s: %s", m, model_list[m])
  # cat(paste0(m, ": ", file_name), sep = "\n")
  # }
  
  for (m in 1:n_models){
    tmp = get(paste0("d", m))
    
    if (names(tmp)[1] != "cpop"){
      names(tmp) = c("cpop", "rep", "gen", "tot_win", "mean_size", "tot_tr", "chr_tr", "prop_tr")
      assign(paste0("d", m), tmp)
    }
  }
  
  for(m in 1:n_models){
    assign(
      paste0("color", m), unlist(color_list[m])
    )
  }
  
  for (m in 1:n_models){
    tmp_data=get(paste0("d", m))
    
    tryCatch({
      tmp_summary = summarySE(tmp_data, measurevar = measure_variable, groupvars = group_variable)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    assign(
      paste0("s", m),
      tmp_summary
    )
    
    
  }
  
  
  if (measure_variable == "tot_tr"){
    mv_label="Total number of troughs"
  } else if (measure_variable == "mean_size"){
    mv_label= "Trough mean size"
  } else if (measure_variable == "prop_tr"){
    mv_label="Proportion of the Chromosome in Troughs"
  } else {
    mv_label=measure_variable
  }
  
  if(f_trans != 1){
    mv_label = paste0(mv_label, " (scaled 1:", formatC(f_trans, format = "g"), ")")
  }
  
  
  if(group_variable == "gen"){
    gv_label="Time (generations)"
  } else if (group_variable== "cpop"){
    gv_label="Deme ID"
  } else {
    gv_label=group_variable
  }
  
  
  poss_y=c()
  min_y=c()
  for (m in 1:n_models){
    tmp_var=get(paste0("s", m))
    poss_y=c(poss_y, max(tmp_var$qt5))
    min_y=c(min_y, min(tmp_var$qt1))
  }
  
  poss_x=c()
  for (m in 1:n_models){
    tmp_var=get(paste0("s", m))
    poss_x=c(poss_x, max(tmp_var[group_variable]), min(tmp_var[group_variable]))
  }
  
  
  yli=max(poss_y)*resc_ymax
  min_yli=min(min_y)
  
  # par(oma=c(0,2,0,0))
  # # layout(matrix(c(1,1,1,1,1, 2), 1, 6))
  # layout(l_matrix)
  mv_title=unlist(strsplit(mv_label, split = "\\("))[1]
  gv_title=unlist(strsplit(gv_label, split = " \\("))[1]
  plot(x=s1[,group_variable], y=s1[,measure_variable]/f_trans, type="l", col=color1, ylim = c(min_yli/f_trans, yli/f_trans), xlim=c(min(poss_x), max(poss_x)), main=sprintf("%s by %s;", mv_title, gv_title), xlab = gv_label, ylab = mv_label, cex.main=cmain, cex.lab=clab, cex.axis=caxis, lwd=main_lwd, mgp=c(mgp_list[1], mgp_list[2], mgp_list[3]), lty=lty_list[[1]])
  
  if (plot_polig == T){
    polygon(c(s1[,group_variable], rev(s1[,group_variable])), c(s1[,'qt5']/f_trans, rev(s1[,'qt1']/f_trans)), col=scales::alpha(color1, unlist(polyg_alpha[1])), border = scales::alpha(color1, 0.35), lty=lty_list[[1]])
  } 
  
  if (plot_polig == T){
    for (m in 2:n_models){
      sdat=get(paste0("s", m))
      acor=get(paste0("color", m))
      polygon(c(sdat[,group_variable], rev(sdat[,group_variable])), c(sdat[,'qt5']/f_trans, rev(sdat[,'qt1']/f_trans)), col=scales::alpha(acor, unlist(polyg_alpha[m])), border = scales::alpha(acor, 0.35), lty=lty_list[[m]])
      
    }
  }
  
  for (m in 2:n_models){
    sdat=get(paste0("s", m))
    acor=get(paste0("color", m))
    
    lines(x=sdat[,group_variable], y=sdat[,measure_variable]/f_trans, col=acor, lwd=main_lwd, lty=lty_list[[m]])
  }
}





# measure_variable = mvar
# group_variable = gvar
# color_list =
# rec_list=c("low", "med")
# f_trans = dividir
# resc_ymax = yscale
# main_lwd=2

plot_recMap_diff_rates_NO_LAYOUT=function(level_trough, rec_list, model_list, measure_variable, group_variable, color_list, f_trans=1, resc_ymax=1, plot_polig, main_lwd, polyg_alpha=NULL, cmain=cmain, clab=clab, caxis=caxis, mgp_list = mgp_list, lty_list = lty_list, summ_type="simple"){
  
  if(is.null(polyg_alpha)){
    polyg_alpha=c(rep(0.01, length(window_list)))
  }
  
  
  n_models=length(model_list)
  for (m in 1:n_models){
    #"r_1d_f10_m10_g5_recMap/summary_simple_recMap_lvl_10_f10_m10_g5_win_10k.txt"
    file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_%s_recMap_lvl_%s_*_win_%s*", summ_type, level_trough*100, win_size)))
    print(sprintf("%s: %s", m, model_list[m]))
    print(file_name)
    
    tmp=fread(sprintf("./r_%s/%s", model_list[m], file_name))
    assign(
      paste0("d", m),
      fread(sprintf("./r_%s/%s", model_list[m], file_name))
    )
    
    med=tmp[tmp$major_rec=="med",]
    high=tmp[tmp$major_rec=="high",]
    low=tmp[tmp$major_rec=="low",]
    
    for (rec in rec_list){
      assign(
        paste0("d", m, "_", rec), get(paste0(rec)))
    }
    
    rm(med)
    rm(high)
    rm(low)
  }

  combinations=c()
  for (m in 1:n_models){
    for (rec in rec_list){
      print(paste0(m, "_", rec))
      combinations=c(combinations, paste0(m, "_", rec))
    }
  }
  
# combinations
  # "d1_low" "d1_med" "d2_low" "d2_med" "d3_low" "d3_med"
  
  # "d1_low"  "d1_med"  "d1_high" 
  # "d2_low"  "d2_med"  "d2_high" 
  # "d3_low"  "d3_med"  "d3_high"
  
  for(m in 1:length(combinations)){
    assign(
      paste0("color", m), unlist(color_list[m])
    )
  }
  
  for (m in 1:length(combinations)){
    print( sprintf("data: %s; summary: %s",   paste0("d", combinations[m]), paste0("s", combinations[m])))
    tmp_data=get(paste0("d", combinations[m]))
    print(tmp_data[1:2,])
    
    tryCatch({
      tmp_summary = summarySE(tmp_data, measurevar = measure_variable, groupvars = group_variable)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    assign(
      paste0("s", combinations[m]), tmp_summary )
    print(tmp_summary[1:2,1:5])
    rm(tmp_summary)
    rm(tmp_data)
  }
  
  
  if (measure_variable == "tot_tr"){
    mv_label="Total number of troughs"
  } else if (measure_variable == "mean_size"){
    mv_label= "Trough mean size"
  } else if (measure_variable == "prop_tr"){
    mv_label="Proportion of the Chromosome in Troughs"
  } else {
    mv_label=measure_variable
  }
  
  if(f_trans != 1){
    mv_label = paste0(mv_label, " (scaled 1:", formatC(f_trans, format = "g"), ")")
  }
  
  
  if(group_variable == "gen"){
    gv_label="Time (generations)"
  } else if (group_variable== "cpop"){
    gv_label="Deme ID"
  } else {
    gv_label=group_variable
  }
  
  
  poss_y=c()
  min_y=c()
  for (m in 1:length(combinations)){
    tmp_var=get(paste0("s", combinations[m]))
    poss_y=c(poss_y, max(tmp_var$qt5))
    min_y=c(min_y, min(tmp_var$qt1))
  }
  
  poss_x=c()
  for (m in 1:length(combinations)){
    tmp_var=get(paste0("s", combinations[m]))
    poss_x=c(poss_x, max(tmp_var[group_variable]), min(tmp_var[group_variable]))
  }
  

  yli=max(poss_y)*resc_ymax
  min_yli=min(min_y)
  
  mv_title=unlist(strsplit(mv_label, split = "\\("))[1]
  gv_title=unlist(strsplit(gv_label, split = " \\("))[1]
  plot(x=s1_low[,group_variable], y=s1_low[,measure_variable]/f_trans, type="l", col=color1, ylim = c(min_yli/f_trans, yli/f_trans), xlim=c(min(poss_x), max(poss_x)), main=sprintf("%s by %s;", mv_title, gv_title), xlab = gv_label, ylab = mv_label,
       cex.main=cmain, cex.lab=clab, cex.axis=caxis, lwd=main_lwd, mgp=c(mgp_list[1], mgp_list[2], mgp_list[3]), lty=lty_list[[1]])
  
  if (plot_polig == T){
    polygon(c(s1_low[,group_variable], rev(s1_low[,group_variable])), c(s1_low[,'qt5']/f_trans, rev(s1_low[,'qt1']/f_trans)), col=scales::alpha(color1, unlist(polyg_alpha[1])), border = scales::alpha(color1, 0.35))
  } 
  
  if (plot_polig == T){
    for (m in 2:length(combinations)){
      sdat=get(paste0("s", combinations[m]))
      acor=get(paste0("color", m))
      polygon(c(sdat[,group_variable], rev(sdat[,group_variable])), c(sdat[,'qt5']/f_trans, rev(sdat[,'qt1']/f_trans)), col=scales::alpha(acor, unlist(polyg_alpha[m])), border = scales::alpha(acor, 0.35))
      
    }
    
    for (m in 2:length(combinations)){
      print(combinations[m])
      sdat=get(paste0("s",  combinations[m]))
      acor=get(paste0("color", m))
      
      lines(x=sdat[,group_variable], y=sdat[,measure_variable]/f_trans, col=acor, lwd=main_lwd)
    }
    
    
  }
  
}


# lposition = "center"
# clegend=1
legend_recMap_diff_rates_NO_LAYOUT=function(model_list, color_list, lposition, clegend){
  n_models=length(model_list)
  combinations=c()
  for (m in 1:n_models){
    for (rec in rec_list){
      combinations=c(combinations, paste0(model_list[m], "_", rec))
    }
  }
  
  legend(lposition, legend = c(combinations), fill=c(unlist(color_list)), cex = clegend, xpd = NA, title = "Models")
}


# Wed Oct 13 23:22:52 2021 ------------------------------
# win_size="10k";
# measure_variable="mean_size";
# group_variable="ploss";
# f_trans=1; resc_ymax=1; main_lwd=4
# m=1
# model_list="f20_m10_g5"
# color_list="pink"
# pi_type="samp"
# level_trough=0.10
# set_y_max=NA

# aff=fread("r_f20_m10_g5/summary_data_lvl_10_f20_m10_g5_10inds_win_10k.txt.gz")

plot_models_DivLoss_NO_LAYOUT=function(level_trough, model_list, win_size, measure_variable, group_variable, color_list, f_trans=1, resc_ymax=1, alpha_list=NULL, cmain, clab, caxis, mgp_list, set_y_max=NA, set_y_min=NA, pi_type=NA, ini_div=NA, cex.symb){
  
  defined <- ls()
  passed <- names(as.list(match.call())[-1])
  
  if (any(!defined %in% passed)) {
    stop(paste("missing values for", paste(setdiff(defined, passed), collapse=", ")))
  }
  
  if(is.na(pi_type)){
    pi_type="samp"
  }
  
  if (is.na(ini_div)){
    ini_div=0.000125
  }
  
  
  n_models=length(model_list)
  for (m in 1:n_models){
    file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("diversity_data_%s_*_wsize_%s*", pi_type, win_size)))
    print(sprintf("%s: %s", m, model_list[m]))
    print(file_name)
    assign(
      paste0("d", m),
      fread(sprintf("./r_%s/%s", model_list[m], file_name))
    )
  }
  
  for (m in 1:n_models){
    file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_*_lvl_%s_*_*_win_%s*", level_trough*100, win_size)))
    print(sprintf("%s: %s", m, model_list[m]))
    print(file_name)
    assign(
      paste0("d", m, m),
      fread(sprintf("./r_%s/%s", model_list[m], file_name))
    )
  }
  
  
  for (m in 1:n_models){
    tmp = get(paste0("d", m))
    
    if (names(tmp)[1] != "cpop"){
      names(tmp) =c("gen",  "rep",  "N", "mean_pi_samp", "sd", "median", "qt2", "qt3",  "qt4",  "qt1", "qt5", "se", "ci")
      assign(paste0("d", m), tmp)
    }
  }
  
  
  for (m in 1:n_models){
    tmp = get(paste0("d", m, m))
    
    if (names(tmp)[1] != "cpop"){
      names(tmp) = c("cpop", "rep", "gen", "tot_win", "mean_size", "tot_tr", "chr_tr", "prop_tr")
      assign(paste0("d", m, m), tmp)
    }
  }
  
  for (m in 1:n_models){
    tmp1 = get(paste0("d", m))
    tmp2 = get(paste0("d", m, m))
  
    tmpmdat=merge(tmp1, tmp2, by= c("gen", "rep"))
    tmpmdat$pleft=tmpmdat[[4]]/ini_div
    tmpmdat$ploss=(1-tmpmdat$pleft)
    assign(paste0("s", m), tmpmdat)
  
  }
  
  for(m in 1:n_models){
    assign(
      paste0("color", m), unlist(color_list[m])
    )
  }
  
  
  if (measure_variable == "tot_tr"){
    mv_label="Total number of troughs"
  } else if (measure_variable == "mean_size"){
    mv_label= "Trough mean size"
  } else if (measure_variable == "prop_tr"){
    mv_label="Proportion of the Chromosome in Troughs"
  } else {
    mv_label=measure_variable
  }
  
  if(f_trans != 1 & measure_variable != "mean_size"){
    mv_label = paste0(mv_label, " (scaled 1:", formatC(f_trans, format = "g"), ")")
  }
  
  if ((f_trans == 1e6) & (measure_variable=="mean_size")){
    mv_label=paste0(mv_label, " (Mb)")
  }
  
  if ((f_trans == 1e3) & (measure_variable=="mean_size")){
    mv_label=paste0(mv_label, " (Kb)")
  }
  
  if(group_variable == "gen"){
    gv_label="Time (generations)"
  } else if (group_variable== "cpop"){
    gv_label="Deme ID"
  } else if (group_variable== "ploss"){
    gv_label="Proportion of diversity LOST"
  }else {
    gv_label=group_variable
  }
  
  
  
  if (is.na(set_y_max)){
    poss_y=c()
    min_y=c()
    for (m in 1:n_models){
      tmp_var=get(paste0("s", m))
      poss_y=c(poss_y, max(tmp_var[[measure_variable]]))
      min_y=c(min_y, min(tmp_var[[measure_variable]]))
    }
    
    yli=max(poss_y)*resc_ymax
    min_yli=min(min_y)
    
  } else {
    yli=set_y_max
    min_yli=set_y_min
  }
  
  poss_x=c()
  for (m in 1:n_models){
    tmp_var=get(paste0("s", m))
    poss_x=c(poss_x, max(tmp_var[[group_variable]]), min(tmp_var[[group_variable]]))
  }
  
  
  
  
  # par(oma=c(0,2,0,0))
  # # layout(matrix(c(1,1,1,1,1, 2), 1, 6))
  # layout(l_matrix)
  mv_title=unlist(strsplit(mv_label, split = "\\("))[1]
  gv_title=unlist(strsplit(gv_label, split = " \\("))[1]
  plot(x=s1[[group_variable]], y=s1[[measure_variable]]/f_trans, type="p",  cex=cex.symb, col=scales::alpha(color1, alpha_list[[1]]), ylim = c(min_yli/f_trans, yli/f_trans), xlim=c(min(poss_x), max(poss_x)), main=sprintf("%s by %s;", mv_title, gv_title), xlab = gv_label, ylab = mv_label, cex.main=cmain, cex.lab=clab, cex.axis=caxis, mgp=c(mgp_list[1], mgp_list[2], mgp_list[3]))
  
  # if (plot_polig == T){
  #   polygon(c(s1[,group_variable], rev(s1[,group_variable])), c(s1[['qt5']]/f_trans, rev(s1[['qt1']]/f_trans)), col=scales::alpha(color1, unlist(polyg_alpha[1])), border = scales::alpha(color1, 0.35), lty=lty_list[[1]])
  # } 
  
  # if (plot_polig == T){
  #   for (m in 2:n_models){
  #     sdat=get(paste0("s", m))
  #     acor=get(paste0("color", m))
  #     polygon(c(sdat[,group_variable], rev(sdat[,group_variable])), c(sdat[,'qt5']/f_trans, rev(sdat[,'qt1']/f_trans)), col=scales::alpha(acor, unlist(polyg_alpha[m])), border = scales::alpha(acor, 0.35), lty=lty_list[[m]])
      
  #   }
  # }
  
  for (m in 2:n_models){
    sdat=get(paste0("s", m))
    acor=get(paste0("color", m))
    
    points(x=sdat[[group_variable]], y=sdat[[measure_variable]]/f_trans, col=scales::alpha(acor, alpha_list[[m]]), cex=cex.symb)
  }
  
  
  
  
  # plot.new()
  # par(oma=c(0,0,0,0))
  # legend("right", legend = c(model_list), fill=c(unlist(color_list)), cex = 2, xpd = NA)
  # 
  
}


legend_different_models_NO_LAYOUT=function(model_list, color_list, clegend, lty_list, lwd_list){
  legend("right", legend = c(model_list), fill=c(unlist(color_list)), cex = clegend, xpd = NA, lty = lty_list, lwd=lwd_list)
}

legend_LINES_models_NO_LAYOUT=function(model_list, color_list, clegend, lty_list=NA, lwd_list=NA, lpos=NA){
  if (is.na(lpos)){
    lpos="right"
  }
  
  if (is.na(lty_list)){
    lty_list=1
  }
  
  if (is.na(lwd_list)){
    lwd_list=2
  }
  
  legend(lpos, legend = c(model_list), fill=c(unlist(color_list)), col=c(unlist(color_list)), cex = clegend, xpd = NA, lty = lty_list, lwd=lwd_list)
}

# Thu Oct 14 02:45:02 2021 ------------------------------



# win_size=ws; measure_variable=mvar; group_variable=gvar;
# color_list=color_list_levels; plot_polig=T;
# main_lwd=4; cmain=0.1; clab=0.1; caxis=0.1;
# lty_list=lty_list_levels; polyg_alpha=polyg_alpha_levels;
# model_list=model_list_levels; set_y_max=NA; f_trans=1;resc_ymax=1
  
plot_different_LEVELS_NO_LAYOUT=function(level_list, model_list, win_size, measure_variable, group_variable, color_list, f_trans=1, resc_ymax=1, plot_polig, main_lwd, polyg_alpha=NULL, cmain, clab, caxis, mgp_list, lty_list, set_y_max=NA, set_y_min=NA){
  
  if(is.null(polyg_alpha)){
    polyg_alpha=c(rep(0.01, length(model_list)))
  }
  
  
  n_models=length(model_list)
  for (m in 1:n_models){
    file_name=list.files(path= sprintf("./r_%s", model_list[m]), pattern = glob2rx(sprintf("summary_*_lvl_%s_*_*_win_%s*", level_list[[m]]*100, win_size)))
    print(sprintf("%s: %s, lvl %s", m, model_list[m], level_list[m]))
    print(file_name)
    assign(
      paste0("d", m),
      fread(sprintf("./r_%s/%s", model_list[m], file_name))
    )
  }

  
  for (m in 1:n_models){
    tmp = get(paste0("d", m))
    
    if (names(tmp)[1] != "cpop"){
      names(tmp) = c("cpop", "rep", "gen", "tot_win", "mean_size", "tot_tr", "chr_tr", "prop_tr")
      assign(paste0("d", m), tmp)
    }
  }
  
  for(m in 1:n_models){
    assign(
      paste0("color", m), unlist(color_list[m])
    )
  }
  
  for (m in 1:n_models){
    tmp_data=get(paste0("d", m))
    
    tryCatch({
      tmp_summary = summarySE(tmp_data, measurevar = measure_variable, groupvars = group_variable)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    assign(
      paste0("s", m),
      tmp_summary
    )
    
    
  }
  
  
  if (measure_variable == "tot_tr"){
    mv_label="Total number of troughs"
  } else if (measure_variable == "mean_size"){
    mv_label= "Trough mean size"
  } else if (measure_variable == "prop_tr"){
    mv_label="Proportion of the Chromosome in Troughs"
  } else {
    mv_label=measure_variable
  }
  
  if(f_trans != 1 & measure_variable != "mean_size"){
    mv_label = paste0(mv_label, " (scaled 1:", formatC(f_trans, format = "g"), ")")
  }
  
  if ((f_trans == 1e6) & (measure_variable=="mean_size")){
    mv_label=paste0(mv_label, " (Mb)")
  }
  
  if ((f_trans == 1e3) & (measure_variable=="mean_size")){
    mv_label=paste0(mv_label, " (Kb)")
  }
  
  if(group_variable == "gen"){
    gv_label="Time (generations)"
  } else if (group_variable== "cpop"){
    gv_label="Deme ID"
  } else {
    gv_label=group_variable
  }
  
  
  
  if (is.na(set_y_max)){
    poss_y=c()
    min_y=c()
    for (m in 1:n_models){
      tmp_var=get(paste0("s", m))
      poss_y=c(poss_y, max(tmp_var$qt5))
      min_y=c(min_y, min(tmp_var$qt1))
    }
    
    yli=max(poss_y)*resc_ymax
    min_yli=min(min_y)
    
  } else {
    yli=set_y_max
    min_yli=set_y_min
  }
  
  poss_x=c()
  for (m in 1:n_models){
    tmp_var=get(paste0("s", m))
    poss_x=c(poss_x, max(tmp_var[group_variable]), min(tmp_var[group_variable]))
  }
  
  
  mv_title=unlist(strsplit(mv_label, split = "\\("))[1]
  gv_title=unlist(strsplit(gv_label, split = " \\("))[1]
  plot(x=s1[,group_variable], y=s1[,measure_variable]/f_trans, type="l", col=color1, ylim = c(min_yli/f_trans, yli/f_trans), xlim=c(min(poss_x), max(poss_x)), main=sprintf("%s by %s;", mv_title, gv_title), xlab = gv_label, ylab = mv_label, cex.main=cmain, cex.lab=clab, cex.axis=caxis, lwd=main_lwd, mgp=c(mgp_list[1], mgp_list[2], mgp_list[3]), lty=lty_list[[1]])
  
  if (plot_polig == T){
    polygon(c(s1[,group_variable], rev(s1[,group_variable])), c(s1[,'qt5']/f_trans, rev(s1[,'qt1']/f_trans)), col=scales::alpha(color1, unlist(polyg_alpha[1])), border = scales::alpha(color1, 0.35), lty=lty_list[[1]])
  } 
  
  if (plot_polig == T){
    for (m in 2:n_models){
      sdat=get(paste0("s", m))
      acor=get(paste0("color", m))
      polygon(c(sdat[,group_variable], rev(sdat[,group_variable])), c(sdat[,'qt5']/f_trans, rev(sdat[,'qt1']/f_trans)), col=scales::alpha(acor, unlist(polyg_alpha[m])), border = scales::alpha(acor, 0.35), lty=lty_list[[m]])
      
    }
  }
  
  for (m in 2:n_models){
    sdat=get(paste0("s", m))
    acor=get(paste0("color", m))
    
    lines(x=sdat[,group_variable], y=sdat[,measure_variable]/f_trans, col=acor, lwd=main_lwd, lty=lty_list[[m]])
  }
  

  
}


legend_different_models_NO_LAYOUT=function(model_list, color_list, clegend, lty_list, lwd_list){
  legend("right", legend = c(model_list), fill=c(unlist(color_list)), cex = clegend, xpd = NA, lty = lty_list[1], lwd=lwd_list[1])
}
