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