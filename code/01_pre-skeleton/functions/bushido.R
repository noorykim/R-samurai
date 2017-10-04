## Filename: bushido.r
## Author: Noory Kim
## Code for SAMURAI version: 1.2.0
## Developed using R version 3.0.0
## Last Updated: 8/22/2013

## References:
## BHHR: Borenstein, Hedges, Higgins & Rothstein (2009)
## CHV: Cooper, Hedges, & Valentine (2009)


#####################
## Tables: Level 0 ##
#####################


ExtractPublishedStudies <- function(table, colname="outlook"){  
  # Extracts the published studies (outlook == "published") from a prepared data set.
  #
  # Args: 
  #   table: The data set, in table form. 
  #   colname: The column containing the outlooks for the studies in the data set.
  #
  # Returns: The subset of published studies. 

  # Note:
  #   Use with(get()) to treat column names (strings) as variables:
  #   The with() statement looks for variables inside a dataframe.
  #   The get() statement accesses the object 'colname'.

  out <- table[which(with(table, get(colname)) == "published" ), ]
  # Rename column as 'outlook'
  names(out)[which(names(out) == colname)] <- "outlook" 
  return(out)  
}
# library(SAMURAI)
# data(greentea)
# load("greentea.rda")
# ExtractPublishedStudies(greentea)


ExtractUnpublishedStudies <- function(table, colname="outlook"){
  # Extracts the unpublished studies (outlook != "published") from a prepared data set.
  #
  # Args: 
  #   table: The data set, in table form. 
  #   colname: The column containing the outlooks for the studies in the data set.
  #
  # Returns: The subset of unpublished studies. 

  out <- table[which(with(table, get(colname)) != "published" ), ]
  ## Rename column as 'outlook'
  names(out)[which(names(out) == colname)] <- "outlook" 
  return(out)
}
# library(SAMURAI)
# data(greentea)
# ExtractUnpublishedStudies(greentea)


FillInMissingEffectSizeDefaults <- function(exhibit, vault){
  # Replace missing values in the 'exhibit' vector 
  # with corresponding entries in the 'vault' vector.
  #
  # Args: 
  #   exhibit: the vector of effect sizes to update.
  #   vault: the vector of effect sizes to draw from.
  #
  # Returns:  The exhibit vector with missing values replaced by 
  #           corresponding entries in the vault vector. 

  # Example: 
  #   u <- c(2,NA)
  #   c <- c(5,3)
  #   FillInMissingEffectSizeDefaults(c,u) 
  #   > 2 3

  # Note: 
  #  There is no easy way to have elementwise addition of vectors while ignoring NA's.
  #  Setting NA's to zeros helped workaround this issue.

  vault.miss <- is.na(vault) ## missingness indicators
  retain <- exhibit * vault.miss
  vault[which(is.na(vault))] <- 0 ## change NA to zeros
  out <- vault + retain
  return(out)
}
# u <- c(2,NA)
# c <- c(5,3)
# FillInMissingEffectSizeDefaults(c,u) 


CompleteOutlooksFactor <- function(table){
  # When importing a data set, R will include in the 'outlooks' vector only those outlooks found in the data set. 
  # This function completes the 'outlooks' vector with other outlooks not found in the data set,
  # so that any outlook may be summoned when imputing unpublished studies with a different outlook. 
  #
  # Args: 
  #   table: The data set with an 'outlook' column. 
  #
  # Returns: The same data set with the factor 'outlook' complete.

  # Example:
  #   greentea$outlook
  #   > Levels: no effect positive published
  #   greentea <- CompleteOutlooksFactor(greentea)
  #   greentea$outlook
  #   > 11 Levels: no effect positive published very positive ... very negative CL 

  # Dependencies:
  #   Callers: AssignSameRustlook()

  outlooks <- c("very positive", "positive", 
                "no effect", 
                "negative", "very negative", 
                "very positive CL", "positive CL", 
                "current effect", 
                "negative CL", "very negative CL")
  levels(table$outlook) <- union(levels(table$outlook),outlooks)
  return(table)
}
# library(SAMURAI)
# data(greentea)
# greentea$outlook
# greentea <- CompleteOutlooksFactor(greentea)
# greentea$outlook


AssignSameRustlook <- function(table,rustlook){
  # Assigns the same outlook ('rustlook') to all unpublished studies in the data set.
  # The assigned rustlook will override any previously assigned outlook. 
  #
  # Args: 
  #   table: The data set.
  #   rustlook: The outlook to be assigned to all unpublished studies.
  #
  # Returns: The data set with the rustlook assigned to all unpublished studies.

  # Examples:
  #   AssignSameRustlook(Hpylori,"negative")
  #   AssignSameRustlook(Hpylori,"current effect")

  # Dependencies:
  #   Calls: CompleteOutlooksFactor()

  table <- CompleteOutlooksFactor(table)
  table$outlook[which(table$outlook != "published")] <- rustlook
  return(table)
}
# library(SAMURAI)
# data(greentea)
# AssignSameRustlook(greentea, "negative")


ConvertBinaryToEffectSize <- function(table, measure="RR"){
  # Converts count data to an effect size (such as log relative risk).  
  #
  # Args:
  #   table: A data set table with binary count data for both the ctrl and expt arms.
  #   measure: which effect size measure to calculate
  #     "RR" = log risk ratio (default)
  #     "OR" = log odds ratio 
  #
  # Returns: Adjusted SMD/Hedges g and its variance under a fixed effects model.
  
  # Dependencies:
  #   Calls: package "metafor" to use function escalc()
  #   Callers: binary.table()

  fixedmodel <- escalc(measure=measure, data=table, 
                  append=TRUE,
                  ai=table$expt.events, n1i=table$expt.n, 
                  ci=table$ctrl.events, n2i=table$ctrl.n)
  return(fixedmodel)
}
# library(SAMURAI)
# data(Hpylori)
# Hpylori.pub <- ExtractPublishedStudies(Hpylori)
# ConvertBinaryToEffectSize(Hpylori.pub)


ConvertMeanSDToSMD <- function(table, measure="SMD"){
  # Converts mean/standard deviation (sd) data to standardized mean differences (SMD).  
  #
  # Args:
  #   table: A data set table with continuous data (mean,sd,n) for both ctrl and expt arms.
  #   measure: which effect size measure to calculate
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #
  # Returns: Adjusted SMD/Hedges g (yi) and its variance (vi) under a fixed effects model.

  # Dependencies: 
  #   Calls: package "metafor" to use function escalc()
  #   Callers: contintuousforest(), binary.table()

  fixedmodel <- escalc(measure=measure, data=table, 
                  append=TRUE,
                  m1i=table$expt.mean, sd1i=table$expt.sd, n1i=table$expt.n,
                  m2i=table$ctrl.mean, sd2i=table$ctrl.sd, n2i=table$ctrl.n)  
  return(fixedmodel)
}
# library(SAMURAI)
# data(greentea)
# greentea.pub <- ExtractPublishedStudies(greentea)
# ConvertMeanSDToSMD(greentea.pub)


CalculateSummaryEffect <- function(table, level=95,
  summary.measure="SMD", method="DL"){
  # Compute a random-effects summary effect (cSMD) from a fixed-effects model with Hedges g
  #
  # Args:
  #   table: Table of Hedges g and their variances under a fixed effects model
  #   summary.measure: 
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #   method: "DL" for the DerSimonian & Laird method (1996) (default)
  #   level: confidence level = 1 - alpha
  #
  # Returns: Random effects model meta-analytic summary

  # Dependencies: 
  #   Calls: package 'metafor' to use the function rma()
  
  # Notes: 
  #   The ... in this function acts as a 'garbage collector' for runaway parameters upstream.
  
  # to avoid R CMD CHECK NOTE: "no visible binding for global variable"
  yi <- NULL 
  vi <- NULL 
  
  randmodel <- rma(yi, vi, data=table, measure=summary.measure, method=method, level=level)  
  
  # binary: summary log RR, or summary log OR
  # continuous: effect size SMD
  m <- randmodel$b[1]  
  m.se <- randmodel$se[1]
  m.lcl <- randmodel$ci.lb[1]  
  m.ucl <- randmodel$ci.ub[1]
  
  # binary: summary RR, or summary OR
  # continuous: effect size SMD
  expm <- exp(m) 
  expm.lcl <- exp(randmodel$ci.lb[1]) 
  expm.ucl <- exp(randmodel$ci.ub[1]) 
  
  # measures of heterogeneity
  tau2  <- randmodel$tau2[1]
  Q     <- randmodel$QE[1]  
  Qpval <- randmodel$QEp[1]

  out <- as.list(c(m, m.se, m.lcl, m.ucl, expm.lcl, expm, expm.ucl,tau2, Q, Qpval))
  names(out) <- c("m","m.se","m.lcl","m.ucl", "exp.m.lcl","exp.m","exp.m.ucl","tau2","Q","Qpval")    
  return(out)
  # return(list(m=m, m.se=m.se, m.lcl=m.lcl, m.ucl=m.ucl, 
  #             expm.lcl=expm.lcl, expm=expm, expm.ucl=expm.ucl,
  #             tau2=tau2, Q=Q, Qpval=Qpval))
}
# library(SAMURAI)
# data(Hpylori)
# Hpylori.pub <- ConvertBinaryToEffectSize(ExtractPublishedStudies(Hpylori))
# Hpylori.pub.summ <- CalculateSummaryEffect(Hpylori.pub)
# Hpylori.pub.summ
# as.numeric(Hpylori.pub.summ[c(3,1,4)])


SetEffectSizesIndependentOfPub <- function(
  binary=NA, higher.is.better=NA, 
  vpos=NA, pos=NA, neg=NA, vneg=NA){
  # Assigns effects that are not based on CI of summary effect across published studies.
  #
  # Args:
  #   binary: True/False - The data is binary/count. 
  #   higher.is.better: T/R - Higher counts/effect sizes are desired.
  #   vpos: The effect size to assign studies with a "very positive" outlook.
  #   pos:  The effect size to assign studies with a "positive" outlook.
  #   neg:  The effect size to assign studies with a "negative" outlook.
  #   vneg: The effect size to assign studies with a "very negative" outlook.
  #
  # Returns: A vector of effect sizes. 

  # Dependencies:
  #   Calls: FillInMissingEffectSizeDefaults()

  # Notes:
  #   This function contains lists of preset default values. 
  
  ## Algorithm defaults
  effects.default.rr <- c(3,2,1,1/2,1/3)
  effects.default.smd <- c(0.8,0.3,0,-0.3,-0.8)
  if(binary==TRUE){
    if(higher.is.better==TRUE){
      effects.default <- effects.default.rr
    } else{
      effects.default <- 1/effects.default.rr
    }  
  } else{
    if(higher.is.better==TRUE){
      effects.default <- effects.default.smd
    } else{
      effects.default <- effects.default.smd*(-1)
    }
  }

  ## User definitions
  effects.user <- c(vpos,pos,NA,neg,vneg)

  # Override algorithm defaults with user definitions
  effect.sizes.list <- FillInMissingEffectSizeDefaults(effects.default,effects.user)  
  return(effect.sizes.list)
}
# SetEffectSizesIndependentOfPub(binary=TRUE, higher.is.better=FALSE)


SetEffectSizesDependentOnPub <- function(table, 
  binary=NA, mean.sd=NA, 
  higher.is.better=NA, 
  level=95, 
  binary.measure="RR", continuous.measure="SMD", 
  summary.measure="SMD", method="DL"){
  # Assigns effects that are based on the CI of the summary effect across published studies.
  #
  # Args:
  #   table: the data set
  #   binary: True/False - The data is binary/count. 
  #   mean.sd: T/F - The data is in the form of means and standard deviations. 
  #   higher.is.better: T/R - Higher counts/effect sizes are desired.
  #   binary.measure:     if(binary == TRUE), 
  #                       which effect size measure do we use for the conversion?
  #     "RR" = relative risk
  #     "OR" = odds ratio 
  #   continuous.measure: if(binary == FALSE && mean.sd == TRUE), 
  #                       which effect size measure do we use for the conversion? 
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #   summary.measure: Which effect size measure do we use for the summary effect?
  #   method: "DL" for the DerSimonian & Laird method (1996) (default)
  #   level: confidence level = 1 - alpha
  # Returns: A list of effect sizes. 

  # Dependencies:
  #   Calls: FillInMissingEffectSizeDefaults(), CalculateSummaryEffect()
  #     ExtractPublishedStudies(), ConvertBinaryToEffectSize(), ConvertMeanSDToSMD(),
  
  ## Extract published studies
  pub <- ExtractPublishedStudies(table)
  
  ## Convert count data to log RR and its variance for each study
  if(binary==TRUE){
    pub <- ConvertBinaryToEffectSize(pub, measure=binary.measure) 
  } else if(mean.sd==TRUE){
    pub <- ConvertMeanSDToSMD(pub, measure=continuous.measure) 
  } 
  
  ## Calculate summary effect across all published studies
  pub.summary <- CalculateSummaryEffect(pub, summary.measure=summary.measure, 
                                        method=method, level=level) 
  
  if(binary == TRUE){
    pub.effect <- pub.summary$exp.m
    pub.effect.lcl <- pub.summary$exp.m.lcl
    pub.effect.ucl <- pub.summary$exp.m.ucl
  } else {
    pub.effect <- pub.summary$m
    pub.effect.lcl <- pub.summary$m.lcl
    pub.effect.ucl <- pub.summary$m.ucl
  }
  
  # Calculate 
  halfdown <- 0.5 * (pub.effect + pub.effect.lcl)
  halfup   <- 0.5 * (pub.effect + pub.effect.ucl)    
  if(higher.is.better == TRUE){
    vposcl <- pub.effect.ucl
    poscl  <- halfup
    negcl  <- halfdown
    vnegcl <- pub.effect.lcl    
  } else {
    vposcl <- pub.effect.lcl
    poscl  <- halfdown
    negcl  <- halfup
    vnegcl <- pub.effect.ucl        
  }
  current <- pub.effect
  effect.sizes.list.pub.ci <- c(vposcl,poscl,current,negcl,vnegcl)
  
  return(effect.sizes.list.pub.ci)
}
# SetEffectSizesDependentOnPub(table=Hpylori,binary=T,higher.is.better=F)


CompileEffectsToBeAssigned <- function(table, 
  binary=TRUE, mean.sd=FALSE, 
  higher.is.better=TRUE, 
  vpos=NA, pos=NA, neg=NA, vneg=NA, 
  level=95,
  binary.measure="RR", continuous.measure="SMD", 
  summary.measure="SMD", method="DL"){
  # Assigns both (1) effect sizes based on the CI of the summary effect across published studies, 
  # and (2) effect sizes not based that CI. 
  #
  # Args:
  #   table: the data set
  #   binary: True/False - The data is binary/count. 
  #   higher.is.better: T/R - Higher counts/effect sizes are desired.
  #   mean.sd: T/F - The data is in the form of means and standard deviations. 
  #   vpos: The effect size to assign studies with a "very positive" outlook.
  #   pos:  The effect size to assign studies with a "positive" outlook.
  #   neg:  The effect size to assign studies with a "negative" outlook.
  #   vneg: The effect size to assign studies with a "very negative" outlook.
  #   binary.measure: 
  #   continuous.measure: The 
  #     "SMD" = standardized mean difference (default)
  #     "SMDH" = standardized mean difference w/o assuming equal population variances 
  #              in the two groups 
  #
  # Returns: A list of effect sizes. 

  # Dependencies:
  #   Calls: ExtractPublishedStudies(), SetEffectSizesIndependentOfPub(), SetEffectSizesDependentOnPub()
  
  ## Extract published studies
  pub <- ExtractPublishedStudies(table)
  
  ## Assign effects (not based on CI of summary effect across published studies)
  eff.nopub <- SetEffectSizesIndependentOfPub(binary=binary, 
                higher.is.better=higher.is.better,
                vpos=vpos, pos=pos, neg=neg, vneg=vneg)
  ## Assign effects based on CI of summary effect across published studies)
  eff.onpub <- SetEffectSizesDependentOnPub(table, 
                binary=binary, mean.sd=mean.sd, 
                higher.is.better=higher.is.better, 
                level=level, 
                binary.measure=binary.measure, continuous.measure=continuous.measure,
                summary.measure=summary.measure, method=method)
  effect.sizes.list <- as.list(c(eff.nopub, eff.onpub))
  names(effect.sizes.list) <- c("vpos", "pos", "noef", "neg", "vneg",
    "vposcl", "poscl", "curr", "negcl", "vnegcl")
  return(effect.sizes.list)
}
# library(SAMURAI)
# data(Hpylori)
# CompileEffectsToBeAssigned(table=Hpylori,
#   binary=TRUE, higher.is.better=TRUE, level=99)
# data(greentea)
# CompileEffectsToBeAssigned(table=greentea, binary=FALSE, mean.sd=TRUE, 
#   higher.is.better=TRUE, level=95)


ImputeControlGroupEvents <- function(table, rounded=FALSE){
  # Calculates the summary relative risk across published studies, 
  # then uses that to impute the number of events in the control arms of unpublished studies
  # (assuming the rate of events in control arms of unpubs are the same as the rate across pubs).
  # 
  # Args: 
  #   table:  A data set with event counts and sample sizes for both arms of published studies
  #           and with sample sizes for both arms of unpublished studies. 
  #   rounded: T/F - Imputed counts are rounded to the nearest integer. 
  #
  # Returns:  The same data set with imputed event counts in control arms of unpublished studies.

  pub                 <- ExtractPublishedStudies(table)
  sum.pub.ctrl.n      <- sum(pub$ctrl.n)
  sum.pub.ctrl.events <- sum(pub$ctrl.events)
  pub.ctrl.propn      <- sum.pub.ctrl.events / sum.pub.ctrl.n

  unpub               <- ExtractUnpublishedStudies(table)
  unpub$ctrl.events   <- pub.ctrl.propn * unpub$ctrl.n    

  if(rounded == TRUE){
    unpub$ctrl.events <- round(unpub$ctrl.events)
  }
  out <- rbind(pub, unpub)
  return(out)
}
# library(SAMURAI)
# data(Hpylori)
# ImputeControlGroupEvents(Hpylori)


ImputeBinaryEvents <- function(table, effect.sizes.list, seed=NA, sims=1){
  # Imputes the number of events in the intervention arms of unpublished studies,
  # depending on the outlook(s) of those studies. 
  # 
  # Args: 
  #   table:  A data set with the following information:
  #             published studies, both arms: event counts, sample sizes
  #             unpublished studies:
  #               outlooks 
  #               control arm:      sample sizes, event counts
  #               intervention arm: sample sizes
  #   effect.sizes.list: A list of effect sizes to be assigned depending on outlook.  
  #   seed: random number seed
  #   sims: number of simulations (over which to average out the summary effect size)
  #
  # Returns:  The same data set with imputed event counts in intervention arms of unpublished studies.

  set.seed(NULL)  ## reset seed - redundant if done within a function call?
  if(!is.na(seed)){
    set.seed(seed)  
  }  

  # table <- Hpylori; sims <- 1       # testing

  ## Impute ctrl.events = proportion of events across published studies
  pub <- ExtractPublishedStudies(table)
  ctrl.propn <- sum(pub$ctrl.events) / sum(pub$ctrl.n)

  table <- ImputeControlGroupEvents(table)
  table <- CompleteOutlooksFactor(table)

  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  varnames <- c("vpos", "pos", "noef", "neg", "vneg", "vposcl", "poscl", "curr", "negcl", "vnegcl")

  # The number of different outlooks in the list above. 
  n <- length(outlooks)
   
  ## go through each of the outlooks in the list 
  for(i in 1:n){  
   
    # i <- 1  # testing

    # does any unpub in data set have that particular outlook?
    if(outlooks[i] %in% table$outlook){  
      
      # how many unpub studies have that outlook? 
      k <- length(which(table$outlook == outlooks[i]))  

      # vector of sample sizes of intervention arms
      expt.n <- table[which(table$outlook == outlooks[i]), ]$expt.n  

      # retrieve effect size for the outlook 
      # and multiply by the proportion of events across control arms
      # to impute the rate of events in the intervention arms
      expt.propn <- with(effect.sizes.list, get(varnames[i])) * ctrl.propn  
      
      # generate random numbers based on binomial distribution 
      # variable 'sims' determines how many simulations to be averaged over 
      table[which(table$outlook == outlooks[i]), ]$expt.events <- rbinom(n=k, size=expt.n*sims, prob=expt.propn) / sims  ## vector
    }  
  }  

  # Round to nearest integer
  table$ctrl.events <- round(table$ctrl.events)
  table$expt.events <- round(table$expt.events)
  # Avoid having # events exceed arm sample size # 
  table$ctrl.events <- pmin(table$ctrl.events, table$ctrl.n)
  table$expt.events <- pmin(table$expt.events, table$expt.n)
  # Avoid having # events negative
  table$ctrl.events <- pmax(table$ctrl.events, 0)
  table$expt.events <- pmax(table$expt.events, 0)
  
  return(table)
}
# library(SAMURAI)
# data(Hpylori)
# effect.sizes.list <- CompileEffectsToBeAssigned(
#   table=Hpylori, binary=TRUE, higher.is.better=TRUE, vpos=NA, pos=NA, neg=NA, vneg=NA,
#   mean.sd=FALSE, binary.measure="RR", continuous.measure="SMD", summary.measure="SMD",
#   method="DL", level=95)
# ImputeBinaryEvents(Hpylori,effect.sizes.list,seed=NA, sims=1)


ImputeSMD <- function(table, effect.sizes.list, seed=NA, noise=0.01){
  # Imputes the standardized mean difference of unpublished studies,
  # depending on the outlook(s) of those studies. 
  # 
  # Args: 
  #   table:  A data set with the following information:
  #             published studies: means, sd's, sample sizes
  #             unpublished studies: outlooks, sample sizes
  #   effect.sizes.list: A list of effect sizes to be assigned depending on outlook.  
  #   seed: random number seed
  #   noise: user added Gaussian random noise; standard deviation from the assigned effect size 
  #
  # Returns:  The same data set with imputed effect sizes for unpublished studies.

  set.seed(NULL)  ## reset seed - redundant if done within a function call?
  if(!is.na(seed)){
    set.seed(seed)  
  }  

  table <- CompleteOutlooksFactor(table)
  
  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  varnames <- c("vpos", "pos", "noef", "neg", "vneg", "vposcl", "poscl", "curr", "negcl", "vnegcl")
  n <- length(outlooks)
  
  for(i in 1:n){
    # i <- 4
    if(outlooks[i] %in% table$outlook){ # if there unpub assigned to that outlook
      k <- length(which(table$outlook == outlooks[i]))
      # retrieve smd to be assigned
      smd.assigned <- with(effect.sizes.list, get(varnames[i]))
      # add some random noise
      table[which(table$outlook==outlooks[i]), ]$yi <- rnorm(k, mean=smd.assigned, sd=noise)
    }  
  }  
  
  return(table)
}
# library(SAMURAI)
# data(greentea)
# effect.sizes.list <- CompileEffectsToBeAssigned(table=greentea,
#                       binary=F, mean.sd=T, 
#                       higher.is.better=F, level=95)
# gt <- ConvertMeanSDToSMD(greentea)
# ImputeSMD(gt, effect.sizes.list, seed=NA, noise=0.01)


ImputeSMDVariance <- function(table
                              # , matchedgroups=FALSE
                         ){
  # Imputes the variance of the standardized mean difference (SMD) of unpublished studies,
  # depending on the SMD and the sample sizes of those studies. 
  # 
  # Args: 
  #   table:  A data set with the following information: 
  #             sample sizes for both arms 
  #             SMD's (Hedges' g, which is a sample estimate of the SMD)
  #
  # Returns:  The same data set with imputed variance of Hedges' g of unpublished studies.
  #           It estimates variance of Hedge's g using a "very good" approximation by Borenstein. 
  #
  # Reference: 
  #   Michael Borenstein, "Effect Sizes for Continuous Data", page 226,
  #   Chapter 12 in Cooper, Hedges, and Valentine, Handbook of Research Synthesis and Meta-analysis
  
  ## testing
  #   table <- table4
  
  table$flagmissing <- as.numeric(is.na(table$vi))
  
  n1 <- table$ctrl.n
  n2 <- table$expt.n
  
  df <- n1+n2-2
  j <- 1 - 3/(4*df-1)  # correction factor between Cohen's d and Hedge's g
  
  ## convert to Cohen's d
  hedgesg <- table$yi
  cohensd <- hedgesg/j
  ## find variance of Cohen's d; convert to variance of Hedges' g
  cohensd.v <- (n1+n2)/(n1*n2) + cohensd^2/(2*(n1+n2))
  table$hedgesg.v <- j^2 * cohensd.v
  ##
  table[table$flagmissing==1,]$vi <- table[table$flagmissing==1,]$hedgesg.v  
  ## drop columns
  table$flagmissing <- NULL
  table$hedgesg.v <- NULL
  return(table)
}
# library(SAMURAI)
# data(greentea)
# effect.sizes.list <- CompileEffectsToBeAssigned(table=greentea,
#                       binary=F, mean.sd=T, 
#                       higher.is.better=F, level=95)
# gt <- ConvertMeanSDToSMD(greentea)
# gt <- ImputeSMD(gt,effect.sizes.list,seed=NA, noise=0.01)
# gt <- ImputeSMDVariance(gt)
# gt


#####################
## Tables: Level 1 ##
#####################


PrepareTableWithBinaryData <- function(table,
  higher.is.better=NA, 
  rustlook=NA,
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95, 
  binary.measure="RR", summary.measure="SMD", method="DL", 
  seed=NA, sims=1){
  # Input a table and impute all events for each study, then calculate Hedges' g and its variance.
  # 
  # Args:
  #
  # Returns: A table with all events imputed, and with Hedges' g and its variance calculated.

  ## testing
  #   binary=TRUE; binary.measure="RR"; continuous.measure="SMD"; mean.sd=FALSE; higher.is.better=TRUE;
  #   vpos=NA; pos=NA; neg=NA; vneg=NA; rustlook=NA; level=95; method="DL"
  
  ## Avoid R CMD CHECK NOTE: "no visible binding for global variable".
  expt.events <- expt.n <- ctrl.events <- ctrl.n <- NULL    
  
  ## Assign effects (not based on CI of summary effect across published studies)
  effect.sizes.list <- CompileEffectsToBeAssigned(table=table,
                        binary=TRUE, 
                        higher.is.better=higher.is.better,
                        vpos=vpos, pos=pos, neg=neg, vneg=vneg,
                        binary.measure=binary.measure, 
                        summary.measure=summary.measure, method=method, 
                        level=level)
    
  ## If desired, assign all unpublished studies the same outlook
  if(!is.na(rustlook)){
    table <- AssignSameRustlook(table,rustlook)
    #     levels(table$outlook)
  }  

  table <- ImputeBinaryEvents(table,effect.sizes.list,sims=sims)
  table <- ConvertBinaryToEffectSize(table,measure=binary.measure)

  return(table)  
}
# library(SAMURAI)
# data(Hpylori)
# library(metafor)
# load("Hpylori.rda")
# effect.sizes.list <- CompileEffectsToBeAssigned(table=Hpylori, 
#                       binary=TRUE, mean.sd=FALSE,
#                       higher.is.better=TRUE, 
#                       level=95,
#                       binary.measure="RR", continuous.measure="SMD", 
#                       summary.measure="SMD", method="DL")
# hp <- ImputeBinaryEvents(table=Hpylori, effect.sizes.list, seed=NA, sims=1)
# PrepareTableWithBinaryData(table=hp, higher.is.better=F, 
#   rustlook="very negative")
# PrepareTableWithBinaryData(table=hp, higher.is.better=F, 
#   rustlook="very negative", seed=52)


PrepareTableWithContinuousData <- function(table,
  mean.sd=TRUE, 
  higher.is.better=TRUE, 
  rustlook=NA, 
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95, 
  continuous.measure="SMD", summary.measure="SMD", method="DL", 
  seed=NA, noise=0.01){ 
  # Input a table and impute Hedges' g and its variance for each study.
  # 
  # Args:
  #
  # Returns: A table with all events imputed, and with Hedges' g and its variance calculated.

  ## testing
  #   table <- greentea
  #   binary=F; continuous.measure="SMD"; mean.sd=TRUE; higher.is.better=TRUE;
  #   vpos=NA; pos=NA; neg=NA; vneg=NA; level=95; method="DL"
  #   rustlook="negative"
  
  if(mean.sd==TRUE){
    table <- ConvertMeanSDToSMD(table)  
  }
  
  ## Assign effects (not based on CI of summary effect across published studies)
  effect.sizes.list <- CompileEffectsToBeAssigned(table=table,
                          binary=FALSE,  # since we have continuous data
                          higher.is.better=higher.is.better,
                          vpos=vpos, pos=pos,neg=neg,vneg=vneg,
                          mean.sd=mean.sd, 
                          binary.measure=NA, continuous.measure=continuous.measure,
                          summary.measure=summary.measure,
                          method=method, level=level)                                        
  
  ## If desired, assign all unpublished studies the same outlook
  if(!is.na(rustlook)){
    table <- AssignSameRustlook(table,rustlook)
#     levels(table$outlook)
  }  

  table <- ImputeSMD(table, effect.sizes.list,
                     seed=seed, noise=noise)
  table <- ImputeSMDVariance(table)
  
  return(table)  
}
# library(SAMURAI)
# data(greentea)
# load("greentea.rda")
# PrepareTableWithContinuousData(table=greentea, higher.is.better=F, 
#                       rustlook="negative",noise=0.01)


CalculateSummaryEffectsForOneTable <- function(table, measure="RR", method="DL", level=95){
  # Given a table, calculate the aggregate estimates for each of the following subsets:
  #   published
  #   unpublished
  #   all (published & unpublished)
  # 
  # Args:
  #   table:  A data set with missing figures already imputed 
  #           and with effect sizes already calculated for each study.
  #
  # Returns a table with 3 rows / summary effects: across pubs, across unpubs, across all.
  #
  # Note: output is formatted for the addpoly() function in the metafor package
  
  # published
  pub <- table[which(table$outlook=="published"),]
  if(nrow(pub)>0){
    pub.agg <- CalculateSummaryEffect(pub, summary.measure="SMD", method=method, level=level) 
  } else { pub.agg <- NA}
  
  # summary effect for unpublished studies only
  unpub <- table[which(table$outlook != "published"),]
  if(nrow(unpub)>0){
    unpub.agg <- CalculateSummaryEffect(unpub, summary.measure="SMD", method=method, level=level) 
  } else { unpub.agg <- NA }

  # overall summary effect
  all.agg <- CalculateSummaryEffect(table, summary.measure="SMD", method=method, level=level) 
  
  out <- as.data.frame(rbind(pub.agg,unpub.agg,all.agg))
  return(out)  
} ## END of CalculateSummaryEffectsForOneTable() ##
# library(SAMURAI)
# data(Hpylori)
# effect.sizes.list <- CompileEffectsToBeAssigned(
#   table=Hpylori, binary=TRUE, higher.is.better=TRUE, vpos=NA, pos=NA, neg=NA, vneg=NA,
#   mean.sd=FALSE, binary.measure="RR", continuous.measure="SMD", summary.measure="SMD",
#   method="DL", level=95)
# hp1 <- ImputeBinaryEvents(Hpylori,effect.sizes.list,seed=NA, sims=1)
# hp2 <-  PrepareTableWithBinaryData(table=hp, higher.is.better=F)
# CalculateSummaryEffectsForOneTable(hp2)


CalculateTauSquared <- function(table, measure="RR", method="DL", level=95){
  # Given a table, calculate tau^2 for each of the following subsets:
  #   published; unpublished; all (published & unpublished)
  #
  # Args:
  #
  # Returns: 

  ## for testing
  #   table=table1; level=95; measure="RR"; method="DL"
  
  # published
  pub <- table[which(table$outlook == "published"), ] 
  if(nrow(pub) > 0){
    pub.agg <-  CalculateSummaryEffect(table=pub, level=level,
                  summary.measure="SMD", method=method) 
    pub.tau2 <- pub.agg$tau2
  } else { 
    pub.agg <- NA
  }

  # unpublished
  unpub <- table[which(table$outlook != "published"), ]
  if(nrow(unpub) > 0){
    unpub.agg <- CalculateSummaryEffect(unpub, level=level,
                  summary.measure="SMD", method=method)  
    unpub.tau2 <- unpub.agg$tau2
  } else { 
    unpub.agg <- NA 
  }
  
  # all
  all.agg <- CalculateSummaryEffect(table=table, level=level,
              summary.measure="SMD", method=method)
  all.tau2 <- all.agg$tau2
  
  return(list(pub=pub.tau2, unpub=unpub.tau2, all=all.tau2))  
} ## Definition of CalculateTauSquared() : END ##
# library(SAMURAI)
# data(Hpylori)
# effect.sizes.list <- CompileEffectsToBeAssigned(table=Hpylori, binary=TRUE, higher.is.better=TRUE)
# hp1 <- ImputeBinaryEvents(Hpylori, effect.sizes.list, seed=NA, sims=1)
# hp2 <-  PrepareTableWithBinaryData(table=hp1, higher.is.better=F)
# hp.summaries <- CalculateSummaryEffectsForOneTable(hp2)
# ( hp.tau2 <- CalculateTauSquared(hp2) )


#####################
## Tables: Level 2 ##
#####################


PrepareTableOfSummaries <- function(table, 
  binary=NA, mean.sd=FALSE,
  higher.is.better=NA, 
  vpos=NA,pos=NA,neg=NA,vneg=NA,
  level=95, 
  binary.measure="RR", continuous.measure="SMD",  
  summary.measure="SMD", method="DL",
  seed=NA, noise=0.01, sims=1,
  
  summarize.published.only=FALSE,
  summarize.unpublished.only=FALSE){
  # For each outlook, assign it to all unpublished studies and calculate the summary effect.
  # Then compile a list of summary effects. 
  # 
  # Args:
  #
  # Returns:  A table of summary effects calculated after all unpublished studies are all
  #           assigned the same outlook. 

  tab.summ <- NULL
  
  if(summarize.published.only == TRUE){
      if(binary == TRUE){
          table <- ExtractPublishedStudies(table)
          table <- ConvertBinaryToEffectSize(table, measure=binary.measure)      
      } else {  # continuous
          table <- ExtractPublishedStudies(table)
          if(mean.sd==TRUE){
            table <- ConvertMeanSDToSMD(table,measure=continuous.measure)
          } else{
            table$yi <- table$smd
            table$vi <- table$smd.v
          }
      }

      tab.summ <- as.data.frame(CalculateSummaryEffect(
                  table, summary.measure=summary.measure,method=method, level=level
                               )                        )   
      return(tab.summ)
  }
  
  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  n <- length(outlooks)
  
  if(binary == TRUE){
      for(i in 1:n){
        ## testing
        #   table=Hpylori; binary=T; binary.measure="RR"; sims=1
        #   continuous.measure="SMD"; mean.sd=TRUE; noise=0.01
        #   summary.measure="SMD"
        #   higher.is.better=F; vpos=NA;pos=NA;neg=NA;vneg=NA
        #   level=95; method="DL"; summarize.unpublished.only=TRUE
        #   i <- 10 
        
        tab2 <- PrepareTableWithBinaryData(table,
                  higher.is.better=higher.is.better, 
                  rustlook=outlooks[i], 
                  vpos=vpos, pos=pos, neg=neg, vneg=vneg, 
                  level=level,
                  binary.measure=binary.measure, 
                  summary.measure=summary.measure, method=method,
                  sims=sims, seed=seed)
        
        if(summarize.unpublished.only==TRUE){ 
          tab2 <- ExtractUnpublishedStudies(tab2) 
        }
        
        tab.summ <- rbind(tab.summ,
                      as.data.frame(
                        CalculateSummaryEffect(tab2, level=level,
                          summary.measure=summary.measure, method=method)))                      
      } ## END of for{} loop
  } else{ ## when binary==FALSE
      for(i in 1:n){
        ## testing
        #   table=greentea; binary=F; binary.measure="RR"; sims=1
        #   continuous.measure="SMD"; mean.sd=TRUE; noise=0.01
        #   summary.measure="SMD"
        #   higher.is.better=F; vpos=NA;pos=NA;neg=NA;vneg=NA
        #   level=95; method="DL"; seed=NA; summarize.unpublished.only=T
        #   i <- 7
          
        tab2 <- PrepareTableWithContinuousData(table, mean.sd=mean.sd,
                  higher.is.better=higher.is.better,
                  rustlook=outlooks[i],
                  vpos=vpos, pos=pos, neg=neg, vneg=vneg, 
                  level=level,
                  continuous.measure=continuous.measure, 
                  summary.measure=summary.measure, method=method,
                  seed=seed, noise=noise)

        if(summarize.unpublished.only==TRUE){
          tab2 <- ExtractUnpublishedStudies(tab2)
        }
        
        tab.summ <- rbind(tab.summ,
                      as.data.frame(
                        CalculateSummaryEffect(tab2, level=level,
                          summary.measure=summary.measure,method=method)))         
      } ## END of for() loop
  }
  
  tab.summ <- cbind(outlooks,tab.summ)
  return(tab.summ)
}
# library(SAMURAI)

# data(Hpylori)
# PrepareTableOfSummaries(Hpylori,binary=T,higher.is.better=F, seed=106)
# ( hp.agg.all   <- PrepareTableOfSummaries(Hpylori,binary=T,higher.is.better=F, seed=55, noise=0.01) )
# ( hp.agg.pub   <- PrepareTableOfSummaries(Hpylori,binary=T,higher.is.better=F, seed=55, noise=0.01,
#                                        summarize.published.only=TRUE) )
# ( hp.agg.unpub   <- PrepareTableOfSummaries(Hpylori,binary=T,higher.is.better=F, seed=55, noise=0.01,
#                                          summarize.unpublished.only=TRUE) )

  
# data(greentea)
# ( gt.agg.all   <- PrepareTableOfSummaries(greentea,
#                     binary=F, mean.sd=T,
#                     higher.is.better=F, seed=55, noise=0.01) )
# ( gt.agg.pub   <- PrepareTableOfSummaries(greentea,binary=F,mean.sd=T,higher.is.better=F, seed=55, noise=0.01,
#                                        summarize.published.only=T) )
# ( gt.agg.unpub <- PrepareTableOfSummaries(greentea,binary=F,mean.sd=T,higher.is.better=F, seed=55, noise=0.01,
#                                        summarize.unpublished.only=T) )


###########################
## Forest Plots: Level 0 ##
###########################


GraphBinaryForestPlot <- function(table, 
  higher.is.better=FALSE, 
  level=95, 
  title=NA, scale=1, digits=3, ...){
  # Given a table, produce a forest plot, which includes summary effects. 
  #
  # Args:
  #   table: A data set that includes the following columns:
  #     study:
  #     year: 
  #     outlook:
  #     expt.n, ctrl.n: sample sizes in the two treatment arms
  #     expt.events, ctrl.events: # events in the two treatment arms
  #     yi: effect size estimate (ex. log risk ratio)
  #     vi: estimated variance of effect size estimate
  #   scale: The relative font size. 
  # 
  # Returns: 

  # Dependencies:
  #   Calls: forest() in metafor package, addpoly() in metafor package  
    
  # for testing
  #   table=table5; rr.vpos=rr[1]; rr.pos=rr[2]; rr.neg=rr[3]; rr.vneg=rr[4]; rr.cur=rr[5]; higher.is.better=FALSE; title=NA
  
  # adjust font sizes 
  scale <-  scale * 0.6
  scale2 <- scale * 1.2 * 0.8
  
  # count number of studies - needed to format forest plot 
  num.studies <- nrow(table)
  
  # set limits of plot
  ymin <- -5
  ymax <- num.studies + 3
  xmin <- -16
  xmax <- 8
    
  metafor::forest(table$yi, table$vi, 
                  atransf = exp,                      # to go from logrr to rr
                  ylim = c(ymin,ymax),       # extra rows needed for labels
                  at = log(c(0.05, 0.25, 1, 4, 20)),  # show axis for RR (log scale)
                  xlim = c(xmin, xmax),                   # horizontal dist relative to the vertical line at rr=1
                  slab = paste(table$study, table$year, table$outlook, sep = ", "),  # print author/year
                  ilab = cbind(table$expt.events, table$expt.n, table$ctrl.events, table$ctrl.n),  # print columns with count data
                  ilab.xpos = c(-9.5, -8, -6, -4.5),  # position columns with count data
                  cex = scale,                        # enlarge/reduce font
                  main = title
  )
  # vertical abline at rr=1
  abline(h=0)  
  # add column labels
  text( c(-9.5,-8,-6,-4.5), y=num.studies+2, rep(c("Event", "Total"),2), cex=scale2 )
  text( x=c(-8.75,-5.25), y=num.studies+3, labels=c("Intervention", "Control"), cex=scale2 )
  text( x=xmin, y=num.studies+2, labels="Study", pos=4 , cex=scale2 )
  text( x=xmax, y=num.studies+2, labels="Relative Risk [95% CI]", pos=2 , cex=scale2 )
  
  if(higher.is.better==TRUE){
    text( x=xmax, y=ymin, labels="(Event is GOOD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=4, cex=scale )
  }
  if(higher.is.better==FALSE){
    text( x=xmax, y=ymin, labels="(Event is BAD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=4, cex=scale )
  }
  text ( x=xmin, y=ymin, labels="All effects are estimated with random effects models", pos=4, cex=scale*0.8 )
  
  aggregates <- CalculateSummaryEffectsForOneTable(table, level=level)
  aggregates <- aggregates[1:3,]
  
  ## generate labels; include tau-squared
  agg.tau2 <- CalculateTauSquared(table)
  l.pub     <- paste("Published  ( tau^2 =",round(agg.tau2$pub,digits),")")
  l.unpub   <- paste("Unpublished with specified outlooks ( tau^2 =",round(agg.tau2$unpub,digits),")")
  l.all     <- paste("Published & Unpublished ( tau^2 =",round(agg.tau2$all,digits),")")
  # l.pub     <- paste("Published")
  # l.unpub   <- paste("Unpublished")
  # l.all     <- paste("Published & Unpublished")
  
  agglabels <- c(l.pub, l.unpub, l.all)
  
  addpoly(as.numeric(aggregates$m), sei=as.numeric(aggregates$m.se), atransf=exp, mlab=agglabels, cex=scale)
} ## Definition GraphBinaryForestPlot() : END ##
# library(SAMURAI)
# library(metafor)
# data(Hpylori)
# ( hp <- PrepareTableWithBinaryData(table=Hpylori, 
#           higher.is.better=F, rustlook="negative") )
# GraphBinaryForestPlot(hp)


GraphContinuousForestPlot <- function(table, 
  higher.is.better=TRUE, 
  level=95,
  title=NA, scale=1, digits=3, ...){
  ## graph individual effects and confidence intervals
  
  # for testing
  #   table=table4b
  #   title=NA
  #   higher.is.better=TRUE
  #   level=95
  
  # adjust font sizes 
  scale <-  scale * 0.6
  scale2 <- scale * 1.2
  
  # count number of studies - needed to format forest plot 
  num.studies <- nrow(table)
  
  # set limits of plot
  ymin <- -5
  ymax <- num.studies + 3
  xmin <- -13
  xmax <- 6
  
  # get SMD and its variance from the table
  table$smd <- table$yi
  table$smd.v <- table$vi
  
  table$smd.round <- sprintf("%.3f", round(table$smd, digits) )
  table$smd.v.round <- sprintf("%.3f", round(table$smd.v, digits) )
  
  # make title
#   main.default <- "Forest Plot"
#   subtitle <- ""
#   #   ifelse(higher.is.better==TRUE, 
#   #          subtitle <- " : Event is GOOD", 
#   #          subtitle <- " : Event is BAD")  
#   if(is.na(title) == T){
#     title <- paste(main.default, subtitle, sep="")
#   }
  
  metafor::forest(x=table$smd, 
                  vi=table$smd.v, 
                  ylim = c(ymin,ymax),       # extra rows needed for labels
                  at = c(-1,-0.5, 0, 0.5,1),  # show axis ticks 
                  xlim = c(xmin, xmax),                    # horizontal dist relative to the vertical line at rr=1
                  slab = paste(table$study, table$year, table$outlook, sep = ", "),  # print author/year
                  ilab = cbind(table$expt.n, table$ctrl.n, table$smd.round, table$smd.v.round),  # print columns with count data
                  ilab.xpos = c(-7.5,-6.5,-5.25,-4),  # position columns with count data
                  cex = scale2,                        # enlarge/reduce font
                  main = title
  )
  # vertical abline at smd=0
  abline(h=0)  
  # add column labels
  text( c(-7.5,-6.5,-5.25,-4), y=num.studies+2, c("Expt", "Ctrl", "SMD", "Variance"), cex=scale )
  text( x=c(-7,-4.625), y=num.studies+3, labels=c("Sample size", "SMD"), cex=scale2 )
  text( x=xmin, y=num.studies+2, labels="Study", pos=4 , cex=scale2 )
  text( x=xmax, y=num.studies+2, labels="SMD [95% CI]", pos=2 , cex=scale2 )
  
  if(higher.is.better==TRUE){
    # text( x=xmax, y=ymin, labels="(Event is GOOD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=4, cex=scale )
  }
  if(higher.is.better==FALSE){
    # text( x=xmax, y=ymin, labels="(Event is BAD)", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Intervention", pos=2, cex=scale )
    text( x=0, y=ymin, labels="Favors Control", pos=4, cex=scale )
  }
  text ( x=xmin, y=ymin, labels="All effects are estimated with random effects models", pos=4, cex=scale*0.8 )
  
  aggregates <- CalculateSummaryEffectsForOneTable(table, level=level)
  aggregates <- aggregates[1:3,]
  
  ## generate labels; include tau-squared
  agg.tau2 <- CalculateTauSquared(table)
  l.pub     <- paste("Published  ( tau^2 =",round(agg.tau2$pub,3),")")
  l.unpub   <- paste("Unpublished with specified outlooks ( tau^2 =",round(agg.tau2$unpub,3),")")
  l.all     <- paste("Published & Unpublished ( tau^2 =",round(agg.tau2$all,3),")")
  agglabels <- c(l.pub, l.unpub,l.all)
  
  addpoly(as.numeric(aggregates$m), sei=as.numeric(aggregates$m.se), mlab=agglabels, cex=scale2)
} ## Definition of GraphContinuousForestPlot() : END ##
# library(SAMURAI)
# data(greentea)
# gt <- PrepareTableWithContinuousData(table=greentea, 
#          higher.is.better=F, rustlook="negative")
# GraphContinuousForestPlot(gt)


###########################
## Forest Plots: Level 1 ##
###########################


CalculateAndGraphForestPlot <- function(table,
  binary=TRUE, mean.sd=FALSE,
  higher.is.better=TRUE, 
  rustlook=NA, vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95,
  binary.measure="RR", continuous.measure="SMD", 
  summary.measure="SMD", method="DL",
  seed=NA, noise=0.01, sims=1, 
  title=NA, digits=3, ...){
  # Given a table with binary event data, 
  # this function imputes data, calculates summary effect sizes, 
  # then graphs a forest plot.
  #
  # Args:
  #   level: confidence level = 1-alpha
  #   vpos: "very positive" outlook
  #
  # Returns: a forest plot for a table of pub & unpub studies with binary outcomes.
  #
  # Notes:
  #   Unlike the function GraphBinaryForestPlot(), this function allows
  #   the user to tweak the default effect sizes to be assigned. 
  
  # Dependencies:
  #   Callers: forestsens()
  
  if(binary == TRUE){
    table1 <- PrepareTableWithBinaryData(table,
                higher.is.better=higher.is.better,
                rustlook=rustlook,
                vpos=vpos, pos=pos, neg=neg, vneg=vneg,
                level=level, binary.measure=binary.measure, 
                summary.measure=summary.measure, method=method, 
                seed=seed, sims=sims)

    GraphBinaryForestPlot(table=table1, 
      level=level, higher.is.better=higher.is.better, 
      digits=digits, title=title, ...)     

  } else {
    table1 <- PrepareTableWithContinuousData(table,
                mean.sd=mean.sd, 
                higher.is.better=higher.is.better, 
                vpos=vpos, pos=pos, neg=neg, vneg=vneg,
                rustlook=rustlook, 
                level=level, continuous.measure=continuous.measure, 
                summary.measure=summary.measure, method=method, 
                seed=seed, noise=noise)

    GraphContinuousForestPlot(table=table1, 
      level=level, higher.is.better=higher.is.better, 
      digits=digits, title=title, ...)     
  }
} ## END of CalculateAndGraphForestPlot() ##
# library(SAMURAI)
# data(Hpylori)
# data(greentea)
# CalculateAndGraphForestPlot(Hpylori, binary=TRUE, rustlook="negative", neg=1.8)
# CalculateAndGraphForestPlot(greentea, binary=FALSE, mean.sd=TRUE, rustlook="negative", neg=1.8)


###########################
## Forest Plots: Level 2 ##
###########################


GraphForestPlotForEveryOutlook <- function(table, 
  binary=TRUE, mean.sd=TRUE,
  higher.is.better=FALSE,
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95,  
  sims=1,
  ...){  
  # Summons a forest plot for each one of 10 outlooks.
  # At each turn, all unpublished studies are assigned the same outlook,
  # and a forest plot is generated complete with summary effect sizes
  # for pubs, unpubs, and pubs with unpubs. 
  # 
  # Args:
  #   table: 

  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  n <- length(outlooks)

  if(binary == TRUE){
    for (i in 1:n){
      outlook <- outlooks[i]
      CalculateAndGraphForestPlot(table,
        binary=TRUE,
        level=level,
        higher.is.better=higher.is.better,
        vpos=vpos, 
        pos=pos, 
        neg=neg, 
        vneg=vneg,
        sims=sims,
        rustlook=outlook,
        ...)
    }      
  } else {
    for (i in 1:n){
      outlook <- outlooks[i]
      CalculateAndGraphForestPlot(table,
        binary=FALSE,
        mean.sd=mean.sd,
        level=level,
        higher.is.better=higher.is.better,                     
        rustlook=outlook,
        ...)
    }  
  }
}
# library(SAMURAI)
# data(Hpylori)
# data(greentea)
# GraphForestPlotForEveryOutlook(Hpylori, binary=T, higher.is.better=F, vpos=0.1)
# GraphForestPlotForEveryOutlook(greentea, binary=FALSE, higher.is.better=F)


##################################
## Wrapper (End User) Functions ##
##################################


forestsens <- function(table, 
  # data set format
  binary = TRUE, # outcomes = c("binary","continuous"),
  mean.sd = FALSE, 

  # data set qualities
  higher.is.better = FALSE,  

  # user specified outlooks and effect sizes
  outlook = NA, all.outlooks = FALSE,
  rr.vpos = NA, rr.pos = NA, rr.neg = NA, rr.vneg = NA,
  smd.vpos = NA, smd.pos = NA, smd.neg = NA, smd.vneg = NA,

  # statistical settings
  level = 95, 
  binary.measure = "RR",
  continuous.measure="SMD", 
  summary.measure="SMD",
  method = "DL",

  # random elements
  random.number.seed = NA, 
  sims = 10, 
  smd.noise=0.01,
  
  # output formatting
  plot.title = "", 
  scale = 1,
  digits = 3                      
){
  # Set random.number.seed
  # (Stating it first relieves us of the need to specify it within each function call.)
  if(is.na(random.number.seed) != T) {set.seed(random.number.seed)} 
  
  if(binary == TRUE){  # for binary outcomes

    if(all.outlooks == TRUE){  # for all outlooks
      GraphForestPlotForEveryOutlook(table, binary=TRUE, 
        higher.is.better=higher.is.better, 
        vpos=rr.vpos, pos=rr.pos, neg=rr.neg, vneg=rr.vneg,
        level=level, binary.measure = binary.measure, 
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims,
        title=plot.title, scale=scale, digits=digits)  
      PrepareTableOfSummaries(table, binary=TRUE, 
        higher.is.better=higher.is.better, 
        vpos=rr.vpos, pos=rr.pos, neg=rr.neg, vneg=rr.vneg,
        level=level, binary.measure = binary.measure, 
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims)
    } else {  # for one outlook
      CalculateAndGraphForestPlot(table, binary=TRUE, 
        higher.is.better=higher.is.better, 
        rustlook=outlook,
        vpos=rr.vpos, pos=rr.pos, neg=rr.neg, vneg=rr.vneg,
        level=level, binary.measure = binary.measure, 
        summary.measure=summary.measure, method=method,
        seed=random.number.seed,  sims=sims,
        title=plot.title, scale=scale, digits=digits)  
    }

  } else { # for continuous outcomes

    if(all.outlooks == TRUE){  # for all outlooks
      GraphForestPlotForEveryOutlook(table, binary=FALSE, mean.sd=mean.sd,  
        higher.is.better=higher.is.better, 
        vpos=smd.vpos, pos=smd.pos, neg=smd.neg, vneg=smd.vneg,
        level=level, continuous.measure=continuous.measure,
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims, noise=smd.noise,
        title=plot.title, scale=scale, digits=digits)  
      PrepareTableOfSummaries(table, binary=FALSE, mean.sd=mean.sd, 
        higher.is.better=higher.is.better, 
        vpos=smd.vpos, pos=smd.pos, neg=smd.neg, vneg=smd.vneg,
        level=level, continuous.measure=continuous.measure,
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims, noise=smd.noise)
    } else {  # for one outlook
      CalculateAndGraphForestPlot(table, binary=FALSE, mean.sd=mean.sd,  
        higher.is.better=higher.is.better, 
        rustlook=outlook,
        vpos=smd.vpos, pos=smd.pos, neg=smd.neg, vneg=smd.vneg,
        level=level, continuous.measure=continuous.measure,
        summary.measure=summary.measure, method=method,
        seed=random.number.seed, sims=sims, noise=smd.noise,
        title=plot.title, scale=scale, digits=digits)  
    }

  }
}
# library(metafor)
# load("Hpylori.rda")
# load("greentea.rda")
# forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE)
# forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, plot.title="Test")
# forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, random.number.seed=52)
# forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, random.number.seed=53)
# forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, outlook="negative")
# forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, all.outlooks=TRUE)
# load("greentea.rda")
# forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE)
# forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE,
#   outlook="negative")
# forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE,
#   outlook="negative", smd.noise=0.3)


funnelplot <- function(table,
  binary=TRUE, mean.sd=TRUE,
  higher.is.better=NA,
  outlook=NA,
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95, 
  binary.measure="RR", continuous.measure="SMD",
  summary.measure="SMD", method="DL",
  random.number.seed=NA, sims=1, smd.noise=0.01,
  title="", pch.pub=19, pch.unpub=0){
  
  # to avoid R CMD CHECK NOTE: "no visible binding for global variable"
  yi <- vi <- NULL

  # Set random.number.seed
  # (Stating it first relieves us of the need to specify it within each function call.)
  if(is.na(random.number.seed) != T) {set.seed(random.number.seed)} 

  if(binary == TRUE){
    table <- PrepareTableWithBinaryData(table=table,  
              higher.is.better=higher.is.better,
              rustlook=outlook,
              vpos=vpos, pos=pos, neg=neg, vneg=vneg,
              level=level, 
              binary.measure=binary.measure, 
              summary.measure=summary.measure, 
              method=method, 
              seed=random.number.seed, sims=sims)
  } else {
    table <- PrepareTableWithContinuousData(table=table, 
              mean.sd=mean.sd,
              higher.is.better=higher.is.better,
              rustlook=outlook,
              vpos=vpos, pos=pos, neg=neg, vneg=vneg,
              level=level, 
              continuous.measure=continuous.measure,
              summary.measure=summary.measure, 
              method=method, 
              seed=random.number.seed, noise=smd.noise)
  }

  table$pch <- ifelse(table$outlook=="published", pch.pub, pch.unpub)
  
  if(binary == TRUE){
    res <- rma(yi, vi, data=table, measure=binary.measure, method=method)
    funnel(res, main=title, pch=table$pch, atransf=exp)  
  } else {
    res <- rma(yi, vi, data=table, measure=continuous.measure, method=method)
    funnel(res, main=title, pch=table$pch)
  }  
  
  abline(v=0)
}
# setwd("/Users/nkim41/Dropbox/2013_samurai/SAMURAI_v1.2/01_pre-skeleton/data/rda")
# load("Hpylori.rda")
# load("greentea.rda")
# library(SAMURAI)
# data(Hpylori)
# funnelplot(Hpylori, binary=TRUE, higher.is.better=FALSE, outlook="very negative")
# data(greentea)
# funnelplot(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE)

