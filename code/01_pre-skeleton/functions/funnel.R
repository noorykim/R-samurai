# library(metafor)
# ### load ETS data
# data(dat.hackshaw1998)
# ### fit fixed-effects model
# res <- rma(yi, vi, data=dat.hackshaw1998, measure="OR", method="FE")
# ### draw funnel plots
# names(res)
# length(res$yi)
# dotstyle <- c(rep(19,20),rep(15,17))
# funnel(res, main="Standard Error",pch=dotstyle)

funnelplot <- function(table,
  binary, 
  higher.is.better,
  measure=NA,
  title="", pch.pub=19, pch.unpub=0,
  ...){
#   table <- Hpylori; binary=T; higher.is.better=F; pch.pub=19; pch.unpub=0
  
  if(binary == TRUE){
    table <- PrepareTableWithBinaryData(table=table,  
              higher.is.better=higher.is.better,
              ...)
  } else {
    table <- PrepareTableWithContinuousData(table, 
              higher.is.better=higher.is.better,
              ...)
  }

  table$pch <- ifelse(table$outlook=="published", pch.pub, pch.unpub)
  
  if(is.na(measure) == TRUE){
    measure <- ifelse(binary == TRUE,"RR","SMD")
  }
  res <- rma(yi, vi, data=table, measure, method="DL")
  
  if(binary==TRUE){
    funnel(res, main=title, pch=table$pch, atransf=exp)  
  } else{
    funnel(res, main=title, pch=table$pch)    
  }
  abline(v=0)
}
# library(metafor)
# load("Hpylori.rda")
# load("greentea.rda")
# funnelplot(Hpylori,binary=T,higher.is.better=F, rustlook="very negative")
# funnelplot(greentea,binary=F,higher.is.better=F)
