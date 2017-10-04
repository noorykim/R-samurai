# setwd("/Users/tamarahowerton/Dropbox/2013_samurai/00CurrentVer/fn")
setwd("/home/nk/Dropbox/2013_samurai/00CurrentVer/fn")
source("fn_wrappers.r")
source("fn_aux_tables.r")
source("fn_aux_plots.r")
source("fn_main_plots.r")
source("fn_main_tables.r")


library(metafor)

## binary
setwd("/Users/tamarahowerton/Dropbox/2013_samurai/00CurrentVer/data/binary")
load("BHHR2009p92.rda")
bor   <- BHHR2009p92
bor.pub <- convertbin2effectsize(bor[which(bor$outlook=="published"),], measure="OR")
names(bor.pub)
bor.summ <- summarizeeffect(bor.pub, confidencelevel=95)

load("Fleiss1993.rda")
convertbin2effectsize(Fleiss1993[which(Fleiss1993$outlook=="published"),])
fl.pub <- convertbin2effectsize(Fleiss1993[which(Fleiss1993$outlook=="published"),])
names(fl.pub)
fl.summ <- summarizeeffect(fl.pub, confidencelevel=95)
fl2 <- impute.ctrl.events(Fleiss1993, roundtable=FALSE) 

fl.rr <- assignrr(fl.summ$exp.m, fl.summ$exp.m.lcl, fl.summ$exp.m.ucl, event.is.good=FALSE)

impute.expt.events(fl2,  fl.rr, simsperstudy=1)


## continuous
# setwd("/Users/tamarahowerton/Dropbox/2013_samurai/00CurrentVer/data/meanSD")
setwd("/home/nk/Dropbox/2013_samurai/00CurrentVer/data/meanSD")
load("SMIdrugD2.rda")
sd2.pub <- convertmeans2smd(SMIdrugD2[which(SMIdrugD2$outlook=="published"),])
sd2.summ <- summarizeeffect(sd2.pub, confidencelevel=95)

setwd("/Users/tamarahowerton/Dropbox/2013_samurai/00CurrentVer/data/meanSD")
load("SMIbehavior.rda")
sb <- SMIbehavior
sb <- convertmeans2smd(sb)
sb.pub <- convertmeans2smd(SMIbehavior[which(SMIbehavior$outlook=="published"),])
sb.summ <- summarizeeffect(sb.pub, confidencelevel=95)
sb.smd <- assignsmd(event.is.good=F, pubsmd=sb.summ$exp.m, pubsmd.lcl=sb.summ$exp.m.lcl, pubsmd.ucl=sb.summ$exp.m.ucl) 
impute.smd(sb,sb.smd, smd.noise=0.1) 

setwd("/Users/tamarahowerton/Dropbox/2013_samurai/00CurrentVer/data/meanSD")
load("SMIbehavior.rda")
sb <- SMIbehavior
forestsens(sb, outcome="continuous", meanssd = TRUE, event.is.good = TRUE, smd.noise=0.02,
           random.number.seed=106)
fst <- tablesens(sb,outcome="continuous", meanssd = TRUE, event.is.good = TRUE, smd.noise=0.02,
                 random.number.seed=14)

fst$effect <- as.numeric(levels(fst$effect))[fst$effect] 
fst$tau2 <- as.numeric(levels(fst$tau2))[fst$tau2] 
plot(x=1:10, y=fst$effect, ylim=c(0,1))
points(x=1:10, y=fst$tau2, col="red")
lines(1:5,fst$effect[1:5])
lines(6:10,fst$effect[6:10])
lines(1:5,fst$tau2[1:5])
lines(6:10,fst$tau2[6:10])
