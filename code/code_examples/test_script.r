setwd("D:\\Documents\\My Dropbox\\2013_samurai\\00CurrentVer\\fn")
source("fn_main_plots.r")
source("fn_aux_plots.r")
source("fn_main_tables.r")
source("fn_aux_tables.r")




setwd("/Users/tamarahowerton/Dropbox/2012_GRA/00CurrentVer/fn")
  source("fn_aux_tables.r")
  source("fn_aux_plots.r")
  source("fn_main_plots.r")
  source("fn_main_tables.r")
  
  setwd("/Users/tamarahowerton/Dropbox/2012_GRA/00CurrentVer/data/meanSD")
  forestsens.importcsv(filename="BHHR2009p88.csv", german.csv=FALSE, binaryoutcome=FALSE, meanssd=TRUE)
  
  library(SAMURAI)
  data(Fleiss1993)
  setwd("/Users/tamarahowerton/Desktop")
  write.csv(Fleiss1993, file="Fleiss.csv")
  forestsens.importcsv(filename="Fleiss.csv", german.csv=FALSE, binaryoutcome=TRUE)