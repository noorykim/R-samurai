setwd("/Users/tamarahowerton/Dropbox/2013_samurai/SAMURAI v1.01/01 pre-skeleton/data")
greentea <- read.csv("greentea/greentea.csv", header=T)
save(greentea, file="greentea.rda")

greentea_original <- read.csv("greentea/greentea_original.csv", header=T)
save(greentea_original, file="greentea_original.rda")

Hpylori <- read.csv("Hpylori/Hpylori.csv", header=T)
save(Hpylori, file="Hpylori.rda")

Hpylori_original <- read.csv("Hpylori/Hpylori_original.csv", header=T)
Hpylori_original <- Hpylori_original[,1:6]
save(Hpylori_original, file="Hpylori_original.rda")