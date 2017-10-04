## Last updated: 8/22/2013

## Outline: 
## 1. Load *.rda data sets and functions. 
## 2. Invoke package skeleton. 
## 3. Flesh out the skeleton (by copying files into it).

## 0. Clear R workspace.

rm(list=ls()) 


## 1. Load *.rda data sets and functions. 

# dir.machine <- "/Users/nkim41/"
dir.machine <- "/Users/tamarahowerton/"
dir.version <- "Dropbox/2013_samurai/SAMURAI_v1.2/"
dir.main <- paste(dir.machine,dir.version,sep="")
setwd(dir.main)
getwd()

dir.data <- paste(dir.main,"01_pre-skeleton/data/rda", sep="")
setwd(dir.data)
load("Hpylori.rda")
load("greentea.rda")
load("Fleiss1993.rda")
load("BHHR2009p92.rda")
# load("Hpylori_original.rda")
# load("greentea_original.rda")

dir.fn <- paste(dir.main,"01_pre-skeleton/functions",sep="")
setwd(dir.fn)
getwd()

library(metafor)
source("bushido.r")


## 2. Invoke skeleton. 

dir.build.platform <- paste(dir.machine,"Desktop/build/",sep="")  # requires desktop to have a 'build' folder
dir.build <- paste(dir.build.platform,"SAMURAI/",sep="")

setwd(dir.build.platform)
getwd()
unlink("SAMURAI", recursive=T) ## delete existing SAMURAI skeleton
package.skeleton("SAMURAI")


## 3. Flesh out the skeleton. 

## Remove the file "Read-and-delete-me"
setwd(dir.build)
getwd()
file.remove("Read-and-delete-me")

## Transfer completed DESCRIPTION and NAMESPACE files into skeleton.
file.copy(from=paste(dir.main,"01_pre-skeleton/NAMESPACE",sep=""),
          to=paste(dir.build,"NAMESPACE",sep=""),overwrite=TRUE)
file.copy(from=paste(dir.main,"01_pre-skeleton/DESCRIPTION",sep=""), 
          to=paste(dir.build,"DESCRIPTION",sep=""),overwrite=TRUE)

## Transfer MAN directory
src.dir <- paste(dir.machine, dir.version,"02_post-skeleton/man/",sep="")
dest.dir <- paste(dir.build,"man/",sep="")
file.names <- dir(src.dir)
unlink(dest.dir, recursive=T) ## delete existing MAN directory
dir.create(dest.dir)  # generate new MAN directory
sapply(file.names, function(x) {file.copy(from=paste(src.dir, x, sep=''), 
                                          to=paste(dest.dir, x, sep=''), 
                                          overwrite = FALSE)})

setwd(paste(dir.build, "data",sep=""))
file.remove("dir.build.platform.rda", "dir.build.rda", "dir.data.rda", "dir.fn.rda")
file.remove("dir.machine.rda", "dir.main.rda", "dir.version.rda")
file.remove("dest.dir.rda", "file.names.rda", "src.dir.rda")