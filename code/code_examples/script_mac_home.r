## Filename: script_mac.r
## Author: Noory Kim
## Last updated: 08 January 2013

#################################################################################
## The following commands only need to be done once per computer. 
install.packages("metafor") ## To install the metafor package
install.packages("SAMURAI") ## To install the SAMURAI package

#################################################################################
## The following commands only need to be done once per R session.
require(SAMURAI) ## To load the SAMURAI package. This will also load the metafor package.

#################################################################################
## To load a dataset in the SAMURAI package: 
data(SMIdrugD2) ## To load the SMIdrugD2 dataset. 
SMIdrugD2       ## To display the SMIdrugD2 dataset.

## Datasets in the SAMURAI package: 
## BHHR2009p92, Fleiss1993, SMIbehavior, SMIdrugD2, SMIdrugD3

#################################################################################
## To change the working directory, where files will be read and/or saved.
getwd() ## To see the current working directory
setwd("/Users/tamarahowerton/Desktop/workdir") ## To change the working directory on a Mac computer. An example. 

#################################################################################
## To load an external dataset saved as a CSV file. 

?read.csv ## To see the help file for the read.csv() and read.csv2() functions. 

setwd("/Users/tamarahowerton/Desktop/workdir/data") ## To change to the directory with the CSV data file. An example.

## To load an external dataset saved as a CSV file (not a German CSV file) as a table in R.
datatable1 <- read.csv(file="datafile1.csv")
datatable1 <- read.csv(file="datafile1.csv", sep=",", dec=".") ## Equivalent to the previous command.

## To load an external dataset saved as a German CSV file 
## (a CSV file on a Windows computer with German language settings). 
datatable2 <- read.csv2(file="datafile2.csv")
datatable2 <- read.csv2(file="datafile1.csv", sep=";", dec=",") ## Equivalent to the previous command.

#################################################################################
## To generate a forest plot. 

## SMIdrugD2 data set:
## -- Outcomes are continuous. (binaryoutcome=FALSE) 
## -- Effects are presented as means and standard deviations. (meanssd=TRUE)
## -- A lower outcome is desired. (event.is.good=FALSE; default setting)
forestsens(SMIdrugD2, binaryoutcome=FALSE, meanssd=TRUE, plot.title="SMI Drug D2")
forestsens(SMIdrugD2, binaryoutcome=FALSE, meanssd=TRUE, plot.title="SMI Drug D2", foralloutlooks=TRUE)

## To generate a forest plot as a PDF file. 
setwd("/Users/JohnDoe/Desktop/workdir/plots") ## To change to the directory where the plot will be saved. An example.
pdf("D2asis.pdf")
forestsens(SMIdrugD2, binaryoutcome=FALSE, meanssd=TRUE, plot.title="SMI Drug D2")
dev.off()