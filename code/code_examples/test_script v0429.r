library(SAMURAI)
sessionInfo() ## check version of packages

## binary
data("greentea")
greentea
forestsens(greentea,binaryoutcome=F,meanssd=T,scale=0.8)

