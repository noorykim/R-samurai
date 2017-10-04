\name{Fleiss1993}
\alias{Fleiss1993}
\docType{data}
\title{Aspirin after Myocardial Infarction}
\description{
A meta-analytic data set that includes 7 published placebo-controlled randomized studies of the effect of aspirin in preventing death after myocardial infarction. The data set also includes 2 (fictional) unpublished studies. 

The defined binary outcome event is death, and is undesired. When using the function \code{forestsens()}, specify the option \code{higher.is.better=FALSE}.
}
\usage{data(Fleiss1993)}
\format{A data frame with 9 observations on the following 8 variables.
\tabular{lll}{ \cr
	\bold{number}   \tab \code{integer}   \tab Study numeric id (optional) \cr
	\bold{study}    \tab \code{character} \tab Name of study or principal investigator \cr
	\bold{year}     \tab \code{integer}   \tab Year (optional) \cr
	\bold{outlook}  \tab \code{factor}    \tab Denotes whether a study is unpublished, and if so, what outlook it has.  \cr
	\bold{ctrl.n}  \tab \code{integer} \tab The sample size of the control arm. \cr  
	\bold{expt.n}  \tab \code{integer} \tab The sample size of the experimental arm. \cr  
	\bold{ctrl.events}  \tab \code{integer} \tab The number of (undesired) events within the control arm. \cr  
	\bold{expt.events}  \tab \code{integer} \tab The number of (undesired) events within the experimental arm. 
	}  
}
\details{
The outlook of a study can be one of the following: \code{published}, \code{very positive}, \code{positive}, 	\code{negative}, \code{very negative}, \code{current effect}, \code{no effect}, \code{very positive CL}, \code{positive CL}, \code{negative CL}, or \code{very negative CL}.

Since the outcome event is undesired, when using the function \code{forestsens()}, specify the option \code{higher.is.better=FALSE}.
}
\source{
Fleiss, JL. (1993) "The statistical basis of meta-analysis." \emph{Stat Methods Med Res.} 2(2):121-45.
}
\references{
Guido Schwartzer. \code{meta} package. 
}
\examples{
library(SAMURAI)
data(Fleiss1993)

forestsens(table=Fleiss1993, binary=TRUE, higher.is.better=FALSE)

# To assign all unpublished studies to each of ten outlooks, one at a time, 
# and then return a table of summary effects, their 95% confidence interval, 
# and tau-squared.
summtab <- forestsens(table=Fleiss1993, binary=TRUE, higher.is.better=FALSE, all.outlooks=TRUE)
summtab
}
\keyword{datasets}
