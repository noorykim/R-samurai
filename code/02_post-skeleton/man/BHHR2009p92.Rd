\name{BHHR2009p92}
\alias{BHHR2009p92}
\docType{data}
\title{Fictional Data Set, with Binary Outcomes}
\description{A fictional meta-analytic data set with 6 published studies and 2 unpublished studies. The binary outcome event is not desired.}
\usage{data(BHHR2009p92)}
\format{A data frame with 8 observations on the following 8 variables.\tabular{lll}{ \cr
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

Since the outcome event is undesired, when using the function \code{forestsens()}, specify the option \code{higher.is.better=FALSE}. Since this is the default setting for \code{forestsens()}, this does not need to be specified explicitly. 
}
\source{
Borenstein, Hedges, Higgins, and Rothstein. Introduction to Meta-analysis. Wiley, 2009, page 92.  
}
\examples{
library(SAMURAI)
data(BHHR2009p92)

forestsens(table=BHHR2009p92, binary=TRUE, higher.is.better=FALSE)

# To assign all unpublished studies to each of ten outlooks, one at a time, 
# and then return a table of summary effects, their 95% confidence interval, 
# and tau-squared.
summtab <- forestsens(table=BHHR2009p92, binary=TRUE, higher.is.better=FALSE, all.outlooks=TRUE)
summtab
}
\keyword{datasets}
