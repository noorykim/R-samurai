\name{forestsens}
\alias{forestsens}
\title{Forest Plot for Sensitivity Analysis}
\description{This function imputes missing effect sizes for unpublished studies and creates a forest plot. A set of forest plots can be generated for multiple imputations.}
\usage{
forestsens(table, binaryoutcome = TRUE, meanssd = FALSE, event.is.good = FALSE,  
	unpub.oneoutlook = NA, foralloutlooks = FALSE,  
	random.number.seed = NA, simsperstudy = 10, 
	sigdigits = 3, confidencelevel = 95, 
	binarymeasure = "RR",  method = "DL",
	rr.vpos = NA, rr.pos = NA, rr.neg = NA, rr.vneg = NA, 
	smd.vpos = NA, smd.pos = NA, smd.neg = NA, smd.vneg = NA,
	plot.title = "Random Effects Forest Plot", scale = 1)
}
\arguments{
  \item{table}{The name of the table containing the meta-analysis data.}
  \item{meanssd}{\code{TRUE} if the data set includes the mean and standard deviation of the both the control and experimental arms of studies with continuous outcomes; \code{FALSE} otherwise }
  \item{confidencelevel}{The percent confidence level.}
  \item{event.is.good}{\code{TRUE} if the event is desired; \code{FALSE} otherwise. For continuous outcomes, set as TRUE if a higher outcome is desired, and set as FALSE if a lower outcome is desired.}
  \item{random.number.seed}{Leave as \code{NA} if results are to be randomized each time. Set this value to a integer between 0 and 255 if results are to be consistent (for purposes of testing and comparison).}
  \item{simsperstudy}{The number of simulations run for each unpublished study.}
  \item{plot.title}{Forest plot title.}
  \item{scale}{Changes the scaling of elements in the forest plot.}
  \item{unpub.oneoutlook}{If you want all unpublished studies to be assigned the same outcome, set this parameter to one of the following values: \code{"very positive"}, \code{"positive"}, \code{"current effect"}, \code{"negative"}, \code{"very negative"}, \code{"no effect"}, \code{"very positive CL"}, \code{"positive CL"}, \code{"negative CL"}, \code{"very negative CL"}.}
  \item{foralloutlooks}{If \code{TRUE}, then a forest plot will be generated for each possible outlook.}
  \item{binaryoutcome}{\code{TRUE} if the event is binary; FALSE if the event is continuous.}
  \item{rr.vpos}{The user-defined relative risk for binary outcomes in unpublished studies with a \code{"very positive"} outlook.}
  \item{rr.pos}{The user-defined relative risk for binary outcomes in unpublished studies with a "positive" outlook.}
  \item{rr.neg}{The user-defined relative risk for binary outcomes in unpublished studies with a "negative" outlook.}
  \item{rr.vneg}{The user-defined relative risk for binary outcomes in unpublished studies with a "very negative" outlook.}
  \item{smd.vpos}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "very positive" outlook.}
  \item{smd.pos}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "positive" outlook.}
  \item{smd.neg}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "negative" outlook.}
  \item{smd.vneg}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "very negative" outlook.}
  \item{binarymeasure}{"RR" for relative risk; "OR" for odds ratios.}
  \item{method}{The same parameter in the escalc() function of the metafor package. "DL" for the DerSimonian-Laird method. }
  \item{sigdigits}{The number of significant digits (decimal places) to appear in the table of summary results which appears if \code{foralloutlooks=TRUE}.}
}
\details{
For unpublished studies with binary outcomes, random numbers are generated from binomial distributions to impute the number of events in the experimental arms of experimental studies. The parameter of these distributions depends out the outlook of the unpublished study and the rate of events in the control arms of published studies. By default, 10 simulations are run and their average is used to impute the number of events in the experimental arm. 

For unpublished studies with continuous outcomes, a 'very good' approximator mentioned by Borenstein is used to impute the variance of the standardized mean difference. See Borenstein et al, 2009, pages 27-28.
}
\references{
Borenstein M, Hedges LV, Higgins JPT, and Rothstein HR (2009). Introduction to Meta-Analysis. Chichester UK: Wiley. 

Cooper HC, Hedges LV, & Valentine JC, eds. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd ed.). New York: Russell Sage Foundation.

DerSimonian R and Laird N (1986). "Meta-analysis in clinical trials." \emph{Controlled Clinical
Trials} 7:177-188 (1986).

Viechtbauer W (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\author{Noory Kim}
\note{The function employs functions in the \code{metafor} package: \code{escalc()} and \code{forest()}.}
\examples{
library(SAMURAI)
data(Fleiss1993)
write.csv(Fleiss1993, file="Fleiss.csv")
forestsens.importcsv(filename="Fleiss.csv", german.csv=FALSE, binaryoutcome=TRUE)
}
\seealso{
\code{\link{Fleiss1993}}, \code{\link{SMIbehavior}}
}