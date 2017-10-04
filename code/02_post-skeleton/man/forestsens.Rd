\name{forestsens}
\alias{forestsens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forest Plot for Sensitivity Analysis}
\description{This function imputes missing effect sizes for unpublished studies and creates a forest plot. A set of forest plots can be generated for multiple imputations.}
\usage{
forestsens(table, 
  binary = TRUE, mean.sd = FALSE, 
  higher.is.better = FALSE, 
  outlook = NA, all.outlooks = FALSE, 
  rr.vpos = NA, rr.pos = NA, rr.neg = NA, rr.vneg = NA, 
  smd.vpos = NA, smd.pos = NA, smd.neg = NA, smd.vneg = NA, 
  level = 95, 
  binary.measure = "RR", continuous.measure="SMD", 
  summary.measure="SMD", method = "DL",
  random.number.seed = NA, sims = 10, smd.noise = 0.01,  
  plot.title = "", scale = 1, digits = 3)
}
\arguments{

  \item{table}{The name of the table containing the meta-analysis data.}

  \item{binary}{\code{TRUE} if the outcomes are binary events; FALSE if the outcome data is continuous.}
  \item{mean.sd}{\code{TRUE} if the data set includes the mean and standard deviation of the both the control and experimental arms of studies with continuous outcomes; \code{FALSE} otherwise.}

  \item{higher.is.better}{\code{TRUE} if higher counts of binary events or higher continuous outcomes are desired; \code{FALSE} otherwise. For continuous outcomes, set as FALSE if a lower outcome (eg. a more negative number) is desired.}

  \item{outlook}{If you want all unpublished studies to be assigned the same outcome, set this parameter to one of the following values: \code{"very positive"}, \code{"positive"}, \code{"current effect"}, \code{"negative"}, \code{"very negative"}, \code{"no effect"}, \code{"very positive CL"}, \code{"positive CL"}, \code{"negative CL"}, \code{"very negative CL"}.}
  \item{all.outlooks}{If \code{TRUE}, then a forest plot will be generated for each possible outlook.}

  \item{rr.vpos}{The user-defined relative risk for binary outcomes in unpublished studies with a \code{"very positive"} outlook.}
  \item{rr.pos}{The user-defined relative risk for binary outcomes in unpublished studies with a "positive" outlook.}
  \item{rr.neg}{The user-defined relative risk for binary outcomes in unpublished studies with a "negative" outlook.}
  \item{rr.vneg}{The user-defined relative risk for binary outcomes in unpublished studies with a "very negative" outlook.}
  \item{smd.vpos}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "very positive" outlook.}
  \item{smd.pos}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "positive" outlook.}
  \item{smd.neg}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "negative" outlook.}
  \item{smd.vneg}{The user-defined standardized mean difference for continuous outcomes in unpublished studies with a "very negative" outlook.}

  \item{level}{The confidence level, as a percent.}
  \item{binary.measure}{The effect size measure used for binary outcomes. "RR" for relative risk; "OR" for odds ratios.}
  \item{continuous.measure}{The effect size measure used for continuous outcomes. "SMD" for standardized mean difference (with the assumption of equal variances).}
  \item{summary.measure}{The measure used for summary effect sizes.}
  \item{method}{The same parameter in the escalc() function of the metafor package. "DL" for the DerSimonian-Laird method.}

  \item{random.number.seed}{Leave as \code{NA} if results are to be randomized each time. Set this value to a integer between 0 and 255 if results are to be consistent (for purposes of testing and comparison).}
  \item{sims}{The number of simulations to run per study when imputing unpublished studies with binary outcomes.}
  \item{smd.noise}{The standard deviation of Gaussian random noise to be added to standardized mean differences when imputing unpublished studies with continuous outcomes.}

  \item{plot.title}{Main title of forest plot.}
  \item{scale}{Changes the scaling of fonts in the forest plot.}
  \item{digits}{The number of significant digits (decimal places) to appear in the table of summary results which appears if \code{all.outlooks=TRUE}.}

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

data(Hpylori)
forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE)
forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, plot.title="Test")
forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, random.number.seed=52)
forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, outlook="negative")
forestsens(Hpylori, binary=TRUE, higher.is.better=FALSE, all.outlooks=TRUE)

data(greentea)
forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE)
forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE,
  outlook="negative")
forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE,
  outlook="negative", smd.noise=0.3)
}

\seealso{\code{\link{Hpylori}}, \code{\link{greentea}}}

\keyword{forest plot}
\keyword{meta-analysis}
\keyword{sensitivity analysis}