\name{funnelplot}
\alias{funnelplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Funnel Plot}
\description{
This function (1) imputes data for a meta-analytic data set with unpublished studies, then (2) generates a funnel plot.}
\usage{
funnelplot(table, 
  binary=TRUE, mean.sd=TRUE,
  higher.is.better=NA,
  outlook=NA,
  vpos=NA, pos=NA, neg=NA, vneg=NA,
  level=95, 
  binary.measure="RR", continuous.measure="SMD",
  summary.measure="SMD", method="DL", 
  random.number.seed=NA, sims=1, smd.noise=0.01,
  title="", pch.pub=19, pch.unpub=0)
}
\arguments{
  \item{table}{The name of the table containing the meta-analysis data.}

  \item{binary}{\code{TRUE} if the outcomes are binary events; FALSE if the outcome data is continuous.}
  \item{mean.sd}{\code{TRUE} if the data set includes the mean and standard deviation of the both the control and experimental arms of studies with continuous outcomes; \code{FALSE} otherwise.}

  \item{higher.is.better}{\code{TRUE} if higher counts of binary events or higher continuous outcomes are desired; \code{FALSE} otherwise. For continuous outcomes, set as FALSE if a lower outcome (eg. a more negative number) is desired.}

  \item{outlook}{If you want all unpublished studies to be assigned the same outcome, set this parameter to one of the following values: \code{"very positive"}, \code{"positive"}, \code{"current effect"}, \code{"negative"}, \code{"very negative"}, \code{"no effect"}, \code{"very positive CL"}, \code{"positive CL"}, \code{"negative CL"}, \code{"very negative CL"}.}

  \item{vpos}{The user-defined effect size for  unpublished studies with a \code{"very positive"} outlook.}
  \item{pos}{The user-defined effect size for  unpublished studies with a \code{"positive"} outlook.}
  \item{neg}{The user-defined effect size for  unpublished studies with a \code{"negative"} outlook.}
  \item{vneg}{The user-defined effect size for  unpublished studies with a \code{"very negative"} outlook.}

  \item{level}{The confidence level, as a percent.}
  \item{binary.measure}{The effect size measure used for binary outcomes. "RR" for relative risk; "OR" for odds ratios.}
  \item{continuous.measure}{The effect size measure used for continuous outcomes. "SMD" for standardized mean difference (with the assumption of equal variances).}
  \item{summary.measure}{The measure used for summary effect sizes.}
  \item{method}{The same parameter in the escalc() function of the metafor package. "DL" for the DerSimonian-Laird method.}

  \item{random.number.seed}{Leave as \code{NA} if results are to be randomized each time. Set this value to a integer between 0 and 255 if results are to be consistent (for purposes of testing and comparison).}
  \item{sims}{The number of simulations to run per study when imputing unpublished studies with binary outcomes.}
  \item{smd.noise}{The standard deviation of Gaussian random noise to be added to standardized mean differences when imputing unpublished studies with continuous outcomes.}

  \item{title}{The title of the funnel plot.}
  \item{pch.pub}{The symbol used to denote a published study.}
  \item{pch.unpub}{The symbol used to denote an unpublished study.}
}

\author{Noory Kim}

\note{The function employs functions in the \code{metafor} package: \code{escalc()} and \code{forest()}.}

\seealso{\code{\link{forestsens}}}

\examples{
library(SAMURAI)

data(Hpylori)
funnelplot(Hpylori, binary=TRUE, higher.is.better=FALSE, 
  outlook="very negative")

data(greentea)
funnelplot(greentea, binary=FALSE, higher.is.better=FALSE)
}

\keyword{funnel plot}
\keyword{meta-analysis}
\keyword{sensitivity analysis}

