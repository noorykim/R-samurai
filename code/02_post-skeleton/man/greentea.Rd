\name{greentea}
\alias{greentea}
\docType{data}
\title{The effect of green tea on weight loss.}
\description{
Randomized clinical trials of at least 12 weeks duration assessing the effect of green tea consumption on weight loss. 
}
\usage{data(greentea)}
\format{
  A data frame with 14 observations on the following 9 variables.
  \tabular{lll}{ \cr
    \bold{study}    \tab \code{character} \tab Name of study or principal investigator \cr
    \bold{year}     \tab \code{numeric (integer)}   \tab Year (optional) \cr
    \bold{outlook}  \tab \code{factor}    \tab Denotes whether a study is unpublished, and if so, what outlook it has.  \cr
    \bold{ctrl.n}  \tab \code{numeric (integer)} \tab The sample size of the control arm. \cr  
    \bold{expt.n}  \tab \code{numeric (integer)} \tab The sample size of the experimental arm. \cr 
    \bold{ctrl.mean}  \tab \code{numeric} \tab The mean effect within the control arm. \cr  
    \bold{expt.mean}  \tab \code{numeric} \tab The mean effect within the experimental arm. \cr
    \bold{ctrl.sd}  \tab \code{numeric} \tab The standard deviation of the outcome within the control arm. \cr  
    \bold{expt.sd}  \tab \code{numeric} \tab The standard deviation of the outcome within the experimental arm. 
  }  
}
\details{
The outlook of a study can be one of the following: \code{published}, \code{very positive}, \code{positive},  \code{negative}, \code{very negative}, \code{current effect}, \code{no effect}, \code{very positive CL}, \code{positive CL}, \code{negative CL}, or \code{very negative CL}.

In this setting, a more negative change in outcome is desired; specify the option \code{higher.is.better=FALSE} for the function \code{forestsens()}.
}
\source{
Jurgens TM, Whelan AM, Killian L, Doucette S, Kirk S, Foy E. "Green tea for weight loss and weight maintenance in overweight or obese adults." \emph{Cochrane Database of Systematic Reviews} 2012, Issue 12. Art. No.: CD008650. DOI: 10.1002/14651858.CD008650.pub2.

Figure 6. Forest plot of comparison: 1 Primary outcomes, outcome: 1.2Weight loss studies conducted in/outside Japan.
}
\examples{
data(greentea)
greentea

forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE)

# To fix the random number seed to make the results reproducible. 
forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE, random.number.seed=52)

# To modify the outlooks of all unpublished studies to, say, "negative". 
forestsens(greentea, binary=FALSE,mean.sd=TRUE,higher.is.better=FALSE,random.number.seed=52,
  outlook="negative")

# To modify the outlooks of all unpublished studies to, say, "negative", and 
# overruling the default standardized mean difference (SMD) assigned to "negative". 
# (In this case, for a negative outlook we might assign a positive SMD, which corresponds to 
# having weight loss under green tea treatment less than weight loss under control treatment, 
# i.e. the green tea treatment is less effective at achieving weight loss than control treatment.
forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE,random.number.seed=52,
  outlook="negative", smd.neg=0.4)

# To generate a forest plot for each of the ten default outlooks defined by forestsens().
forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE, random.number.seed=52,
  all.outlooks=TRUE)
}
\keyword{datasets}
