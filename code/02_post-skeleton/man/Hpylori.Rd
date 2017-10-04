\name{Hpylori}
\alias{Hpylori}
\docType{data}
\title{Healing of duodenal ulcers by Helicobacter pylori eradication therapy}
\description{
Randomized clinical trials comparing duodenal ulcer acute healing among (1) patients on ulcer healing drug + Helicobacter pylori eradication therapy vs. (2) patients ulcer healing drug alone. The event counts represent the numbers of patients \emph{not} healed.
}
\usage{data(Hpylori)}
\format{
  A data frame with 33 observations on the following 7 variables.
  \tabular{lll}{ \cr
    \bold{study}    \tab \code{character} \tab Name of study or principal investigator \cr
    \bold{year}     \tab \code{numeric}   \tab Year (optional) \cr
    \bold{outlook}  \tab \code{factor}    \tab Denotes whether a study is unpublished, and if so, what outlook it has.  \cr
    \bold{ctrl.n}  \tab \code{numeric} \tab The sample size of the control arm. \cr  
    \bold{expt.n}  \tab \code{numeric} \tab The sample size of the experimental arm. \cr  
    \bold{ctrl.events}  \tab \code{numeric} \tab The number of (undesired) events within the control arm. \cr  
    \bold{expt.events}  \tab \code{numeric} \tab The number of (undesired) events within the experimental arm. 
  }  
}
\details{
The outlook of a study can be one of the following: \code{published}, \code{very positive}, \code{positive},  \code{negative}, \code{very negative}, \code{current effect}, \code{no effect}, \code{very positive CL}, \code{positive CL}, \code{negative CL}, or \code{very negative CL}.

Since the outcome event is undesired, when using the function \code{forestsens()}, specify the option \code{higher.is.better=FALSE}. 
}
\source{
Ford AC, Delaney B, Forman D, Moayyedi P. "Eradication therapy for peptic ulcer disease in \emph{Helicobacter pylori} positive patients." \emph{Cochrane Database of Systematic Reviews} 2006, Issue 2. Art No.: CD003840. DOI: 10.1002/14651858.CD003840.pub4. 

Figure 3. Forest plot of comparison: 1 duodenal ulcer acute healing hp eradication + ulcer healing drug vs. ulcer healing drug alone, outcome: 1.1 Proportion not healed.
}
\examples{
data(Hpylori)
Hpylori

forestsens(table=Hpylori, binary=TRUE, higher.is.better=FALSE, scale=0.8)

# To fix the random number seed to make the results reproducible.
forestsens(table=Hpylori, binary=TRUE, higher.is.better=FALSE, scale=0.8, 
  random.number.seed=106)

# To modify the outlooks of all unpublished studies to, say, "very negative". 
forestsens(table=Hpylori, binary=TRUE, higher.is.better=FALSE, scale=0.8, 
  random.number.seed=106, outlook="very negative")

# To modify the outlooks of all unpublished studies to, say, "very negative", 
# and overruling the default relative risk assigned to "very negative".
forestsens(table=Hpylori, binary=TRUE, higher.is.better=FALSE, scale=0.8, 
  random.number.seed=106, outlook="very negative", rr.vneg=2.5)

# To generate a forest plot for each of the ten default outlooks 
# defined by forestsens().
forestsens(table=Hpylori, binary=TRUE, higher.is.better=FALSE, scale=0.8, 
  random.number.seed=106, all.outlooks=TRUE)
}
\keyword{datasets}
