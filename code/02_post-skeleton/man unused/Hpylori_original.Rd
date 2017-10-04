\name{Hpylori_original}
\alias{Hpylori_original}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Hpylori_original)}
\format{
  A data frame with 33 observations on the following 10 variables.
  \describe{
    \item{\code{study}}{a factor with levels \code{Asaka} \code{Avsar} \code{Bardhan} \code{Bayrerdorffer} \code{Bianchi Porro} \code{Carpintero} \code{Figueroa} \code{Graham} \code{Harford} \code{HentscheI} \code{Hosking} \code{Kato} \code{Katoh} \code{Kepecki} \code{Lin} \code{Logan} \code{Mantzaris} \code{Mones} \code{O'Morain} \code{Parente} \code{Pinero} \code{Porro} \code{Pounder} \code{Rauws} \code{Schwartz} \code{Shirotani} \code{Sobhani} \code{Spinzi} \code{van Zanten} \code{Wang} \code{Wong}}
    \item{\code{year}}{a numeric vector}
    \item{\code{expt.events}}{a numeric vector}
    \item{\code{expt.n}}{a numeric vector}
    \item{\code{ctrl.events}}{a numeric vector}
    \item{\code{ctrl.n}}{a numeric vector}
    \item{\code{X}}{a logical vector}
    \item{\code{X.1}}{a logical vector}
    \item{\code{X.2}}{a logical vector}
    \item{\code{X.3}}{a logical vector}
  }
}
\details{
Figure 3. Forest plot of comparison: 1 duodenal ulcer acute healing hp eradication + ulcer healing drug vs. ulcer healing drug alone, outcome: 1.1 Proportion not healed.
}
\source{
Ford AC, Delaney B, Forman D, Moayyedi P. Eradication therapy for peptic ulcer disease in Helicobacter pylori positive patients. Cochrane Database of Systematic Reviews 2006, Issue 2. Art No.: CD003840. DOI: 10.1002/14651858.CD003840.pub4. 
}
\examples{
data(Hpylori_original)
## maybe str(Hpylori_original) ; plot(Hpylori_original) ...
}
\keyword{datasets}
