\name{tablesens}
\alias{tablesens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tablesens(table, binaryoutcome = TRUE, outcome = c("binary", "continuous"), binarymeasure = "RR", rr.vpos = NA, rr.pos = NA, rr.neg = NA, rr.vneg = NA, meanssd = FALSE, smd.vpos = NA, smd.pos = NA, smd.neg = NA, smd.vneg = NA, event.is.good = FALSE, simsperstudy = 10, confidencelevel = 95, method = "DL", sigdigits = 3, scale = 1, random.number.seed = NA, smd.noise = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{binaryoutcome}{
%%     ~~Describe \code{binaryoutcome} here~~
}
  \item{outcome}{
%%     ~~Describe \code{outcome} here~~
}
  \item{binarymeasure}{
%%     ~~Describe \code{binarymeasure} here~~
}
  \item{rr.vpos}{
%%     ~~Describe \code{rr.vpos} here~~
}
  \item{rr.pos}{
%%     ~~Describe \code{rr.pos} here~~
}
  \item{rr.neg}{
%%     ~~Describe \code{rr.neg} here~~
}
  \item{rr.vneg}{
%%     ~~Describe \code{rr.vneg} here~~
}
  \item{meanssd}{
%%     ~~Describe \code{meanssd} here~~
}
  \item{smd.vpos}{
%%     ~~Describe \code{smd.vpos} here~~
}
  \item{smd.pos}{
%%     ~~Describe \code{smd.pos} here~~
}
  \item{smd.neg}{
%%     ~~Describe \code{smd.neg} here~~
}
  \item{smd.vneg}{
%%     ~~Describe \code{smd.vneg} here~~
}
  \item{event.is.good}{
%%     ~~Describe \code{event.is.good} here~~
}
  \item{simsperstudy}{
%%     ~~Describe \code{simsperstudy} here~~
}
  \item{confidencelevel}{
%%     ~~Describe \code{confidencelevel} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{sigdigits}{
%%     ~~Describe \code{sigdigits} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{random.number.seed}{
%%     ~~Describe \code{random.number.seed} here~~
}
  \item{smd.noise}{
%%     ~~Describe \code{smd.noise} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (table, binaryoutcome = TRUE, outcome = c("binary", 
    "continuous"), binarymeasure = "RR", rr.vpos = NA, rr.pos = NA, 
    rr.neg = NA, rr.vneg = NA, meanssd = FALSE, smd.vpos = NA, 
    smd.pos = NA, smd.neg = NA, smd.vneg = NA, event.is.good = FALSE, 
    simsperstudy = 10, confidencelevel = 95, method = "DL", sigdigits = 3, 
    scale = 1, random.number.seed = NA, smd.noise = 0.01) 
{
    if (is.na(random.number.seed) != T) {
        set.seed(random.number.seed)
    }
    if (binaryoutcome == TRUE) {
        summarize.eachoutlook(table, binaryoutcome = TRUE, confidencelevel = confidencelevel, 
            event.is.good = event.is.good, rr.vpos = rr.vpos, 
            rr.pos = rr.pos, rr.neg = rr.neg, rr.vneg = rr.vneg, 
            simsperstudy = simsperstudy, measure = binarymeasure, 
            method = method, sigdigits = sigdigits, smd.noise = smd.noise)
    }
    else {
        summarize.eachoutlook(table, binaryoutcome = FALSE, meanssd = meanssd, 
            confidencelevel = confidencelevel, event.is.good = event.is.good, 
            smd.vpos = smd.vpos, smd.pos = smd.pos, smd.neg = smd.neg, 
            smd.vneg = smd.vneg, method = method, sigdigits = sigdigits, 
            smd.noise = smd.noise, exp = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
