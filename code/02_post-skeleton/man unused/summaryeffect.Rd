\name{summaryeffect}
\alias{summaryeffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summaryeffect(table, measure = "RR", method = "DL", confidencelevel = 95, exp = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{measure}{
%%     ~~Describe \code{measure} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{confidencelevel}{
%%     ~~Describe \code{confidencelevel} here~~
}
  \item{exp}{
%%     ~~Describe \code{exp} here~~
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
function (table, measure = "RR", method = "DL", confidencelevel = 95, 
    exp = TRUE) 
{
    summ <- summarizeeffect(table, method = method, confidencelevel = confidencelevel)
    if (exp == TRUE) {
        effectsize <- summ$exp.m
        ci <- c(summ$exp.m.lcl, summ$exp.m.ucl)
    }
    else {
        effectsize <- summ$m
        ci <- c(summ$m.lcl, summ$m.ucl)
    }
    tau2 <- summ$tau2
    out <- as.list(c(effectsize, ci, tau2))
    names(out) <- c("effect", "lcl", "ucl", "tau2")
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
