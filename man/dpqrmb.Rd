\name{dmb}
\alias{dmb}
\alias{pmb}
\alias{qmb}
\alias{rmb}
\title{The MBBEFD Distribution}
\description{
Density, distribution function, quantile function and random generation for the
MBBEFD distribution with parameters \code{g} and \code{b}.
}
\usage{
dmb(x, g, b, log = FALSE)
pmb(q, g, b, lower.tail = TRUE, log.p = FALSE)
qmb(p, g, b, lower.tail = TRUE, log.p = FALSE)
rmb(n, g, b)
}
\arguments{
  \item{x, q}{\strong{numeric}; vector of quantiles.}
  \item{p}{\strong{numeric}; vector of probabilities.}
  \item{n}{\strong{numeric}; number of observations.  If \code{length(n) > 1},
  the length is taken to be the number required.}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
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
##--	or standard data sets, see data().

## The function is currently defined as
function (x, g, b, log = FALSE)
{
    .Call(dmb_c, as.double(x), as.double(g), as.double(b), as.logical(log))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
