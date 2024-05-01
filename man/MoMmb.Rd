\name{MoMmb}
\alias{MoMmb}
\title{
Method of Moments Parameter Estimation for the MBBEFD distribution.
}
\description{
Attempts to find the best \code{g} and \code{b} parameters which are consistent
with the first and second moments of the supplied data.
}
\usage{
MoMmb(x, maxit = 100L, tol = .Machine$double.eps ^ 0.5, na.rm = TRUE)
}
\arguments{
  \item{x}{\strong{numeric}; vector of observations between 0 and 1.}
  \item{maxit}{\strong{integer}; maximum number of iterations.}
  \item{tol}{\strong{numeric}; tolerance. If too tight, algorithm may fail.
  Defaults to the square root of \code{.Machine$double.eps} or roughly
  \eqn{1.49\times 10^{-8}}{1.49e-8}.}
  \item{na.rm}{\strong{logical}; if TRUE (default) \code{NA}s are removed. If
  FALSE and there are \code{NA}s, the algorithm will stop with an error.}
}
\details{
The algorithm is based on sections 4.1 and 4.2 of Bernegger (1997). With rare
exceptions, the fitted \eqn{g} and \eqn{b} parameters must conform to:
\deqn{\mu = \frac{\ln(gb)(1-b)}{\ln(b)(1-gb)}}{
\mu = ln(gb)(1-b) / ln(b)(1-gb)}

where \eqn{\mu} is the empirical mean.

However, in step 2 of section 4.2, the \eqn{p} component is estimated as the
difference between the numerical integration of \eqn{x^2 f(x)} and the empirical
second moment, per equation (4.3). This is converted to \eqn{g} by
reciprocation, and convergence is tested by the difference between this
\dQuote{new} \eqn{g} and its prior value.
}
\value{
Returns a \link{list} containing:
    \item{g}{The fitted \code{g} parameter.}
    \item{b}{The fitted \code{b} parameter.}
    \item{iter}{The number of iterations used.}
    \item{err}{The absolute error between fitted the empirical mean and the
    theoretical mean given \code{g} and \code{b}.}
}

\references{
Bernegger, S. (1997) The Swiss Re Exposure Curves and the MBBEFD Distribution
    Class. \emph{ASTIN Bulletin} \strong{27}(1), 99--111.
    \doi{10.2143/AST.27.1.563208}
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}

\seealso{
\code{\link{rmb}} for random variate generation.
}

\examples{
set.seed(85L)
x <- rmb(1000, 25, 4)
MoMmb(x)
}

\keyword{distribution}
