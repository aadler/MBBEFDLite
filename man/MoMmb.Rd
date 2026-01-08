\name{mommb}
\alias{mommb}
\title{
Method of Moments Parameter Estimation for the MBBEFD distribution
}
\description{
Attempts to find the best \code{g} and \code{b} parameters which are consistent
with the first and second moments of the supplied data.
}
\usage{
mommb(x, m = FALSE, maxit = 100L, tol = NULL, na.rm = TRUE, trace = FALSE)
}
\arguments{
  \item{x}{\strong{numeric}; If \code{m} is \code{FALSE}, a vector of
  observations between 0 and 1. If \code{m} is \code{TRUE}, then a vector of
  length 2, where the first element is the first central moment (mean) of the
  \acronym{MBBEFD} distribution and the second element is the second central
  moment (variance) of the \acronym{MBBEFD} distribution.}
  \item{m}{\strong{logical}; When \code{FALSE}---the default---\code{x} is
  treated as a vector of observations. When \code{TRUE}, \code{x} is treated as
  the couplet of the distribution's first two central \strong{m}oments---
  \eqn{E[X]} and \eqn{Var[X]}.}
  \item{maxit}{\strong{integer}; maximum number of iterations.}
  \item{tol}{\strong{numeric}; tolerance. If too tight, algorithm may fail.
  Defaults to the square root of \code{.Machine$double.eps} or roughly
  \eqn{1.49\times 10^{-8}}{1.49e-8}.}
  \item{na.rm}{\strong{logical}; if \code{TRUE} (default) \code{NA}s are
  removed. If \code{FALSE}, and there are \code{NA}s, the algorithm will stop
  with an error.}
  \item{trace}{\strong{logical}; if \code{TRUE}, the fitting routine will
  print the values of \code{g} and \code{b} at each iteration \code{i}. The
  default is \code{FALSE}.}
}
\details{
The algorithm is based on sections 4.1 and 4.2 of Bernegger (1997). With rare
exceptions, the fitted \eqn{g} and \eqn{b} parameters must conform to:
\deqn{\mu = \frac{\ln(gb)(1-b)}{\ln(b)(1-gb)}}{
\mu = ln(gb)(1-b) / ln(b)(1-gb)}

subject to:

\deqn{\mu^2 \le E[x^2]\le\mu\\ p\le E[x^2]}{\mu^2 <= E[x^2] <= \mu\\
p<= E[x^2]}

where \eqn{\mu} and \eqn{\mu^2} are the \dQuote{true} first and second raw
moments, \eqn{E[x^2]} is the empirical second raw moment, and \eqn{p} is the
mass point probability of a maximal loss: \eqn{1 - F(1^{-})}.

The algorithm starts with the estimate \eqn{p = E[x^2]} as an upper bound.
However, in step 2 of section 4.2, the \eqn{p} component is estimated as the
difference between the numerical integration of \eqn{x^2 f(x)} and the empirical
second moment---\eqn{p = E[x^2] - \int x^2 f(x) dx}---as seen in equation (4.3).
This is converted to \eqn{g} by reciprocation and convergence is tested by the
difference between this new \eqn{g} and its prior value. If the new
\eqn{p \le 0}, the algorithm stops with an error.
}
\value{
Returns a \link{list} containing:
    \item{g}{The fitted \code{g} parameter.}
    \item{b}{The fitted \code{b} parameter.}
    \item{iter}{The number of iterations used.}
    \item{sqerr}{The squared error between the empirical mean and the
    theoretical mean given the fitted \code{g} and \code{b}.}
}

\note{Anecdotal evidence indicates that the results of this fitting algorithm
can be volatile, especially with fewer than a few hundred observations.}

\references{
Bernegger, S. (1997) The Swiss Re Exposure Curves and the \acronym{MBBEFD}
    Distribution Class. \emph{\acronym{ASTIN} Bulletin} \strong{27}(1), 99--111.
    \doi{10.2143/AST.27.1.563208}
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}

\seealso{
\code{\link{rmb}} for random variate generation.
}

\examples{
set.seed(85L)
x <- rmb(1000, 25, 4)
mommb(x)
}

\keyword{distribution}
