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
mommb(x, m = FALSE, tol = NULL, na.rm = TRUE, opts = list())
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
  \item{tol}{\strong{numeric}; tolerance of the expectation-maximization
  algorithm. If too tight, algorithm may fail. Defaults to the square root of
  \code{.Machine$double.eps} or roughly \eqn{1.49\times 10^{-8}}{1.49e-8}---see
  Details.}
  \item{na.rm}{\strong{logical}; if \code{TRUE} (default) \code{NA}s are
  removed. If \code{FALSE}, and there are \code{NA}s, the algorithm will stop
  with an error.}
  \item{opts}{\link{list}; Configuration options including:
    \itemize{
      \item \code{alg}: \strong{character}; either "\acronym{EM}", the
      expectation-maximization algorithm of the package author or
      "\acronym{LS}", a form of the grid-search algorithm of Bernegger (2026),
      implemented as a nested set of 1D line-search optimizations.
      \item \code{maxit}: \strong{integer}; maximum number of iterations for the
      \acronym{EM} algorithm. Ignored for the \acronym{LS} algorithm.
      \item \code{maxb}: \strong{numeric}; the upper bound of the \code{b}
      parameter for fitting purposes. Used in both algorithms. Defaults to 1e6.
      \item \code{minb}: \strong{numeric}; the lower bound of the \code{b}
      parameter for fitting purposes. Only used in \acronym{LS} algorithm. Must
      be positive and less than \code{maxb}. Defaults to 1e-10.
      \item \code{maxg}: \strong{numeric}; the upper bound of the \code{g}
      parameter for fitting purposes. Only used in \acronym{LS} algorithm.
      Must be positive. Defaults to 1e6.
      \item \code{ming}: \strong{numeric}; the lower bound of the \code{g}
      parameter for fitting purposes. Only used in \acronym{LS} algorithm.
      Must be strictly greater than 1 and less than \code{maxg}.
      Defaults to 1 + 1e-10.
      \item \code{trace}: \strong{logical}; if \code{TRUE}, the \acronym{EM}
      algorithm will print the values of \code{g} and \code{b} at each iteration
      \code{i}. Ignored with a message for the \acronym{LS} algorithm. The
      default is \code{FALSE}.
    }
  }
}
\details{
There are two fitting algorithms.
\subsection{Expectation-Maximization}{
The default is an expectation-maximization form based on sections 4.1 and 4.2 of
Bernegger (1997). With rare exceptions, the fitted \eqn{g} and \eqn{b}
parameters must conform to:
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
\eqn{p \le 0}, the algorithm stops with an error.}
\subsection{Line Search}{
Bernegger (2026) in Algorithm 3 (Appendix C) describes a grid-search algorithm
for converging on \eqn{g} and \eqn{b}. The original algorithm looks for the sets
of parameters which return the appropriate mean and coefficient of variation
(\acronym{CV}), the latter of which can be expressed in closed form using the
dilogarithm function. However, instead of the 50,000 point grid suggested in the
paper, this package implements the algorithm as a nested set of calls to the
one-dimensional optimization algorithm of \code{optimize}. This is significantly
faster and has the benefit of returning a value even when no zeros can be found.
The algorithm basically uses the ranges suggested in the paper, which for
\eqn{g} caps at 1e6. Sometimes, the algorithm converges to the upper bound. In
this case, the package implementation will try once more using a max of the
square root of the original \code{maxg} instead. This can allow
convergence---often to the same point of the \acronym{EM} algorithm, given that
converged. When this happens, the \code{iter} value will be \code{2} instead of
\code{1}.}
}
\value{
Returns a \link{list} containing:
    \item{g}{The fitted \code{g} parameter.}
    \item{b}{The fitted \code{b} parameter.}
    \item{iter}{For the \acronym{EM} algorithm, the number of iterations used.
        For the \acronym{LS} algorithm, the number of attempts (1 or 2 if the
        retry was needed; see \strong{Details}).}
    \item{sqerr}{The squared error between the empirical mean and the
    theoretical mean given the fitted \code{g} and \code{b}. This does not make
    sense for some of the special branch cases, like \eqn{b = 0,\; g = 1,} etc.}
}

\note{Anecdotal evidence indicates that parameter estimates from either
fitting algorithm can be volatile when sample sizes are small (fewer than
a few hundred observations).}

\references{
Bernegger, Stefan. (1997) The Swiss Re Exposure Curves and the \acronym{MBBEFD}
    Distribution Class. \emph{\acronym{ASTIN} Bulletin} \strong{27}(1), 99--111.
    \doi{10.2143/AST.27.1.563208}

Bernegger, Stefan. (2026) Properties of the \acronym{MBBEFD} Distribution
    Classes.
    \url{https://www.researchgate.net/publication/400516019_Properties_of_the_MBBEF_D_Distribution_Classes}
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}

\seealso{
\code{\link{rmb}} for random variate generation.
}

\examples{
set.seed(85L)
x <- rmb(1000, 25, 4)
mommb(x)
mommb(x, opts = list(alg = "LS"))
}

\keyword{distribution}
