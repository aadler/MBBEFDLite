% Copyright (c) 2024, Avraham Adler All rights reserved
% SPDX-License-Identifier: MPL-2.0+

\name{dmb}
\alias{dmb}
\alias{pmb}
\alias{qmb}
\alias{rmb}
\title{The MBBEFD Distribution}
\description{
Density, distribution function, quantile function and random generation for the
\acronym{MBBEFD} distribution with parameters \code{g} and \code{b}.
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
  \item{g}{\strong{numeric}; (vector of) Bernegger's \code{g} parameter, which
  is also the reciprocal of the probability of a maximum loss.}
  \item{b}{\strong{numeric}; (vector of) Bernegger's \code{b} parameter.}
  \item{log, log.p}{\strong{logical}; if TRUE, probabilities p are given as
  log(p).}
  \item{lower.tail}{\strong{logical}; if TRUE (default), probabilities are
  \eqn{P[X\leq x]}{P[X ≤ x]} otherwise \eqn{P[X > x]}}
}
\details{
The \acronym{MBBEFD}  class of curves are defined in Bernegger (1997) and are
often used to model insurance risk. The density is defined on the semi-open
interval \eqn{[0, 1)} and the distribution and quantile functions are defined on
the closed interval \eqn{[0, 1]}.
}
\value{
\code{dnorm} gives the density, \code{pnorm} gives the distribution function,
\code{qnorm} gives the quantile function, and \code{rnorm} generates random
deviates.

The length of the result is determined by \code{n} for \code{rnorm}, and is the
length of \code{x}, \code{p}, or \code{q} as appropriate for the other
functions.

Numerical arguments other than \code{n} are recycled to the length of the
result. Logical arguments should be of length 1.
}
\references{
Bernegger, S. (1997) The Swiss Re Exposure Curves and the MBBEFD Distribution
    Class. \emph{\acronym{ASTIN} Bulletin} \strong{27}(1), 99--111.
    \doi{10.2143/AST.27.1.563208}
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}

\note{
This package follows Bernegger's convention that the density function does not
exists at 1. This differs from the \pkg{mbbefd} package.
}

\seealso{
\code{\link{MoMmb}} for parameter estimation.
}
\examples{
all.equal(dmb(0.5, 1, 0), 0)
dmb(0.2, 20, 5)
pmb(0.98, 25, 4)
qmb(0.98, 25, 4) == 1
all.equal(qmb(pmb(0.98, 25, 4), 25, 4), 0.98)
set.seed(45)
rmb(3, 4, 12)
set.seed(45)
rmb(99:101, 4, 12) # Should equal previous call
}

\keyword{distribution}
