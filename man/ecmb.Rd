\name{ecmb}
\alias{ecmb}
\title{Exposure Curve for the MBBEFD Distribution}
\description{
Returns the limited average severity at \code{x} of a random variable with an
\acronym{MBBEFD} distribution with parameters \code{g} and \code{b}.
}
\usage{
ecmb(x, g, b, c = NULL, lower.tail = TRUE)
}
\arguments{
  \item{x}{\strong{numeric}; vector of quantiles.}
  \item{g}{\strong{numeric}; (vector of) the \code{g} parameter, which is also
  the reciprocal of the probability of a maximum loss.}
  \item{b}{\strong{numeric}; (vector of) the \code{b} parameter.}
  \item{c}{\strong{numeric}; (vector of) the optional \code{c} parameter. Should
  be NULL if \code{g} and \code{b} are passed. Otherwise,
  \eqn{g = e^{(0.78 + 0.12c)c)}}{exp((0.78 + 0.12 * c) * c)} and
  \eqn{b = e^{3.1 - 0.15(1+c)c}}{b = exp(3.1 - 0.15 * (1 + c) * c))}.}
  \item{lower.tail}{\strong{logical}; if TRUE (default), percentages are of the
  total loss being less than or equal to \code{x}. Otherwise they are the
  percentage of total loss greater than x.}
}
\details{
Given random variable \eqn{X} with an \acronym{MBBEFD} distribution with
parameters \eqn{g} and \eqn{b}, the \emph{exposure curve} is defined as the
ratio of the limited average severity (\acronym{LAS}) of the variable at \eqn{x}
divided by the unlimited expected severity of the distribution. If one considers
\eqn{x} as a policy limit, then the value of \code{ecmb(x, g, b)} is the
percentage of the total expected loss which will contained within the policy
limit. If \code{lower.tail} is \code{FALSE}, the return value is the percentage
of total loss which will exceed the limit.
}
\value{
A numeric vector containing the values of the exposure curve for the passed
\code{x}, \code{b}, and \code{g} vectors. If \code{lower.tail} is \code{FALSE},
the return value is the complement of the exposure curve.
}
\references{
Bernegger, S. (1997) The Swiss Re Exposure Curves and the \acronym{MBBEFD}
    Distribution Class. \emph{\acronym{ASTIN} Bulletin} \strong{27}(1), 99--111.
    \doi{10.2143/AST.27.1.563208}
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}

\seealso{
\code{\link{dmb}} and \code{\link{pmb}} for the density and distribution.
}
\examples{
all.equal(ecmb(c(0, 1), 20, 5), c(0, 1))
ecmb(0.25, 100, 34)
}

\keyword{distribution}
