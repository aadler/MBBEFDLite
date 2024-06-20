# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

mommb <- function(x, maxit = 100L, tol = .Machine$double.eps ^ 0.5,
                  na.rm = TRUE) {

  if (!na.rm && anyNA(x)) {
    stop("There are NAs in the data yet na.rm was passed as FALSE.")
  }

  findb <- function(mu, g, tol) {
    if (abs(mu - 1) <= tol) {
      return(0)
    } else if (abs(mu - (g - 1) / (log(g) * g)) <= tol) {
      return(1 / g)
    } else if (abs(mu - log(g) / (g - 1)) <= tol) {
      return(1)
    } else if (abs(mu - 1 / g) <= tol) {
      return(Inf)
    } else {
      errf <- function(b) {
        (log(g * b) * (1 - b) / (log(b) * (1 - g * b)) - mu) ^ 2
      }
      return(optimize(errf, c(0, 1e100))$minimum)
    }
  }

  mu <- mean(x, na.rm = na.rm)
  mu2 <- mean(x ^ 2, na.rm = na.rm)
  g <- 1 / mu2
  b <- findb(mu, g, tol)
  converged <- FALSE
  i <- 0L
  while (!converged && i < maxit) {
    i <- i + 1L
    oldg <- g
    intxsqrd <- tryCatch(integrate(\(x) x ^ 2 * dmb(x, g, b), lower = 0,
                                   upper = 1, subdivisions = 1000L,
                                   rel.tol = tol)$value,
                         error = \(cond) simpleError(trimws(cond$message)))
    if (inherits(intxsqrd, "simpleError")) {
      stop(trimws(intxsqrd$message), "; perhaps try a looser tolerance.")
    }

    newp <- mu2 - intxsqrd
    if (newp > 0) {
      g <- 1 / newp
    } else {
      g <- g * 3 # Force restart keeping 1 / mu2 as upper limit of p
    }
    converged <- abs(oldg - g) <= tol
    b <- findb(mu, g, tol)
  }

  if ((i >= maxit && !converged) || b == 1e100) {
    stop("Algorithm has insufficient data to converge to a method of ",
         "moments solution.")
  }

  list(g = g, b = b, iter = i,
       sqerr = (log(g * b) * (1 - b) / (log(b) * (1 - g * b)) - mu) ^ 2)
}
