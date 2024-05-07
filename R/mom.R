# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

mommb <- function(x, maxit = 100L, tol = .Machine$double.eps ^ 0.5,
                  na.rm = TRUE) {

  if (anyNA(x) && !na.rm) {
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
        abs(log(g * b) * (1 - b) / (log(b) * (1 - g * b)) - mu)
      }
      bFit <- optimize(errf, c(0, 1e10))
    }

    bFit$minimum
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
    mu2part <- tryCatch(integrate(\(x) x ^ 2 * dmb(x, g, b), lower = 0,
                                  upper = 1, subdivisions = 1000L,
                                  rel.tol = tol)$value,
                        error = \(cond) simpleError(trimws(cond$message)))
    if (inherits(mu2part, "simpleError")) {
      stop(trimws(mu2part$message), "; please try a looser tolerance")
    }

    ppart <- mu2 - mu2part
    if (ppart <= 0) {
      stop("Algorithm has insufficient data to converge to a method of ",
           "moments solution.")
    }

    g <- 1 / ppart
    converged <- abs(oldg - g) <= tol
    b <- findb(mu, g, tol)
  }

  list(g = g, b = b, iter = i,
       err = abs(log(g * b) * (1 - b) / (log(b) * (1 - g * b)) - mu))
}
