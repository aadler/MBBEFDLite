# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

findbIntMax <- 1e50

findb <- function(mu, g, tol = NULL) {

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  mm1 <- mu - 1
  gm1 <- g - 1
  lg <- log(g)

  if (abs(mm1) <= tol) {
    return(0)
  } else if (abs(mu - 1 / g) <= tol) {
    return(Inf)
  } else if (abs(mu - lg / gm1) <= tol) {
    return(1)
  } else if (abs(mu - gm1 / (lg * g)) <= tol) {
    return(1 / g)
  } else {
    errf <- function(b) {
      gb <- g * b
      (log(gb) * (1 - b) / (log(b) * (1 - gb)) - mu) ^ 2
    }
    return(optimize(errf, c(.Machine$double.eps, findbIntMax))$minimum)
  }
}

mommb <- function(x, m = FALSE, maxit = 100L, tol = NULL, na.rm = TRUE,
                  trace = FALSE) {

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  if (m) {
    if (length(x) != 2L) {
      stop("Was expecting the first and second moments but something other ",
      "than 2 parameters was passed.")
    }

    mu <- x[1L]
    mu2 <- mu ^ 2 + x[2L]
  } else {
    if (!na.rm && anyNA(x)) {
      stop("There are NAs in the data yet na.rm was passed as FALSE.")
    }

    mu <- mean(x, na.rm = na.rm)
    mu2 <- mu ^ 2 + var(x, na.rm = na.rm)
  }

  g <- 1 / mu2
  b <- findb(mu, g, tol)
  converged <- FALSE
  i <- 0L
  while (!converged && i < maxit) {
    if (trace) message("i: ", i, "\tg: ", g, "\tb: ", b)
    i <- i + 1L
    oldg <- g
    intxsqrd <- tryCatch(integrate(function(x) {x ^ 2 * dmb(x, g, b)},
                                   lower = 0, upper = 1, subdivisions = 1000L,
                                   rel.tol = tol)$value,
                         error = function(cond) {
                           simpleError(trimws(cond$message))
                         })
    if (inherits(intxsqrd, "simpleError")) {
      stop(trimws(intxsqrd$message), "; perhaps try a looser tolerance.")
    }

    newp <- mu2 - intxsqrd

    if (newp <= 0) {
      stop("Algorithm has insufficient data to converge to a method of ",
           "moments solution.")
    }

    g <- 1 / newp
    converged <- abs(oldg - g) <= tol
    b <- findb(mu, g, tol)
  }

  if ((i >= maxit && !converged) || b >= findbIntMax) {
    stop("Algorithm has insufficient data to converge to a method of ",
         "moments solution.")
  }

  list(g = g, b = b, iter = i,
       sqerr = (log(g * b) * (1 - b) / (log(b) * (1 - g * b)) - mu) ^ 2)
}
