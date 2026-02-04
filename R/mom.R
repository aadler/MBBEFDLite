# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

findb <- function(mu, g, tol = NULL, maxb) {

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  if (abs(mu - 1) <= tol) {
    return(0)
  } else if (abs(mu - 1 / g) <= tol) {
    return(Inf)
  }

  gm1 <- g - 1
  lg <- log(g)

  if (abs(mu - lg / gm1) <= tol) {
    return(1)
  } else if (abs(mu - gm1 / (lg * g)) <= tol) {
    return(1 / g)
  } else {
    errf <- function(b) {
      gb <- g * b
      (log(gb) * (1 - b) / (log(b) * (1 - gb)) - mu) ^ 2
    }
    return(optimize(errf, c(.Machine$double.eps, maxb))$minimum)
  }
}

mommb <- function(x, m = FALSE, maxit = 100L, tol = NULL, na.rm = TRUE,
                  trace = FALSE, maxb = 1e3) {

  if (!is.numeric(maxit) || maxit < 1L) {
    stop("maxit must be a positive integer.")
  }

  if (!is.numeric(maxb) || maxb <= 0) {
    stop("maxb must be positive.")
  }

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  if (m) {
    if (length(x) != 2L) {
      stop("Was expecting the mean and variance but something other than 2 ",
      "parameters was passed.")
    }

    mu <- x[1L]
    mu2 <- mu ^ 2 + x[2L] # x[2L] is the variance; mu2 is E[X^2]
  } else {
    if (!na.rm && anyNA(x)) {
      stop("There are NAs in the data yet na.rm was passed as FALSE.")
    }

    mu <- mean(x, na.rm = na.rm)
    mu2 <- mu ^ 2 + var(x, na.rm = na.rm)
  }

  if (!is.finite(mu) || mu < 0 || mu > 1) {
    stop("The mean must be in [0, 1] for the MBBEFD distribution")
  }

  if (!is.finite(mu2) || mu2 > mu) {
    stop("The variance of an MBBEFD distribution must be less than or equal ",
    "to its mean.")
  }

  g <- 1 / mu2
  b <- findb(mu, g, tol, maxb)
  converged <- FALSE
  i <- 1L
  while (!converged && i <= maxit) {
    if (trace) message("i: ", i, "\tg: ", g, "\tb: ", b)
    if (is.infinite(b)) warning("Parameter b is Inf. Mean must be = 1 / g.")
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
    b <- findb(mu, g, tol, maxb)
    i <- i + 1L
  }

  if (i >= maxit && !converged) {
    stop("Algorithm failed to converge after ", maxit, " iterations. ",
         "Final change in g: ", abs(oldg - g), ". Try increasing maxit or tol.")
  }

  if (b >= 0.999 * maxb) {
    stop("Parameter b approaching maximum bound (", maxb, "). ",
         "Algorithm may not have converged properly. Try increasing maxb.")
  }

  list(g = g, b = b, iter = i,
       sqerr = (log(g * b) * (1 - b) / (log(b) * (1 - g * b)) - mu) ^ 2)
}
