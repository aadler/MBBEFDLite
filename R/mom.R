# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

dilog <- function(x) {
  .Call(dilog_c, as.double(x))
}

# findb, mugb, and cvgb are internal only so no standard checks are needed as
# that will happen inside of mommb.
findb <- function(mu, g, maxb) {
  eps <- .Machine$double.eps

  if (abs(mu - 1) < eps) {
    return(0)
  } else if (abs(mu - 1 / g) < eps) {
    return(Inf)
  }

  gm1 <- g - 1
  lg <- log(g)

  if (abs(mu - lg / gm1) < eps) {
    return(1)
  } else if (abs(mu - gm1 / (lg * g)) < eps) {
    return(1 / g) # nocov With mean restricted to [0, 1], this may be impossible
  } else {
    errf <- function(b) {
      gb <- g * b
      (log(gb) * (1 - b) / (log(b) * (1 - gb)) - mu) ^ 2
    }
    return(optimize(errf, c(.Machine$double.eps, maxb))$minimum)
  }
}

mugb <- function(g, b) {
  eps <- .Machine$double.eps
  gm1 <- g - 1
  bc <- 0.5 - b + 0.5
  bm1 <- -bc
  gb <- g * b
  gbc <- 0.5 - gb + 0.5
  lb <- log(b)

  suppressWarnings(
    ifelse(b < 0 | gm1 < 0, NaN,
         ifelse(abs(b * gm1) < eps, 1,
                ifelse(abs(bm1) < eps, log(g) / gm1,
                       ifelse(abs(gbc) < eps, bm1 / lb,
                              log(gb) * bc / (lb * gbc)))))
  )
}

cvgb <- function(g, b) {
  eps <- .Machine$double.eps
  gm1 <- g - 1
  bc <- 0.5 - b + 0.5   # b-complement (1 - b)
  bg <- b * g
  bgm1 <- bg - b
  lb <- log(b)
  lg <- log(g)
  lbg <- lb + lg

  cvp <- suppressWarnings(
         ifelse(abs(b * gm1) < eps, 0.5,
           ifelse(abs(bc) < eps, (gm1 - log(g)) / log(g) ^ 2,
             ifelse(abs(bg - 1) < eps, (lb * b + bc) / bc ^ 2,
               ifelse(b > 0 & b < 1, (lb * log(g * bc / gm1) -
                      dilog(1 - bc / bgm1) + dilog(b - bc / gm1)) /
                        (-bc / (bg - 1) * lbg ^ 2),
                  ifelse(b > 1, (lbg * log(1 + bc / (bg -1)) +
                                   dilog(bc / bgm1) -
                                   dilog(1 - (bg - 1) / gm1)) /
                           (-bc / (bg - 1) * lbg ^ 2), NaN)))))
  )
  suppressWarnings(sqrt(2 * cvp - 1))
}

mommb <- function(x, m = FALSE, tol = NULL, na.rm = TRUE, opts = list()) {

  nopts <- names(opts)

  if (!("alg" %in% nopts)) {
    opts$alg <- "EM"
  } else if (!(toupper(opts$alg) %in% c("EM", "LS"))) {
    stop("Algorithm must be either 'EM' (expectation-maximization) or 'LS' ",
         "(line search). See documentation.")
  }

  if (!("maxit" %in% nopts)) {
    opts$maxit <- 100L
  } else if (!is.numeric(opts$maxit) || opts$maxit < 1L) {
    stop("maxit must be a positive integer.")
  }

  if (!("maxb" %in% nopts)) {
    opts$maxb <- 1e6
  } else if (!is.numeric(opts$maxb) || opts$maxb <= 0 || !is.finite(opts$maxb)) {
    stop("maxb must be positive and finite.")
  }

  if (!("minb" %in% nopts)) {
    opts$minb <- 1e-10
  } else if (!is.numeric(opts$minb) || opts$minb <= 0 || !is.finite(opts$minb)) {
    stop("minb must be positive and finite.")
  }

  if (!("maxg" %in% nopts)) {
    opts$maxg <- 1e6
  } else if (!is.numeric(opts$maxg) || opts$maxg <= 0 || !is.finite(opts$maxg)) {
    stop("maxg must be positive and finite.")
  }

  if (!("ming" %in% nopts)) {
    opts$ming <- 1 + 1e-10
  } else if (!is.numeric(opts$ming) || opts$ming <= 1 || !is.finite(opts$ming)) {
    stop("ming must be finite and strictly greater than 1.")
  }

  if (opts$minb >= opts$maxb) {
    stop("minb must be strictly less than maxb.")
  }

  if (opts$ming >= opts$maxg) {
    stop("ming must be strictly less than maxg.")
  }

  if (!("trace" %in% nopts)) {
    opts$trace <- FALSE
  } else if (!is.logical(opts$trace)) {
      stop("trace must be a logical (TRUE or FALSE).")
  }

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  if (m) {
    if (length(x) != 2L) {
      stop("Was expecting the mean and variance but something other than 2 ",
           "parameters was passed.")
    }
    mu <- x[1L]
    mu2 <- mu ^ 2 + x[2L] # x[2L] is the variance; mu2 is E[X^2]
    s <- sqrt(x[2L])
  } else {
    if (!na.rm && anyNA(x)) {
      stop("There are NAs in the data yet na.rm was passed as FALSE.")
    }
    mu <- mean(x, na.rm = na.rm)
    mu2 <- mu ^ 2 + var(x, na.rm = na.rm)
    s <- sd(x, na.rm = na.rm)
  }

  if (!is.finite(mu) || mu < 0 || mu > 1) {
    stop("The mean must be in [0, 1] for the MBBEFD distribution")
  }

  if (!is.finite(mu2) || mu2 > mu || s > sqrt(mu * (1 - mu))) {
    stop("The variance of an MBBEFD distribution must be less than or equal ",
          "to its mean.")
  }

  if (opts$alg == "LS") {                # Line Search
    if (opts$trace) message("trace is ignored for the Bernegger algorithm")

    cv <- s / mu

    getb <- function(g) {
      f <- function(b) (mugb(g, b) - mu) ^ 2
      optimize(f, lower = opts$minb, upper = opts$maxb, tol = 1e-12)$minimum
    }

    getg <- function(g) (cvgb(g, getb(g)) - cv) ^ 2

    g <- optimize(getg, lower = opts$ming, upper = opts$maxg,
                  tol = 1e-12)$minimum
    b <- getb(g)

    gTries <- 1
    if (g >= 0.999 * opts$maxg) {
      gTries <- gTries + 1
      newG <- sqrt(opts$maxg)
      message("Parameter g approaching maximum bound (", opts$maxg, "). ",
      "Trying again with square root: ", newG)
      g <- optimize(getg, lower = opts$ming, upper = newG, tol = 1e-12)$minimum
      b <- getb(g)
      if (g >= 0.999 * newG) {
        stop("Algorithm has insufficient data to converge to a method of ",
             "moments solution.")
      }
    }
    iter <- gTries
  } else {                               # Expectation-Maximization
    g <- 1 / mu2
    b <- findb(mu, g, opts$maxb)
    converged <- FALSE
    i <- 0L
    while (!converged && i <= opts$maxit) {
      i <- i + 1L
      if (opts$trace) message("i: ", i, "\tg: ", g, "\tb: ", b)
      if (is.infinite(b)) warning("Parameter b is Inf. Mean must be = 1 / g.") # nolint nonportable_path_linter
      oldg <- g

      # Decompose second moment: E[X²] = p + ∫x²f(x)dx (eq. 4.3)
      # where p is the point mass probability at x=1
      intxsqrd <- tryCatch(integrate(function(x) {x ^ 2 * dmb(x, g, b)},
                                     lower = 0, upper = 1, subdivisions = 1001L,
                                     rel.tol = tol)$value,
                           error = function(cond) {
                             simpleError(trimws(cond$message))
                           })

      if (inherits(intxsqrd, "simpleError")) {
        stop(trimws(intxsqrd$message), "; perhaps try a looser tolerance.")
      }

      # Extract point mass probability p from residual
      newp <- mu2 - intxsqrd

      if (newp <= 0) {
        stop("Algorithm has insufficient data to converge to a method of ",
             "moments solution.")
      }

      # Update g = 1/p (section 4.1)
      g <- 1 / newp
      converged <- abs(oldg - g) <= tol
      b <- findb(mu, g, opts$maxb)
    }
    iter <- i
    if (i >= opts$maxit && !converged) {
      stop("Algorithm failed to converge after ", opts$maxit, " iterations. ",
           "Final change in g: ", abs(oldg - g), ". Try increasing maxit ",
           "or tol.")
    }
  }

  if (b >= 0.999 * opts$maxb) {
    stop("Parameter b approaching maximum bound (", opts$maxb, "). ",
         "Algorithm may not have converged properly. Try increasing maxb.")
  }

  list(g = g, b = b, iter = iter,
       sqerr = (log(g * b) * (1 - b) / (log(b) * (1 - g * b)) - mu) ^ 2)
}
