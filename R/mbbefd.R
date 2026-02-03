# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

dmb <- function(x, g, b, c = NULL, log = FALSE) {
  gb <- checkgcb(g, b, c)
  .Call(dmb_c, as.double(x), as.double(gb$g), as.double(gb$b), as.logical(log))
}

pmb <- function(q, g, b, c = NULL, lower.tail = TRUE, log.p = FALSE) {
  gb <- checkgcb(g, b, c)
  .Call(pmb_c, as.double(q), as.double(gb$g), as.double(gb$b),
        as.logical(lower.tail), as.logical(log.p))
}

qmb <- function(p, g, b, c = NULL, lower.tail = TRUE, log.p = FALSE) {
  gb <- checkgcb(g, b, c)
  .Call(qmb_c, as.double(p), as.double(gb$g), as.double(gb$b),
        as.logical(lower.tail), as.logical(log.p))
}

rmb <- function(n, g, b, c = NULL) {
  gb <- checkgcb(g, b, c)
  if (length(n) > 1) n <- length(n)
  .Call(rmb_c, as.double(n), as.double(gb$g), as.double(gb$b))
}

ecmb <- function(x, g, b, c = NULL, lower.tail = TRUE) {
  gb <- checkgcb(g, b, c)
  .Call(ecmb_c, as.double(x), as.double(gb$g), as.double(gb$b),
        as.logical(lower.tail))
}

checkgcb <- function(g, b, c) {
  if (is.null(c)) {
    return(list(g = g, b = b))
  } else if (missing(g) && missing(b)) {
    return(c2gb(c))
  }
  stop("A c parameter was passed together with either a g or b parameter.")
}

c2gb <- function(c) {
  list(g = exp((0.78 + 0.12 * c) * c), b = exp(3.1 - 0.15 * (1 + c) * c))
}
