# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

dmb <- function(x, g, b, log = FALSE) {
  .Call(dmb_c, as.double(x), as.double(g), as.double(b), as.logical(log))
}

pmb <- function(q, g, b, lower.tail = TRUE, log.p = FALSE) {
  .Call(pmb_c, as.double(q), as.double(g), as.double(b), as.logical(lower.tail),
        as.logical(log.p))
}

qmb <- function(p, g, b, lower.tail = TRUE, log.p = FALSE) {
  .Call(qmb_c, as.double(p), as.double(g), as.double(b), as.logical(lower.tail),
        as.logical(log.p))
}

rmb <- function(n, g, b) {
  if (length(n) > 1) n <- length(n)
  .Call(rmb_c, as.double(n), as.double(g), as.double(b))
}

ecmb <- function(x, g, b, lower.tail = TRUE) {
  .Call(ecmb_c, as.double(x), as.double(g), as.double(b),
        as.logical(lower.tail))
}
