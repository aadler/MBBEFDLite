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
  stopifnot(exprs = {
    length(p) >= 1
    is.logical(log.p)
  })

  n <- length(p)
  g <- rep(g, n)
  b <- rep(b, n)
  gm1 <- g - 1
  gb <- g * b
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  pc <- 1 - p
  ret <- ifelse(g < 1 | b < 0 | p < 0 | p > 1, NaN,
                ifelse(p >= 1 - 1 / g, 1,
                       ifelse(g == 1 | b == 0 | p <= 0, 0,
                              ifelse(b == 1, (1 - 1 / pc) / gm1,
                                     ifelse(gb == 1, log(pc) / log(b), 1 -
                                              (log((1 - b) / pc - 1 + gb) -
                                                 log(gm1)) / log(b)
                                     )
                              )
                       )
                )
  )

  ret
}

rmb <- function(n, g, b) {
  qmb(runif(n), g, b)
}
