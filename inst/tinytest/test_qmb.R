# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 20 * .Machine$double.eps

p <- c(NA, NaN, -0.2, 0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1, 1.2)
omp <- 1 - p
logp <- suppressWarnings(log(p))

# Standard b & g and passing log.p and lower.tail
g <- 20
b <- 0.5
control <- c(NA, NaN, NaN, 0, 0.00076677920534173882, 0.0084122398161999845,
             0.025090980962830467, 0.074000581443776747, 0.21150410519371177,
             0.55942740861401874, 1, NaN)

expect_equal(qmb(p, g, b), control, tolerance = tol)
expect_equal(qmb(logp, g, b, log.p = TRUE), control, tolerance = tol)
expect_equal(qmb(p, g, b, lower.tail = FALSE), qmb(omp, g, b), tolerance = tol)
expect_equal(qmb(logp, g, b, lower.tail = FALSE, log.p = TRUE), qmb(omp, g, b),
             tolerance = tol)

# Nonstandard g & b
## g < 1 and b < 0
expect_true(is.nan(qmb(0.5, 0.2, 6)))
expect_true(is.nan(qmb(0.5, 1.2, -0.3)))

## g == 1 and b == 0
expect_identical(qmb(0.5, 1, 1), 0)
expect_identical(qmb(0.5, 1.3, 0), 0)

## b == 1
expect_equal(qmb(0.25, 4, 1), 0.1111111111111111, tolerance = tol)

## bg == 1
expect_equal(qmb(0.25, 4, 0.25), 0.20751874963942191, tolerance = tol)

# Test vectorized b & g
g <- c(1.2, 4, 100)
b <- c(0.001, 0.17)
control <- c(qmb(p[1L], g[1L], b[1L]),
             qmb(p[2L], g[2L], b[2L]),
             qmb(p[3L], g[3L], b[1L]),
             qmb(p[4L], g[1L], b[2L]),
             qmb(p[5L], g[2L], b[1L]),
             qmb(p[6L], g[3L], b[2L]))

expect_identical(qmb(p, g, b)[1:6], control)
