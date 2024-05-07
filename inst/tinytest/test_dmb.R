# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps
x <- c(NaN, NA_real_, seq(-0.05, 1.05, 0.05))

# Standard b & g and passing log
g <- 20
b <- 0.5
control <- c(NaN, NA, 0, 13.169796430638961, 4.8885512023122937,
             2.5263813046917969, 1.5385762756474961, 1.0338586295606644,
             0.74178454647271053, 0.55773712910156925, 0.43433663080658902,
             0.34759405724528358, 0.28430720436213047, 0.23672338994222369,
             0.20004716240988563, 0.1711831300348825, 0.14806084939852232,
             0.12925298240848238, 0.1137496715734246, 0.10082014212059009,
             0.089925016722583739, 0.080659166235059657, 0.072713521352461077,
             0, 0)

expect_equal(dmb(x, g, b), control, tolerance = tol)
expect_equal(dmb(x, g, b, log = TRUE), log(control), tolerance = tol)

# Nonstandard g & b
## g < 1 and b < 0
expect_true(is.nan(dmb(0.5, 0.2, 6)))
expect_true(is.nan(dmb(0.5, 1.2, -0.3)))

## g == 1 and b == 0
expect_identical(dmb(0.5, 1, 1), 0)
expect_identical(dmb(0.5, 1.3, 0), 0)

## b == 1
expect_equal(dmb(0.5, 1.2, 1), 0.16528925619834711, tolerance = tol)

## bg == 1
expect_equal(dmb(0.5, 5, 0.2), 0.71976251555360038, tolerance = tol)

# Test vectorized b & g
g <- c(1.2, 4, 100)
b <- c(0.001, 0.17)
control <- c(dmb(x[1L], g[1L], b[1L]),
             dmb(x[2L], g[2L], b[2L]),
             dmb(x[3L], g[3L], b[1L]),
             dmb(x[4L], g[1L], b[2L]),
             dmb(x[5L], g[2L], b[1L]),
             dmb(x[6L], g[3L], b[2L]))

expect_identical(dmb(x, g, b)[1:6], control)
