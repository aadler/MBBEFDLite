# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps
x <- c(NaN, NA_real_, seq(-0.05, 1.05, 0.05))
g <- 20
b <- 0.5

control <- c(NaN, NA, 0, 0, 0.20767369695142743, 0.34348858177475161,
             0.44365276558828731,  0.5224559528656193, 0.58702189080549916,
             0.64142025322000318,  0.68819605267843664, 0.72904736047448293,
             0.76516384953514172,  0.79741144637789729, 0.82644008446002115,
             0.85275005076754173,  0.87673466414677226, 0.89870874607807527,
             0.91892820692329635,  0.93760388008162843, 0.95491151921550443,
             0.97099916872996539,  0.9859926945290427, 1, 1)

expect_equal(ecmb(x, g, b), control, tolerance = tol)
expect_equal(ecmb(x, g, b, lower.tail = FALSE), 1 - control, tolerance = tol)

# Nonstandard g & b
## g < 1 and b < 0
expect_true(is.nan(ecmb(0.5, 0.2, 6)))
expect_true(is.nan(ecmb(0.5, 1.2, -0.3)))

## g == 1 and b == 0
expect_identical(ecmb(0.37, 1, 1), 0.37)
expect_identical(ecmb(0.9, 1.3, 0), 0.9)

## b == 1
expect_equal(ecmb(0.25, 4, 1), 0.40367746102880203, tolerance = tol)

## bg == 1
expect_equal(ecmb(0.5, 5, 0.2), 0.69098300562505266, tolerance = tol)

# Test vectorized b & g
g <- c(1.2, 4, 100)
b <- c(0.001, 0.17)

control <- c(ecmb(x[3L], g[1L], b[1L]),
             ecmb(x[4L], g[2L], b[2L]),
             ecmb(x[5L], g[3L], b[1L]),
             ecmb(x[6L], g[1L], b[2L]),
             ecmb(x[7L], g[2L], b[1L]),
             ecmb(x[8L], g[3L], b[2L]))

expect_equal(ecmb(x[3:8], g, b), control, tolerance = tol)
