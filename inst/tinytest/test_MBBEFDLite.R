# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps

############################### Testing dmb ####################################
x <- c(0.069360915804281831, 0.81777519872412086, 0.94262173236347735,
       0.26938187633641064, 0.16934812325052917)
y <- c(-1.2, -1.0, 0.0,  0.5,  0.9,  1.0,  1.2, NaN, NA)

# Standard b & g and passing log
g <- 20
b <- 0.5
controlx <- c(3.6876104357149244, 0.096737198074403216, 0.073812924647130657,
             0.6609736077485624, 1.3068098071845686)
controly <- c(0, 0, 13.169796430638961, 0.23672338994222369,
              0.080659166235059657, 0, 0, NaN, NA)

expect_equal(dmb(x, g, b), controlx, tolerance = tol)
expect_equal(dmb(x, g, b, TRUE), log(controlx), tolerance = tol)
expect_equal(dmb(y, g, b), controly, tolerance = tol)
expect_equal(dmb(y, g, b, TRUE), log(controly), tolerance = tol)

# Nonstandard g & b
## g < 1 and b < 0
expect_true(is.nan(dmb(0.5, 0.2, 6)))
expect_true(is.nan(dmb(0.5, 1.2, -.3)))

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
controlx <- c(dmb(x[1L], g[1L], b[1L]),
              dmb(x[2L], g[2L], b[2L]),
              dmb(x[3L], g[3L], b[1L]),
              dmb(x[4L], g[1L], b[2L]),
              dmb(x[5L], g[2L], b[1L]))

expect_identical(dmb(x, g, b), controlx)

############################### Testing pmb ####################################
# Standard b & g
g <- 20
b <- 0.5
controlx <- c(0.48341340954506551, 0.93544650077891778, 0.94599947678770424,
              0.79594064206174842, 0.7029515920261491)
controly <- c(0, 0, 0, 0.88726116159525426, 0.94271065786700514, 1, 1, NaN, NA)

expect_identical(pmb(x, g, b), controlx)
expect_identical(pmb(x, g, b, lower.tail = FALSE), 1 - controlx)
expect_identical(pmb(x, g, b, log.p = TRUE), log(controlx))
expect_identical(pmb(x, g, b, lower.tail = FALSE, log.p = TRUE),
                 log(1 - controlx))
expect_identical(pmb(y, g, b), controly)
expect_identical(pmb(y, g, b, lower.tail = FALSE), 1 - controly)
expect_identical(pmb(y, g, b, log.p = TRUE), log(controly))
expect_identical(pmb(y, g, b, lower.tail = FALSE, log.p = TRUE),
                 log(1 - controly))

# Nonstandard g & b
## g < 1 and b < 0
expect_true(is.nan(pmb(0.5, 0.2, 6)))
expect_true(is.nan(pmb(0.5, 1.2, -.3)))

## g == 1 and b == 0
expect_identical(pmb(0.5, 1, 1), 0)
expect_identical(pmb(0.5, 1.3, 0), 0)

## b == 1
expect_equal(pmb(0.5, 1.2, 1), 0.090909090909090939, tolerance = tol)

## bg == 1
expect_equal(pmb(0.5, 5, 0.2), 0.55278640450004213, tolerance = tol)

############################### Testing qmb ####################################
p <- c(-0.2, 0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1, 1.2, NA, NaN)

# Standard b & g
g <- 20
b <- 0.5
controlp <- c(NaN, 0, 0.00076677920534173882, 0.0084122398161999845,
              0.025090980962830467, 0.074000581443776747, 0.21150410519371177,
              0.55942740861401874, 1, NaN, NA, NaN)

expect_equal(qmb(p, g, b), controlp, tolerance = tol)
expect_equal(qmb(p, g, b, lower.tail = FALSE), qmb(1 - p, g, b),
             tolerance = tol)
expect_equal(qmb(log(p), g, b, log.p = TRUE), controlp, tolerance = tol)
# Ask MBBEFD team why they return NaNs in case below
expect_equal(qmb(log(p), g, b, lower.tail = FALSE, log.p = TRUE),
             qmb(1 - p, g, b), tolerance = tol)

############################### Testing rmb ####################################
set.seed(9712L)
u <- runif(10L)
## Standard b & g
g <- 20
b <- 0.5
control <- qmb(u, g, b)
set.seed(9712L)
expect_identical(rmb(10L, g, b), control)
