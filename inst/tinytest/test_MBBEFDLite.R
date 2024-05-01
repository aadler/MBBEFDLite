# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps
x <- c(NaN, NA_real_, seq(-0.05, 1.05, .05))

############################### Testing dmb ####################################
# Standard b & g and passing log
g <- 20
b <- 0.5
dcontrolx <- c(NaN, NA, 0, 13.169796430638961, 4.8885512023122937,
               2.5263813046917969, 1.5385762756474961, 1.0338586295606644,
               0.74178454647271053, 0.55773712910156925, 0.43433663080658902,
               0.34759405724528358, 0.28430720436213047, 0.23672338994222369,
               0.20004716240988563, 0.1711831300348825, 0.14806084939852232,
               0.12925298240848238, 0.1137496715734246, 0.10082014212059009,
               0.089925016722583739, 0.080659166235059657, 0.072713521352461077,
               0, 0)

expect_equal(dmb(x, g, b), dcontrolx, tolerance = tol)
expect_equal(dmb(x, g, b, TRUE), log(dcontrolx), tolerance = tol)

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
dcontrolx <- c(dmb(x[1L], g[1L], b[1L]),
              dmb(x[2L], g[2L], b[2L]),
              dmb(x[3L], g[3L], b[1L]),
              dmb(x[4L], g[1L], b[2L]),
              dmb(x[5L], g[2L], b[1L]),
              dmb(x[6L], g[3L], b[2L]))

expect_identical(dmb(x, g, b)[1:6], dcontrolx)

############################### Testing pmb ####################################
# Standard b & g
g <- 20
b <- 0.5

pcontrolx <- c(NaN, NA, 0, 0, 0.40120963545198918, 0.57693371329906551,
               0.67551641239099069, 0.73858045887054991, 0.78236907383302978,
               0.81453124867180904,  0.83914170955380674, 0.85857029312602218,
               0.8742891819815175, 0.88726116159525426, 0.89814240493843012,
               0.90739552971536974, 0.91535606623643084, 0.92227329508477796,
               0.92833628856856798,  0.93369105420497067, 0.93845213614173995,
               0.94271065786700514,  0.94654001782806174, 1, 1)

expect_equal(pmb(x, g, b), pcontrolx, tolerance = tol)
expect_equal(pmb(x, g, b, lower.tail = FALSE), 1 - pcontrolx, tolerance = tol)
expect_equal(pmb(x, g, b, log.p = TRUE), log(pcontrolx), tolerance = tol)
expect_equal(pmb(x, g, b, lower.tail = FALSE, log.p = TRUE), log(1 - pcontrolx),
             tolerance = tol)

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

# Test vectorized b & g
g <- c(1.2, 4, 100)
b <- c(0.001, 0.17)
pcontrolx <- c(pmb(x[1L], g[1L], b[1L]),
               pmb(x[2L], g[2L], b[2L]),
               pmb(x[3L], g[3L], b[1L]),
               pmb(x[4L], g[1L], b[2L]),
               pmb(x[5L], g[2L], b[1L]),
               pmb(x[6L], g[3L], b[2L]))

expect_identical(pmb(x, g, b)[1:6], pcontrolx)


############################### Testing qmb ####################################
p <- c(NA, NaN, -0.2, 0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1, 1.2)
omp <- 1 - p
logp <- suppressWarnings(log(p))
# Standard b & g
g <- 20
b <- 0.5
qcontrolp <- c(NA, NaN, NaN, 0, 0.00076677920534173882, 0.0084122398161999845,
              0.025090980962830467, 0.074000581443776747, 0.21150410519371177,
              0.55942740861401874, 1, NaN)

expect_equal(qmb(p, g, b), qcontrolp, tolerance = tol)
expect_equal(qmb(logp, g, b, log.p = TRUE), qcontrolp, tolerance = tol)

expect_equal(qmb(p, g, b, lower.tail = FALSE), qmb(omp, g, b), tolerance = tol)
# Ask MBBEFD team why they return NaNs in case below
expect_equal(qmb(logp, g, b, lower.tail = FALSE, log.p = TRUE), qmb(omp, g, b),
             tolerance = tol)

# Edge case b & g and vector

############################### Testing rmb ####################################
set.seed(9712L)
u <- runif(100L)
## Standard b & g
g <- 20
b <- 0.5
rcontrol <- qmb(u, g, b)
set.seed(9712L)
expect_identical(rmb(100L, g, b), rcontrol)

# Edge case b & g and vector
