# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps

# Test Functionality: m = FALSE
set.seed(138L)
x1 <- rmb(1e5, 150, 250)
testx1 <- mommb(x1)
expect_true(abs(testx1$g / 150 - 1) <= 0.02)
expect_true(abs(testx1$b / 250 - 1) <= 0.05)
expect_equal(testx1$sqerr,
             (log(testx1$g * testx1$b) * (1 - testx1$b) /
               (log(testx1$b) * (1 - testx1$g * testx1$b)) - mean(x1)) ^ 2,
             tolerance = tol)

# Test Functionality: m = TRUE
z <- c(mean(x1), var(x1))
testz <- mommb(z, m = TRUE)
expect_equal(testz$g, testx1$g, tolerance = tol)
expect_equal(testz$b, testx1$b, tolerance = tol)

# Testing simple findb branch
expect_identical(mommb(c(1, 1, 1))$b, 0)                    # findb => 0
x <- c(0.9, 0.9, 0.9, 0.93785026012074624)
expect_silent(mommb(x, tol = 0.002))                        # findb => 1

# Malformed input created solely to test branches in findb. Values found using
# ChatGPT
expect_error(suppressWarnings(mommb(c(0, 0.5, 1))))         # findb => Inf
expect_warning(mommb(c(0, 0.5, 1)), "Parameter b is Inf")

a <- 2 * log(2)
d <- sqrt(2 - a ^ 2)
expect_error(mommb(c(a - d, a, a + d)))    # mean must be in [0, 1]

# Test Error Trapping
expect_error(mommb(x, maxit = "3L"), "maxit must be a positive integer.")
expect_error(mommb(x, maxb = "3L"), "maxb must be positive.")
expect_error(mommb(x, maxit = 3L), "Try increasing maxit or tol")
expect_error(mommb(c(0.2, 0.3), m = TRUE), "less than or equal")
expect_error(mommb(x, maxb = 2e-16), "Parameter b approaching")

set.seed(76L)
expect_error(mommb(rmb(10, 10, 9)), "insufficient data")
expect_error(mommb(rmb(10, 10, 9), tol = 1e-16), "looser tolerance")
expect_error(mommb(NA_real_, na.rm = FALSE), "passed as FALSE")
expect_error(mommb(rmb(10, 10, 9), m = TRUE), "other than 2")

# Test trace
expect_message(mommb(x, trace = TRUE), "i: 1")

# Test tol in finddb
expect_equal(MBBEFDLite:::findb(0.7, 4, maxb = 6),
             MBBEFDLite:::findb(0.7, 4, tol = sqrt(.Machine$double.eps),
                                maxb = 6), tolerance = tol)
