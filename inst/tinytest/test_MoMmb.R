# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps

# Test Functionality: m = FALSE
set.seed(138L)
x <- rmb(1e5, 150, 250)
testx <- mommb(x)
expect_true(abs(testx$g / 150 - 1) <= 0.02)
expect_true(abs(testx$b / 250 - 1) <= 0.05)

# Test Functionality: m = TRUE
z <- c(mean(x), mean(x ^ 2))
testz <- mommb(z, m = TRUE)
expect_equal(testz$g, testx$g, tolerance = tol)
expect_equal(testz$b, testx$b, tolerance = tol)

# Testing simple findb branch
expect_identical(mommb(c(1, 1, 1))$b, 0)                    # findb => 0
x <- c(0.9, 0.9, 0.9, 0.93785026012074624)
expect_silent(mommb(x, tol = 0.002))                        # findb => 1 / g
# Malformed input created solely to test branches in findb
expect_error(mommb(c(0.2, 0.6, (sqrt(2.6) + 1) / 2)))       # findb => Inf
expect_error(mommb(c(0.2, 0.6, 0.1770747858956514)))        # findb => 1

# Test Error Trapping
set.seed(76L)
expect_error(mommb(rmb(10, 10, 9)), "insufficient data")
expect_error(mommb(rmb(10, 10, 9), tol = 1e-16), "looser tolerance")
expect_error(mommb(NA_real_, na.rm = FALSE), "passed as FALSE")
expect_error(mommb(rmb(10, 10, 9), m = TRUE), "first and second raw moments")

# Test trace
expect_message(mommb(x, trace = TRUE), "i: 1")
