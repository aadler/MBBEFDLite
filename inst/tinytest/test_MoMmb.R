# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps

# Test Functionality
set.seed(138L)
x <- rmb(1e5, 150, 250)
testx <- mommb(x)
expect_true(abs(testx$g / 150 - 1) <= 0.02)
expect_true(abs(testx$b / 250 - 1) <= 0.05)

# Testing simple findb branch
expect_identical(mommb(c(1, 1, 1))$b, 0)                 # findb returns 0

# Test Error Trapping
set.seed(76L)
expect_error(mommb(rmb(10, 10, 9)), "insufficient data")
expect_error(mommb(rmb(10, 10, 9), tol = 1e-16), "looser tolerance")
expect_error(mommb(NA_real_, na.rm = FALSE), "passed as FALSE")
