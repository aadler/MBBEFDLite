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
expect_identical(mommb(c(1, 1, 1))$b, 0)                     # findb returns 0
# Below is invalid input but findb returns Inf
expect_error(mommb(c(0.6798293601721539, 1.1838579062372445), tol = 1e-5))

# Below is invalid input but findb returns 1 / g
expect_error(mommb(c(1.0032341003417968, 1.004322052001954), tol = 1e-5))
expect_error(mommb(c(0, 0.31824773751438351), tol = 1e-5))   # findb returns 1

# Test Error Trapping
set.seed(76L)
expect_error(mommb(rmb(10, 10, 9)), "insufficient data")
expect_error(mommb(rmb(10, 10, 9), tol = 1e-16), "looser tolerance")
expect_error(mommb(NA_real_, na.rm = FALSE), "passed as FALSE")

