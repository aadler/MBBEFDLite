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

testx1b <-  mommb(x1, opts = list(alg = "LS"))
expect_true(abs(testx1b$g / 150 - 1) <= 0.02)
expect_true(abs(testx1b$b / 250 - 1) <= 0.05)
expect_equal(testx1b$sqerr,
             (log(testx1b$g * testx1b$b) * (1 - testx1b$b) /
                (log(testx1b$b) * (1 - testx1b$g * testx1b$b)) - mean(x1)) ^ 2,
             tolerance = tol)

# Test second-iteration Bernegger
set.seed(14)
x <- rmb(10, 6, 0.2)
fit <- mommb(x, opts = list(alg = "LS"))
expect_identical(fit$iter, 2)

# Test Functionality: m = TRUE
z <- c(mean(x1), var(x1))
testz <- mommb(z, m = TRUE)
expect_equal(testz$g, testx1$g, tolerance = tol)
expect_equal(testz$b, testx1$b, tolerance = tol)

# Testing simple findb branch
expect_identical(mommb(c(1, 1, 1))$b, 0)                    # findb => 0

expect_equal(MBBEFDLite:::findb(log(5) / 4, 5, 1e3), 1,
             tolerance = tol)                        # g = 5; findb => 1

expect_error(suppressWarnings(mommb(c(0, 0.5, 1))))         # findb => Inf
expect_warning(mommb(c(0, 0.5, 1)), "Parameter b is Inf")

# Malformed input created solely to test branches in findb. Values found using
# ChatGPT

a <- 2 * log(2)
d <- sqrt(2 - a ^ 2)
expect_error(mommb(c(a - d, a, a + d)))    # mean must be in [0, 1]

# Test mommb Error Trapping
x <- c(0.9, 0.9, 0.9, 0.93785026012074624)
expect_error(mommb(x, opts = list(maxit = "3L")),
             "maxit must be a positive integer.")
expect_error(mommb(x, opts = list(maxb = "3L")), "maxb must be positive.")
expect_error(mommb(x, opts = list(maxit = 3L)), "Try increasing maxit or tol")
expect_error(mommb(c(0.2, 0.3), m = TRUE), "less than or equal")
expect_error(mommb(x, opts = list(maxb = 2e-16)), "Parameter b approaching")

set.seed(76L)
expect_error(mommb(rmb(10, 10, 9)), "insufficient data")
set.seed(76L)
expect_error(mommb(rmb(10, 10, 9), opts = list(alg = "LS")),
             "insufficient data")
expect_error(mommb(rmb(10, 10, 9), tol = 1e-16), "looser tolerance")
expect_error(mommb(NA_real_, na.rm = FALSE), "passed as FALSE")
expect_error(mommb(rmb(10, 10, 9), m = TRUE), "other than 2")
expect_error(mommb(x, opts = list(alg = "sqrt")), "Algorithm must be")

# Test trace
expect_message(mommb(x, opts = list(trace = TRUE)), "i: 1")
expect_error(mommb(x, opts = list(trace = "EM")), "trace must be")
expect_message(mommb(x, opts = list(alg = "LS", trace = TRUE)),
               "trace is ignored")

# Test mugb
expect_true(is.nan(suppressWarnings(MBBEFDLite:::mugb(4, -1))))
expect_true(is.nan(MBBEFDLite:::mugb(0.2, 1)))
expect_equal(MBBEFDLite:::mugb(10, 0), 1, tolerance = tol)
expect_equal(MBBEFDLite:::mugb(10, 1 + 1e-36), log(10) / 9, tolerance = tol)
expect_equal(MBBEFDLite:::mugb(5, 0.2), -0.8 / log(0.2), tolerance = tol)

# Test cvgb
cvp <- function(x) sqrt(2 * x - 1)
expect_equal(MBBEFDLite:::cvgb(10, 0), cvp(0.5), tolerance = tol)
expect_equal(MBBEFDLite:::cvgb(3, 1 + 1e-36),
             cvp((2 - log(3)) / log(3) ^ 2), tolerance = tol)
expect_equal(MBBEFDLite:::cvgb(5, 0.2),
             cvp((log(0.2) * 0.2 + 0.8) / 0.64), tolerance = tol)

