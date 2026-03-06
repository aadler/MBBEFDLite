# Copyright Avraham Adler (c) 2026
# SPDX-License-Identifier: MPL-2.0+

tol <- 20 * .Machine$double.eps

# Eight closed form values of Li₂(x)
pisq_6 <- pi ^ 2 / 6
phi <- (1 + sqrt(5)) / 2
lphi <- log(phi) ^ 2
expect_equal(MBBEFDLite:::dilog(0), 0, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(1), pisq_6, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(-1), -0.5 * pisq_6, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(0.5), 0.5 * (pisq_6 - log(2) ^ 2),
             tolerance = tol)
expect_equal(MBBEFDLite:::dilog(0.5 * (3 - sqrt(5))),
             pi ^ 2 / 15 - lphi, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(0.5 * (-1 + sqrt(5))),
             pi ^ 2 / 10 - lphi, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(0.5 * (1 - sqrt(5))),
             -pi ^ 2 / 15 + 0.5 * lphi, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(-phi), -pi ^ 2 / 10 - lphi, tolerance = tol)

# From Wolfram Mathworld
A <- MBBEFDLite:::dilog(sqrt(2) - 1) - MBBEFDLite:::dilog(1 - sqrt(2))
B <- pi ^ 2 / 8 - 0.5 * log(1 + sqrt(2)) ^ 2
expect_equal(A, B, tolerance = tol)

# Unclosed forms using Wolfram Alpha
expect_equal(MBBEFDLite:::dilog(0.75), 0.978469392930306, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(2), 2.46740110027234, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(-2), -1.43674636688368, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(-200), -15.6760237613305, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(200), -10.7512215885639, tolerance = tol)

# Edge cases using Wolfram Alpha
expect_equal(MBBEFDLite:::dilog(1.000000001), 1.64493408857149, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(0.999999999), 1.64493404512496, tolerance = tol)
expect_equal(MBBEFDLite:::dilog(-0.999999999), -0.822467032730966,
             tolerance = tol)
expect_equal(MBBEFDLite:::dilog(-1.000000001), -0.82246703411726,
             tolerance = tol)

# Inversion identity: Li₂(x) + Li₂(1/x) = -π²/6 - 0.5*log(-x)^2 for x < 0
x <- -3.7
expect_equal(MBBEFDLite:::dilog(x) + MBBEFDLite:::dilog(1 / x),
             -pisq_6 - 0.5 * log(-x)^2, tolerance = tol)

# Euler reflection: Li₂(x) + Li₂(1-x) = π²/6 - log(x)*log(1-x)
x <- 0.3
expect_equal(MBBEFDLite:::dilog(x) + MBBEFDLite:::dilog(1 - x),
             pisq_6 - log(x) * log(1 - x), tolerance = tol)

# Li₂(x) + Li₂(x/(x-1)) = -0.5 * log(1-x)^2
x <- -0.6
expect_equal(MBBEFDLite:::dilog(x) + MBBEFDLite:::dilog(x / (x - 1)),
             -0.5 * log1p(-x)^2, tolerance = tol)

# Vectorized
expect_equal(MBBEFDLite:::dilog(c(-10, -1, -0.8, -0.2, 0, 0.2, 0.8, 1, 10)),
             c(MBBEFDLite:::dilog(-10),
               MBBEFDLite:::dilog(-1),
               MBBEFDLite:::dilog(-0.8),
               MBBEFDLite:::dilog(-0.2),
               MBBEFDLite:::dilog(0),
               MBBEFDLite:::dilog(0.2),
               MBBEFDLite:::dilog(0.8),
               MBBEFDLite:::dilog(1),
               MBBEFDLite:::dilog(10)),
             tolerance = tol)
