# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps

set.seed(9712L)
u <- runif(100L)
## Standard b & g
g <- 20
b <- 0.5
control <- qmb(u, g, b)
set.seed(9712L)
expect_identical(rmb(100L, g, b), control)

# Vectored parameters and input (other edge cases handed in test_qmb.R).
g <- c(1.2, 4, 100)
b <- c(0.001, 0.17)
vv <- 4:9
set.seed(9712L)
control <- c(rmb(1L, g[1L], b[1L]),
             rmb(1L, g[2L], b[2L]),
             rmb(1L, g[3L], b[1L]),
             rmb(1L, g[1L], b[2L]),
             rmb(1L, g[2L], b[1L]),
             rmb(1L, g[3L], b[2L]))

set.seed(9712L)
expect_identical(rmb(vv, g, b), control)

# Test c
## Scalar
set.seed(9712L)
control <- rmb(5L, MBBEFDLite:::c2gb(4)$g, MBBEFDLite:::c2gb(4)$b)
set.seed(9712L)
expect_identical(rmb(5L, c = 4), control)
## Scalar
c <- c(2, 3)
set.seed(9712L)
control <- rmb(5L, MBBEFDLite:::c2gb(c)$g, MBBEFDLite:::c2gb(c)$b)
set.seed(9712L)
expect_identical(rmb(5L, c = c), control)
