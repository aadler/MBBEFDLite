# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

tol <- 10 * .Machine$double.eps
x <- c(NaN, NA_real_, seq(-0.05, 1.05, 0.05))

# Standard b & g and passing log.p and lower.tail
g <- 20
b <- 0.5

control <- c(NaN, NA, 0, 0, 0.40120963545198918, 0.57693371329906551,
             0.67551641239099069, 0.73858045887054991, 0.78236907383302978,
             0.81453124867180904,  0.83914170955380674, 0.85857029312602218,
             0.8742891819815175, 0.88726116159525426, 0.89814240493843012,
             0.90739552971536974, 0.91535606623643084, 0.92227329508477796,
             0.92833628856856798,  0.93369105420497067, 0.93845213614173995,
             0.94271065786700514,  0.94654001782806174, 1, 1)

expect_equal(pmb(x, g, b), control, tolerance = tol)
expect_equal(pmb(x, g, b, lower.tail = FALSE), 1 - control, tolerance = tol)
expect_equal(pmb(x, g, b, log.p = TRUE), log(control), tolerance = tol)
expect_equal(pmb(x, g, b, lower.tail = FALSE, log.p = TRUE), log(1 - control),
             tolerance = tol)

# Nonstandard g & b
## g < 1 and b < 0
expect_true(is.nan(pmb(0.5, 0.2, 6)))
expect_true(is.nan(pmb(0.5, 1.2, -0.3)))

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
control <- c(pmb(x[1L], g[1L], b[1L]),
             pmb(x[2L], g[2L], b[2L]),
             pmb(x[3L], g[3L], b[1L]),
             pmb(x[4L], g[1L], b[2L]),
             pmb(x[5L], g[2L], b[1L]),
             pmb(x[6L], g[3L], b[2L]))

expect_identical(pmb(x, g, b)[1:6], control)
