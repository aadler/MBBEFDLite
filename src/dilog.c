// Copyright Avraham Adler (c) 2026
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

// Developed with help of ChatGPT. Apparently similar to Cephes and Boost

static double dilog_s(double x) {
  // uses the fact that the ratio of two successive terms is x * ((k - 1) / k)²
  const int Maxk = 1000;
  const double EPS = 2.2204460492503131e-16;
  double sum = x;
  double term = x;
  int k = 2;
  bool CONVERGED = false;

  while (!CONVERGED) {
    const double tk = (k - 1.0) / k;
    term *= x * tk * tk;
    sum += term;
    CONVERGED = (fabs(term) < EPS * fabs(sum) || k > Maxk);
    k++;
  }

  return(sum);
}

SEXP dilog_c(SEXP x_) {
  double x = REAL(x_)[0];
  const double EPS = 2.2204460492503131e-16;
  const double PISQ_6 = M_PI * M_PI / 6.0;
  SEXP ret = PROTECT(allocVector(REALSXP, 1));

  if (fabs(x - 1.0) < EPS) {
    REAL(ret)[0] = PISQ_6;
  } else if (fabs(x) < EPS) {
    REAL(ret)[0] = 0.0;
  } else if (x < -1) {
    double lx = log(-x);
    REAL(ret)[0] = -dilog_s(1.0 / x) - PISQ_6 - 0.5 * lx * lx;
  } else if (x <= 0.5) {
    REAL(ret)[0] = dilog_s(x);
  } else if (x < 1.0) {
    REAL(ret)[0] = PISQ_6 - log(x) * log1p(-x) - dilog_s(1.0 - x);
  } else if (x > 1.0) {
    double lx = log(x);
    REAL(ret)[0] = PISQ_6 - 0.5 * lx * lx - dilog_s(1.0 / x);
  } else {
    REAL(ret)[0] = R_NaN;
  }

  UNPROTECT(1);
  return(ret);
}
