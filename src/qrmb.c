// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

double quantilemb (double p, double g, double b) {
  double gm1 = g - 1.0;
  double gb = g * b;
  double pc = 0.5 - p + 0.5; // p-complement
  if (ISNA(p) || ISNA(g) || ISNA(b)) {
    return(NA_REAL);
  } else if (g < 1.0 || b < 0.0 || ISNAN(p + g + b) || p < 0.0 || p > 1.0) {
    return(R_NaN);
  } else if (g == 1.0 || b == 0.0 || p == 0.0) {
    return(0.0);
  } else if (p >= 1.0 - 1.0 / g) {
    return(1.0);
  } else if (b == 1.0) {
    return(p / (pc * gm1));
  } else if (gb == 1.0) {
    return(log(pc) / log(b));
  } else {
    return(1.0 - (log((1.0 - b) / pc - 1.0 + gb) - log(gm1)) / log(b));
  }
}

extern SEXP qmb_c(SEXP p, SEXP g, SEXP b, SEXP lower_tail, SEXP log_p) {
  const R_xlen_t n = xlength(p);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *pp = REAL(p);
  double *pg = REAL(g);
  double *pb = REAL(b);
  Rboolean lt = asLogical(lower_tail);
  Rboolean lp = asLogical(log_p);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  Memzero(pret, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    double x = lp ? exp(pp[i]) : pp[i];
    x = lt ? x : 0.5 - x + 0.5; // See dpq.h
    pret[i] = quantilemb(x, pg[i % gg], pb[i % bb]);
  }

  UNPROTECT(1);
  return(ret);
}

extern SEXP rmb_c(SEXP n_, SEXP g, SEXP b) {
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  const R_xlen_t n = asReal(n_);
  double *pg = REAL(g);
  double *pb = REAL(b);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  Memzero(pret, n);

  GetRNGstate();
  for (R_xlen_t i = 0; i < n; ++i) {
    pret[i] = quantilemb(unif_rand(), pg[i % gg], pb[i % bb]);
  }
  PutRNGstate();

  UNPROTECT(1);
  return(ret);
}
