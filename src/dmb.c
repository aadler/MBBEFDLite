// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+


#include <Rmath.h>
// #include <Rconfig.h>

#include "dmb.h"

extern SEXP dmb_c(SEXP x, SEXP g, SEXP b, SEXP lg) {
  const R_xlen_t n = xlength(x);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *px = REAL(x);
  double *pg = REAL(g);
  double *pb = REAL(b);
  Rboolean pl = asLogical(lg);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  memset(pret, 0, n * sizeof(double));

  for (R_xlen_t i = 0; i < n; ++i) {
    R_xlen_t igg = i % gg;
    R_xlen_t ibb = i % bb;
    if (pg[igg] < 1.0 | pb[ibb] < 0.0 | ISNAN(px[i])) {
      pret[i] = R_NaN;
    } else if (pg[igg] == 1.0 | pb[ibb] == 0.0 | px[i] < 0.0 | px[i] >= 1.0) {
      pret[i] = 0.0;
    } else if (pb[ibb] == 1.0) {
      pret[i] = (pg[igg] - 1) / pow1p((pg[igg] - 1) * px[i], 2);
    } else if (pb[ibb] * pg[igg] == 1.0) {
      pret[i] = -log(pb[ibb]) * R_pow(pb[ibb], px[i]);
    } else {
      double b1x = R_pow(pb[ibb], 1 - px[i]);
      pret[i] = (pb[ibb] - 1) * (pg[igg] - 1) * b1x * log(pb[ibb]) /
        R_pow((pg[igg] - 1) * b1x + 1 - pb[ibb] * pg[igg], 2);
    }
  }

  if (pl == 1) {
    for (R_xlen_t i = 0; i < n; ++i) {
      pret[i] = log(pret[i]);
    }
  }

    UNPROTECT(1);
    return(ret);
}
