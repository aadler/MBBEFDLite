// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+


#include <Rmath.h>

#include "dmb.h"

extern SEXP dmb_c(SEXP x, SEXP g, SEXP b, SEXP lg) {
  const R_xlen_t n = xlength(x);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *px = REAL(x);
  double *pg = REAL(g);
  double *pb = REAL(b);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  memset(pret, 0, n * sizeof(double));

  for (R_xlen_t i = 0; i < n; ++i) {
    R_xlen_t ig = i % gg;
    R_xlen_t ib = i % bb;
    double gm1 = pg[ig] - 1.0;
    if (pg[ig] < 1.0 | pb[ib] < 0.0 | ISNAN(px[i])) {
      pret[i] = R_NaN;
    } else if (pg[ig] == 1.0 | pb[ib] == 0.0 | px[i] < 0.0 | px[i] >= 1.0) {
      pret[i] = 0.0;
    } else if (pb[ib] == 1.0) {
      pret[i] = gm1 / R_pow_di(gm1 * px[i] + 1.0, 2);
    } else if (pb[ib] * pg[ig] == 1.0) {
      pret[i] = -log(pb[ib]) * R_pow(pb[ib], px[i]);
    } else {
      double gm1b1x = gm1 * R_pow(pb[ib], 1 - px[i]);
      pret[i] = (pb[ib] - 1) * gm1b1x * log(pb[ib]) /
        R_pow_di(gm1b1x + (1 - pb[ib] * pg[ig]), 2);
    }
  }

  if (asLogical(lg) == 1) {
    for (R_xlen_t i = 0; i < n; ++i) {
      pret[i] = log(pret[i]);
    }
  }

    UNPROTECT(1);
    return(ret);
}
