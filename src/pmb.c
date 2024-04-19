// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+


#include <Rmath.h>

#include "pmb.h"

extern SEXP pmb_c(SEXP q, SEXP g, SEXP b, SEXP lt, SEXP lg) {
  const R_xlen_t n = xlength(q);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *pq = REAL(q);
  double *pg = REAL(g);
  double *pb = REAL(b);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  memset(pret, 0, n * sizeof(double));

  for (R_xlen_t i = 0; i < n; ++i) {
    R_xlen_t ig = i % gg;
    R_xlen_t ib = i % bb;
    double gm1 = pg[ig] - 1.0;
    double gb = pg[ig] * pb[ib];
    if (ISNA(pq[i])) {
      pret[i] = NA_REAL;
    } else if (pg[ig] < 1.0 | pb[ib] < 0.0 | ISNAN(pq[i])) {
      pret[i] = R_NaN;
    } else if (pq[i] >= 1.0) {
      pret[i] = 1.0;
    } else if (pg[ig] == 1.0 | pb[ib] == 0.0 | pq[i] < 0.0) {
      pret[i] = 0.0;
    } else if (pb[ib] == 1.0) {
      pret[i] = 1.0 - 1.0 / (1.0 + gm1 * pq[i]);
    } else if (gb == 1.0) {
      pret[i] = 1.0 - R_pow(pb[ib], pq[i]);
    } else {
      double gm1b1x = gm1 * R_pow(pb[ib], 1 - pq[i]);
      pret[i] = 1.0 - (1.0 - pb[ib]) / (gm1b1x + 1.0 - gb);
    }
  }

  if (asLogical(lt) == 0) {
    for (R_xlen_t i = 0; i < n; ++i) {
      pret[i] = 1.0 - (pret[i]);
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
