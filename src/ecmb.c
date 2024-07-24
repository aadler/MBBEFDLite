// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

extern SEXP ecmb_c(SEXP x, SEXP g, SEXP b, SEXP lower_tail) {
  const R_xlen_t n = xlength(x);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *px = REAL(x);
  double *pg = REAL(g);
  double *pb = REAL(b);
  Rboolean lt = asLogical(lower_tail);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  memset(pret, 0, n * sizeof(double));

  for (R_xlen_t i = 0; i < n; ++i) {
    double gi = pg[i % gg];
    double bi = pb[i % bb];
    double gm1 = gi - 1.0;
    double gb = bi * gi;
    if (ISNA(px[i]) || ISNA(gi) || ISNA(bi)) {
      pret[i] = NA_REAL;
    } else if (gi < 1.0 || bi < 0.0 || ISNAN(px[i] + gi + bi)) {
      pret[i] = R_NaN;
    } else if (px[i] <= 0.0) {
      pret[i] = 0.0;
    } else if (px[i] >= 1.0) {
      pret[i] = 1.0;
    } else if (gi == 1.0 || bi == 0.0) {
      pret[i] = px[i];
    } else if (bi == 1.0) {
      pret[i] = log1p(gm1 * px[i]) / log(gi);
    } else if (gb == 1.0) {
      pret[i] = (1.0 - R_pow(bi, px[i])) / (1.0 - bi);
    } else {
      pret[i] = log((gm1 * bi + (1.0 - gb) * R_pow(bi, px[i])) / (1.0 - bi)) /
        log(gb);
    }
    pret[i] = lt ? pret[i] : 0.5 - pret[i] + 0.5; // See dpq.h
  }

    UNPROTECT(1);
    return(ret);
}
