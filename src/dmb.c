// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

extern SEXP dmb_c(SEXP x, SEXP g, SEXP b, SEXP give_log) {
  const R_xlen_t n = xlength(x);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *px = REAL(x);
  double *pg = REAL(g);
  double *pb = REAL(b);
  Rboolean gl = asLogical(give_log);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  Memzero(pret, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    double gi = pg[i % gg];
    double bi = pb[i % bb];
    double gm1 = gi - 1.0;
    double gb = gi * bi;
    if (ISNA(px[i]) || ISNA(gi) || ISNA(bi)) {
      pret[i] = NA_REAL;
    } else if (gi < 1.0 || bi < 0.0 || ISNAN(px[i] + gi + bi)) {
      pret[i] = R_NaN;
    } else if (gi == 1.0 || bi == 0.0 || px[i] < 0.0 || px[i] >= 1.0) {
      pret[i] = 0.0;
    } else if (bi == 1.0) {
      pret[i] = gm1 / R_pow_di(gm1 * px[i] + 1.0, 2);
    } else if (gb == 1.0) {
      pret[i] = -log(bi) * R_pow(bi, px[i]);
    } else {
      double gm1b1x = gm1 * R_pow(bi, 1 - px[i]);
      pret[i] = (bi - 1.0) * gm1b1x * log(bi) /
        R_pow_di(gm1b1x + (1.0 - gb), 2);
    }

    pret[i] = gl ? log(pret[i]) : pret[i];
  }

    UNPROTECT(1);
    return(ret);
}
