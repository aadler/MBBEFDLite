// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

SEXP dmb_c(SEXP x, SEXP g, SEXP b, SEXP give_log) {
  const R_xlen_t n = xlength(x);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *px = REAL(x);
  double *pg = REAL(g);
  double *pb = REAL(b);
  Rboolean gl = asLogical(give_log);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);

  R_xlen_t ig = 0;
  R_xlen_t ib = 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    double gi = pg[ig];
    double bi = pb[ib];
    if (++ig == gg) ig = 0;
    if (++ib == bb) ib = 0;

    if (ISNA(px[i]) || ISNA(gi) || ISNA(bi)) {
      pret[i] = NA_REAL;
      continue;
    }

    else if (gi < 1.0 || bi < 0.0 || !R_FINITE(px[i]) || !R_FINITE(gi) ||
             !R_FINITE(bi)) {
      pret[i] = R_NaN;
      continue;
    }

    if (gi == 1.0 || bi == 0.0 || px[i] < 0.0 || px[i] >= 1.0) {
      pret[i] = 0.0;
      continue;
    }

    double gm1 = gi - 1.0;

    if (bi == 1.0) {
      double gm1px = gm1 * px[i] + 1.0;
      pret[i] = gm1 / (gm1px * gm1px);
      continue;
    }

    double gb = gi * bi;
    double lb = log(bi);
    double bix = exp(px[i] * lb);

    if (gb == 1.0) {
      pret[i] = -lb * bix;
    } else {
      double gm1b1x = gm1 * bi / bix;
      double gm1b1x1mgb = gm1b1x + (1.0 - gb);
      pret[i] = (bi - 1.0) * gm1b1x * lb / (gm1b1x1mgb * gm1b1x1mgb);
    }
  }

  if (gl) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (pret[i] > 0.0) {
        // log positive values
        pret[i] = log(pret[i]);
      } else if (pret[i] == 0.0) {
        // convert 0 to NegInf
        pret[i] = R_NegInf;
      }
        // Allow NA and NaN to flow through
    }
  }

  UNPROTECT(1);
  return(ret);
}
