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

    if (gi < 1.0 || bi < 0.0 || !R_finite(px[i]) || !R_finite(gi) ||
        !R_finite(bi)) {
      pret[i] = R_NaN;
      continue;
    }

    if (px[i] <= 0.0) {
      pret[i] = 0.0;
      continue;
    }

    if (px[i] >= 1.0) {
      pret[i] = 1.0;
      continue;
    }

    if (gi == 1.0 || bi == 0.0) {
      pret[i] = px[i];
      continue;
    }

    double gm1 = gi - 1.0;

    if (bi == 1.0) {
      pret[i] = log1p(gm1 * px[i]) / log1p(gm1);
      continue;
    }

    double ombi = 1.0 - bi;
    double gb = bi * gi;
    double lbi = log(bi);
    double bix = exp(px[i] * lbi);

    if (gb == 1.0) {
      pret[i] = (1.0 - bix) / ombi;
    } else {
      pret[i] = log((gm1 * bi + (1.0 - gb) * bix) / ombi) / log(gb);
    }
  }

  if (!lt) {
    for (R_xlen_t i = 0; i < n; ++i) {
      pret[i] =  0.5 - pret[i] + 0.5; // Avoid cancellation; see dpq.h
    }
  }

  UNPROTECT(1);
  return(ret);
}
