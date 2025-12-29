// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

extern SEXP pmb_c(SEXP q, SEXP g, SEXP b, SEXP lower_tail, SEXP log_p) {
  const R_xlen_t n = xlength(q);
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *pq = REAL(q);
  double *pg = REAL(g);
  double *pb = REAL(b);
  Rboolean lt = asLogical(lower_tail);
  Rboolean lp = asLogical(log_p);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  Memzero(pret, n);

  for (R_xlen_t i = 0; i < n; ++i) {
    double gi = pg[i % gg];
    double bi = pb[i % bb];
    double gm1 = gi - 1.0;
    double gb = gi * bi;
    if (ISNA(pq[i]) || ISNA(gi) || ISNA(bi)) {
      pret[i] = NA_REAL;
    } else if (gi < 1.0 || bi < 0.0 || ISNAN(pq[i] + gi + bi)) {
      pret[i] = R_NaN;
    } else if (pq[i] >= 1.0) {
      pret[i] = 1.0;
    } else if (gi == 1.0 || bi == 0.0 || pq[i] < 0.0) {
      pret[i] = 0.0;
    } else if (bi == 1.0) {
      pret[i] = 1.0 - 1.0 / (1.0 + gm1 * pq[i]);
    } else if (gb == 1.0) {
      pret[i] = 1.0 - R_pow(bi, pq[i]);
    } else {
      pret[i] = 1.0 - (1.0 - bi) / (gm1 * R_pow(bi, 1.0 - pq[i]) + 1.0 - gb);
    }

    pret[i] = lt ? pret[i] : 0.5 - pret[i] + 0.5; // See dpq.h
    pret[i] = lp ? log(pret[i]) : pret[i];
  }

    UNPROTECT(1);
    return(ret);
}
