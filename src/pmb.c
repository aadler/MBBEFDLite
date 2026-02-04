// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

SEXP pmb_c(SEXP q, SEXP g, SEXP b, SEXP lower_tail, SEXP log_p) {
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

  R_xlen_t ig = 0;
  R_xlen_t ib = 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    double gi = pg[ig];
    double bi = pb[ib];

    if (++ig == gg) ig = 0;
    if (++ib == bb) ib = 0;

    if (ISNA(pq[i]) || ISNA(gi) || ISNA(bi)) {
      pret[i] = NA_REAL;
      continue;
    }

    if (gi < 1.0 || bi < 0.0 || !R_FINITE(pq[i]) || !R_FINITE(gi) ||
        !R_FINITE(bi)) {
      pret[i] = R_NaN;
      continue;
    }

    if (pq[i] >= 1.0) {
      // pret[i] = 1.0;
      pret[i] = 1.0;
      continue;
    }

    if (gi == 1.0 || bi == 0.0 || pq[i] < 0.0) {
      // pret[i] = 0.0;
      pret[i] = 0.0;
      continue;
    }

    double gm1 = gi - 1.0;
    double tmp;

    if (bi == 1.0) {
      tmp = 1.0 / (1.0 + gm1 * pq[i]);
      pret[i] = 0.5 - tmp + 0.5;  // 0.5 - x + 0.5 is more accurate than 1.0 - x
      continue;
    }

    double gb = gi * bi;
    double lbi = log(bi);
    double biq = exp(pq[i] * lbi);

    if (gb == 1.0) {
      pret[i] = 0.5 - biq + 0.5;
    } else {
      tmp = (1.0 - bi) / (gm1 * bi / biq + 1.0 - gb);
      pret[i] = 0.5 - tmp + 0.5;
    }
  }

  if (lt) { // Lower tail section
    if (lp) {
      for (R_xlen_t i = 0; i < n; ++i) {
        if (pret[i] > 0.0) { // log positive values
          pret[i] = log(pret[i]);
        } else if (pret[i] == 0.0) {  // log(0) is NegInf
          pret[i] = R_NegInf;
        }
        // Allow NA and NaN to flow through
      }
    } // Lower tail non-log is default so no else here
  } else { // Upper tail section
    if (lp) { // Upper tail Log
      for (R_xlen_t i = 0; i < n; ++i) {
        if (pret[i] < 1.0) { // log complementary positive values
          pret[i] = log1p(-pret[i]);
        } else if (pret[i] == 1.0) { // convert log(1 - 1) to NegInf
          pret[i] = R_NegInf;
        }
      }
    } else { // Upper tail Normal
      for (R_xlen_t i = 0; i < n; ++i) {
        if (R_FINITE(pret[i])) {
          pret[i] = 0.5 - pret[i] + 0.5;
        }
      }
    }
  }

  UNPROTECT(1);
  return(ret);
}
