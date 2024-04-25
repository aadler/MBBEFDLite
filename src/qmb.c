// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "dpqrmb.h"

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
  memset(pret, 0, n * sizeof(double));

  for (R_xlen_t i = 0; i < n; ++i) {
    R_xlen_t ig = i % gg;
    R_xlen_t ib = i % bb;
    double gm1 = pg[ig] - 1.0;
    double gb = pg[ig] * pb[ib];
    pp[i] = lp ? exp(pp[i]) : pp[i];
    pp[i] = lt ? pp[i] : 0.5 - pp[i] + 0.5; // See dpq.h
    double pc = 0.5 - pp[i] + 0.5; // p-complement
    if (ISNA(pp[i]) || ISNA(pg[ig]) || ISNA(pb[ib])) {
      pret[i] = NA_REAL;
    } else if (pg[ig] < 1.0 || pb[ib] < 0.0 || ISNAN(pp[i] + pg[ig] + pb[ib]) ||
      pp[i] < 0 || pp[i] > 1.0) {
      pret[i] = R_NaN;
    } else if (pp[i] >= 1.0 - 1.0 / pg[ig]) {
      pret[i] = 1.0;
    } else if (pg[ig] == 1.0 | pb[ib] == 0.0 | pp[i] == 0.0) {
      pret[i] = 0.0;
    } else if (pb[ib] == 1.0) {
      pret[i] = (1.0 - 1.0 / pc) / gm1 ;
    } else if (gb == 1.0) {
      pret[i] = log(pc) / log(pb[ib]);
    } else {
      pret[i] = 1.0 -
        (log((1.0 - pb[ib]) / pc - 1.0 + gb) - log(gm1)) / log(pb[ib]);
    }
  }

  UNPROTECT(1);
  return(ret);
}
