// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "dpqrmb.h"

extern SEXP rmb_c(SEXP n_, SEXP g, SEXP b) {
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  const R_xlen_t n = asReal(n_);
  double *pg = REAL(g);
  double *pb = REAL(b);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  memset(pret, 0, n * sizeof(double));

  GetRNGstate();
  for (R_xlen_t i = 0; i < n; ++i) {
    double x = unif_rand();
    R_xlen_t ig = i % gg;
    R_xlen_t ib = i % bb;
    double gm1 = pg[ig] - 1.0;
    double gb = pg[ig] * pb[ib];
    double xc = 0.5 - x + 0.5; // x-complement
    if (ISNA(pg[ig]) || ISNA(pb[ib])) {
      pret[i] = NA_REAL;
    } else if (pg[ig] < 1.0 || pb[ib] < 0.0 || ISNAN(pg[ig] + pb[ib])) {
      pret[i] = R_NaN;
    } else if (x >= 1.0 - 1.0 / pg[ig]) {
      pret[i] = 1.0;
    } else if (pg[ig] == 1.0 || pb[ib] == 0.0 || x == 0.0) {
      pret[i] = 0.0;
    } else if (pb[ib] == 1.0) {
      pret[i] = (1.0 - 1.0 / xc) / gm1 ;
    } else if (gb == 1.0) {
      pret[i] = log(xc) / log(pb[ib]);
    } else {
      pret[i] = 1.0 -
        (log((1.0 - pb[ib]) / xc - 1.0 + gb) - log(gm1)) / log(pb[ib]);
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return(ret);
}
