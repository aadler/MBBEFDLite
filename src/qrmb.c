// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

static double quantilemb (double p, double g, double b) {

  if (ISNA(p) || ISNA(g) || ISNA(b)) {
    return(NA_REAL);
  }

  if (!R_FINITE(p) || !R_FINITE(g) || !R_FINITE(b)) {
    return(R_NaN);
  }

  if (g < 1.0 || b < 0.0 || p < 0.0 || p > 1.0) {
    return(R_NaN);
  }

  if (g == 1.0 || b == 0.0 || p == 0.0) {
    return(0.0);
  }

  if (p >= 1.0 - 1.0 / g) {
    return(1.0);
  }

  double gm1 = g - 1.0;
  double pc = 0.5 - p + 0.5; // p-complement; avoid cancellation

  if (b == 1.0) {
    return(p / (pc * gm1));
  }

  double gb = g * b;
  double lb = log(b);

  if (gb == 1.0) {
    return(log(pc) / lb);
  } else {
    return(1.0 - log(((1.0 - b) / pc - 1.0 + gb) / gm1) / lb);
  }
}

SEXP qmb_c(SEXP p, SEXP g, SEXP b, SEXP lower_tail, SEXP log_p) {
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

  R_xlen_t ig = 0;
  R_xlen_t ib = 0;
  double x;

  for (R_xlen_t i = 0; i < n; ++i) {

    if (lp) {
      x = lt ? exp(pp[i]) : -expm1(pp[i]); // -expm1(x) = 1 - exp(x)
    } else {
      x = lt ? pp[i] : 0.5 - pp[i] + 0.5; // Use for precision (see R's dpq.h)
    }

    pret[i] = quantilemb(x, pg[ig], pb[ib]);
    if (++ig == gg) ig = 0;
    if (++ib == bb) ib = 0;
  }

  UNPROTECT(1);
  return(ret);
}

SEXP rmb_c(SEXP n_, SEXP g, SEXP b) {
  const R_xlen_t gg = xlength(g);
  const R_xlen_t bb = xlength(b);
  double *pg = REAL(g);
  double *pb = REAL(b);

  // Convert long integer to R_xlen_t via cast to and from REAL
  double dn = asReal(n_);
  if (!R_FINITE(dn) || dn < 0 || dn > R_XLEN_T_MAX) error("invalid 'n'");
  R_xlen_t n = (R_xlen_t) dn;

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);

  R_xlen_t ig = 0;
  R_xlen_t ib = 0;

  GetRNGstate();

  for (R_xlen_t i = 0; i < n; ++i) {
    pret[i] = quantilemb(unif_rand(), pg[ig], pb[ib]);
    if (++ig == gg) ig = 0;
    if (++ib == bb) ib = 0;
  }

  PutRNGstate();

  UNPROTECT(1);
  return(ret);
}
