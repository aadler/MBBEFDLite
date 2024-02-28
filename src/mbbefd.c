// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rconfig.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP dmb_c(SEXP x, SEXP g, SEXP b, SEXP lg) {
  const R_xlen_t n = xlength(x);
  double *px = REAL(x);
  double *pg = REAL(g);
  double *pb = REAL(b);
  Rboolean pl = asLogical(lg);

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);
  memset(pret, 0, n * sizeof(double));

  if (pg[0] < 1.0 | pb[0] < 0.0) {
    for (R_xlen_t i = 0; i < n; ++i) {
      pret[i] = R_NaN;
    }
    UNPROTECT(1);
    return(ret);
  }

  if (pg[0] == 1.0 | pb[0] == 0.0) {
    for (R_xlen_t i = 0; i < n; ++i) {
      pret[i] = 0.0;
    }

    UNPROTECT(1);
    return(ret);
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    if (px[i] < 0.0 | px[i] >= 1.0) {
      pret[i] = 0.0;
    } else if (pb[0] == 1.0) {
      pret[i] = (pg[0] - 1) / R_pow(1 + (pg[0] - 1) * px[i], 2);
    } else if (pb[0] * pg[0] == 1.0) {
      pret[i] = -log(pb[0]) * R_pow(pb[0], px[i]);
    } else {
      double b1x = R_pow(pb[0], 1 - px[i]);
      pret[i] = (pb[0] - 1) * (pg[0] - 1) * b1x * log(pb[0]) /
        R_pow((pg[0] - 1) * b1x + 1 - pb[0] * pg[0], 2);
    }
  }

  if (pl == 1) {
    for (R_xlen_t i = 0; i < n; ++i) {
      pret[i] = log(pret[i]);
    }
  }

    UNPROTECT(1);
    return(ret);
}

static const R_CallMethodDef CallEntries[] = {
  {"dmb_c",    (DL_FUNC) &dmb_c,  4},
  {NULL,       NULL,              0}
};

void R_init_MBBEFDLite(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
