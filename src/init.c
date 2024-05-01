// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "MBBEFDLite.h"

static const R_CallMethodDef CallEntries[] = {
  {"dmb_c",     (DL_FUNC) &dmb_c,   4},
  {"pmb_c",     (DL_FUNC) &pmb_c,   5},
  {"qmb_c",     (DL_FUNC) &qmb_c,   5},
  {"rmb_c",     (DL_FUNC) &rmb_c,   3},
  {"ecmb_c",    (DL_FUNC) &ecmb_c,  4},
  {NULL,        NULL,               0}
};

void R_init_MBBEFDLite(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
  R_RegisterCCallable("MBBEFDLite", "dmb_c",   (DL_FUNC) &dmb_c);
  R_RegisterCCallable("MBBEFDLite", "pmb_c",   (DL_FUNC) &pmb_c);
  R_RegisterCCallable("MBBEFDLite", "qmb_c",   (DL_FUNC) &qmb_c);
  R_RegisterCCallable("MBBEFDLite", "rmb_c",   (DL_FUNC) &rmb_c);
  R_RegisterCCallable("MBBEFDLite", "ecmb_c",  (DL_FUNC) &ecmb_c);
}
