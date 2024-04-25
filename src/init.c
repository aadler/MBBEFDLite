// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "dpqrmb.h"

static const R_CallMethodDef CallEntries[] = {
  {"dmb_c",    (DL_FUNC) &dmb_c,  4},
  {"pmb_c",    (DL_FUNC) &pmb_c,  5},
  {"qmb_c",    (DL_FUNC) &qmb_c,  5},
  {NULL,       NULL,              0}
};

void R_init_MBBEFDLite(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
