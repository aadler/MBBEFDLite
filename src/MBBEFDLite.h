// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#ifndef MBBEFDLite_H
#define MBBEFDLite_H

#include <R.h>
#include <Rinternals.h>

extern SEXP dmb_c(SEXP x, SEXP g, SEXP b, SEXP lg);
extern SEXP pmb_c(SEXP q, SEXP g, SEXP b, SEXP lt, SEXP lg);
extern SEXP qmb_c(SEXP q, SEXP g, SEXP b, SEXP lt, SEXP lg);
extern SEXP rmb_c(SEXP n_, SEXP g, SEXP b);
extern SEXP ecmb_c(SEXP x, SEXP g, SEXP b, SEXP lower_tail);

#endif
