// Copyright Avraham Adler (c) 2026
// SPDX-License-Identifier: MPL-2.0+

#include <Rmath.h>
#include "MBBEFDLite.h"

// Developed with help of Claude. Apparently similar to Cephes and Boost.
SEXP dilog_c(SEXP x_) {
  const R_xlen_t n = xlength(x_);
  double *px = REAL(x_);
  const double PISQ_6 = M_PI * M_PI / 6.0;
  double ans;

  SEXP ret = PROTECT(allocVector(REALSXP, n));
  double *pret = REAL(ret);

  for (R_xlen_t i = 0; i < n; ++i) {
    if (fabs(px[i]) < DBL_EPSILON) {
      ans = 0.0;
    } else if (fabs(px[i] - 1.0) < DBL_EPSILON) {
      ans = PISQ_6;
    } else if (fabs(px[i] + 1.0) < DBL_EPSILON) {
      ans = -PISQ_6 * 0.5;
    } else {
      double y;
      double r;
      double s;
      // Move x into "proper" region for series
      if (px[i] > 1.0) {
        // Region 1: Re(Li₂(x)) = π²/3 − ½(ln x)² − Li₂(1/x)
        y = 1.0 / px[i];
        r = 2 * PISQ_6 - 0.5 * log(px[i]) * log(px[i]);
        s = -1.0;
        if (y > 0.5) {
          // Use Euler reflection in y: Li₂(y) = π²/6 − ln(y)ln(1−y) - Li₂(1−y)
          // Substituting: ans = r - [π²/6  ln(y)ln(1−y) - Li₂(1−y)]
          r = r - PISQ_6 + log(y) * log1p(-y);
          y = 0.5 - y + 0.5;   // = (x-1)/x, tiny for x near 1
          s = 1.0;
        }
      } else if (px[i] > 0.5) {
        // Region 2: Li₂(x) = π²/6 − ln(x)ln(1−x) - Li₂(1−x)
        y = 0.5 - px[i] + 0.5;
        r = PISQ_6 - log(px[i]) * log1p(-px[i]);
        s = -1.0;
      } else if (px[i] >= -1.0) {
        if (px[i] < 0.0) {
          // Region 3a: x ∈ [-1.0, 0]
          // Use Li₂(x) = -½ln²(1-x) - Li₂(x/(x-1)), maps to y ∈ [0, 0.5]
          y = px[i] / (px[i] - 1.0);
          r = -0.5 * log1p(-px[i]) * log1p(-px[i]);
          s = -1.0;
        } else {
          // Region 3b: x ∈ [0, 0.5] Use Li₂(x) = Sum(x^k / k²)
          y = px[i];
          r = 0.0;
          s = 1.0;
        }
      } else {
        // Region 4: Li₂(x) = −π²/6 − ½ln²(−x) − Li₂(1/x)
        // y = 1/x is always ∈ (−1, 0), so always apply:
        // Li₂(y) = −½ln²(1−y) − Li₂(y/(y−1))
        // Combined: r absorbs both log terms, z = 1/(1−x) ∈ (0, 0.5)
        double ly = log1p(-1.0 / px[i]);
        r = -PISQ_6 - 0.5 * log(-px[i]) * log(-px[i]);
        r += 0.5 * ly * ly;
        y = 1.0 / (1.0 - px[i]);
        s = 1.0;
      }
      double term = y;
      ans = y;
      for (int k = 2; k < 1000; ++k) {
        double tk = (double)(k - 1) / (double)k;
        term *= y * tk * tk;
        ans += term;
        if (fabs(term) < DBL_EPSILON * fabs(ans)) break;
      }
      ans = s * ans + r;
    }
    pret[i] = ans;
  }
  UNPROTECT(1);
  return ret;
}
