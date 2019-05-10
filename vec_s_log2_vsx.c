/* Double-precision vector log(x) function.
   Copyright (C) 2019 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
#include <altivec.h>
#include <stdint.h>

//#include <sysdeps/ieee754/dbl-64/math_config.h>
#define __FP_FAST_FMA 1
#include "data.h"

typedef vector long long unsigned v64u;
typedef vector long long v64i;

typedef union {
  vector double d;
  v64u l;
} u;

typedef union {
  double d;
  uint64_t l;
  unsigned u;
  float f;
} us;

#define T __log_data.tab
#define T2 __log_data.tab2
#define B __log_data.poly1
#define A __log_data.poly
#define Ln2hi __log_data.ln2hi
#define Ln2lo __log_data.ln2lo
#define N (1 << LOG_TABLE_BITS)
#define OFF 0x3fe6000000000000

#define INF 0x7ff0000000000000

vector double _ZGVbN2v_log(vector double x) {
  v64u inf = {INF, INF};
  v64u ninf = inf | (1ULL << 63);
  v64u zero = {0, 0};

  u un;
  un.d = x;
  v64u xi = un.l;

  us lo;
  us hi;
  lo.d = 1.0 - 0x1p-4;
  hi.d = 1.0 + 0x1.09p-4;
  v64u lov = {lo.l, lo.l};
  v64u hiv = {hi.l, hi.l};
  u res;
  v64u is_close_to_one = (v64u)vec_cmplt(xi - lov, hiv - lov);
  if (!vec_all_eq(is_close_to_one, zero)) {
    vector double r = x - 1.0;
    vector double r2 = r * r;
    vector double r3 = r * r2;
    vector double b0 = {B[0], B[0]};
    vector double b12 = {B[1], B[2]};
    vector double b1 = vec_splat(b12, 0);
    vector double b2 = vec_splat(b12, 1);
    vector double b34 = {B[3], B[4]};
    vector double b3 = vec_splat(b34, 0);
    vector double b4 = vec_splat(b34, 1);
    vector double b56 = {B[5], B[6]};
    vector double b5 = vec_splat(b56, 0);
    vector double b6 = vec_splat(b56, 1);
    vector double b78 = {B[7], B[8]};
    vector double b7 = vec_splat(b78, 0);
    vector double b8 = vec_splat(b78, 1);
    vector double b910 = {B[9], B[10]};
    vector double b9 = vec_splat(b910, 0);
    vector double b10 = vec_splat(b910, 1);

    res.d = r3 * (b1 + r * b2 + r2 * b3 + r3 *
      (b4 + r * b5 + r2 * b6 + r3 *
        (b7 + r * b8 + r2 * b9 + r3 * b10)));
    /* Worst-case error is around 0.507 ULP.  */
    vector double w = r * 0x1p27;
    vector double rhi = r + w - w;
    vector double rlo = r - rhi;
    w = rhi * rhi * b0; /* B[0] == -0.5 */
    vector double hi = r + w;
    vector double lo = r - hi + w;
    lo += b0 * rlo * (rhi + r);
    res.d += lo;
    res.d += hi;

#if WANT_ROUNDING
    /* 1.0 -> 0 */
    u oned;
    vector double one = {1.0, 1.0};
    oned.d = one;
    v64u is_one = (v64u)vec_cmpeq(xi, oned.l);
    res.l = vec_sel(res.l, zero, is_one);
#endif
  } else
    res.l = zero;

  v64u infexp = {0x7ff0000000000000, 0x7ff0000000000000};
  v64u is_special_cases = (v64u)vec_cmpge(xi - 0x0010000000000000, infexp - 0x0010000000000000);
  if (!vec_all_eq(is_special_cases, zero)) {
    v64u is_zero = (v64u)vec_cmpeq(xi << 1, zero);
    res.l = vec_sel(res.l, ninf, is_zero);

    v64u is_inf = (v64u)vec_cmpeq(xi, inf);
    res.l = vec_sel(res.l, inf, is_inf);

    v64u is_neg = (v64u)vec_cmpne(xi >> 63, zero);
    v64u is_nan = (v64u)vec_cmpeq(xi & infexp, infexp) & ~is_inf;
    res.l = vec_sel(res.l, infexp + 1/*NaN*/, is_nan | (is_neg & ~is_zero));

    /* subnormals: normalize, remove from is_special_cases */
    v64u is_not_subnormal = is_nan | is_neg | is_zero | is_inf;
    un.d = x * 0x1p52;
    xi = vec_sel(xi, un.l - (52ULL << 52), is_special_cases & ~is_not_subnormal);
    is_special_cases = is_not_subnormal;
  }
  /* x = 2^k z; where z is in range [OFF,2*OFF) and exact.
	   The range is split into N subintervals.
	   The ith subinterval contains z and c is near its center.  */
  v64u tmp = xi - OFF;
  v64u i = (tmp >> (52 - LOG_TABLE_BITS)) % N;
  v64i k = ((v64i)tmp >> 52);
  v64u iz = xi - (tmp & 0xfffULL << 52);

  vector double invc = {T[i[0]].invc, T[i[1]].invc};
  vector double logc = {T[i[0]].logc, T[i[1]].logc};

  u z;
  z.l = iz;
/* VSX has a fast fma, but we do this so that we get the same result as log()*/
#ifdef __FP_FAST_FMA
  vector double neg1 = {-1.0, -1.0};
  vector double r = vec_madd(z.d, invc, neg1);
#else
  vector double chi = {T2[i[0]].chi, T2[i[1]].chi};
  vector double clo = {T2[i[0]].clo, T2[i[1]].clo};
  vector double r = (z - chi - clo) * invc;
#endif
  vector double ln2hi = {Ln2hi, Ln2hi};
  vector double ln2lo = {Ln2lo, Ln2lo};
  vector double kd = {(double)k[0], (double)k[1]};
  vector double w2 = kd * ln2hi + logc;
  vector double hid = w2 + r;
  vector double lod = w2 - hid + r + kd * ln2lo;


  /* log(x) = lo + (log1p(r) - r) + hi.  */
  vector double r2 = r * r; /* rounding error: 0x1p-54/N^2.  */

  /* Worst case error if |y| > 0x1p-5:
  0.5 + 4.13/N + abs-poly-error*2^57 ULP
  Worst case error if |y| > 0x1p-4:
  0.5 + 2.06/N + abs-poly-error*2^56 ULP.  */
  vector double a0 = {A[0], A[0]};
  vector double a1 = {A[1], A[1]};
  vector double a2 = {A[2], A[2]};
  vector double a3 = {A[3], A[3]};
  vector double a4 = {A[4], A[4]};
  vector double y = lod + r2 * a0 + r * r2 *
    (a1 + r * a2 + r2 *
      (a3 + r * a4)) + hid;
  res.d = vec_sel(res.d, y, ~is_close_to_one & ~is_special_cases);
  return res.d;
}
