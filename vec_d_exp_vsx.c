/* Single-precision vector expf(x) function.
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
#include <math.h>

#define TOINT_INTRINSICS 0
#define WANT_ROUNDING 1

#include "exp_data.h"
#include "math_config_flt.h"

//#include <sysdeps/ieee754/dbl-64/math_config.h>

typedef vector long long unsigned v64u;
typedef union {
  vector unsigned u;
  vector float f;
  v64u l;
  vector double d;
} u;
typedef union {
  double d;
  uint64_t l;
  unsigned u;
  float f;
} us;

#define N (1 << EXP_TABLE_BITS)
#define InvLn2N __exp_data.invln2N
#define NegLn2hiN __exp_data.negln2hiN
#define NegLn2loN __exp_data.negln2loN
#define Shift __exp_data.shift
#define T __exp_data.tab
#define C2 __exp_data.poly[5 - EXP_POLY_ORDER]
#define C3 __exp_data.poly[6 - EXP_POLY_ORDER]
#define C4 __exp_data.poly[7 - EXP_POLY_ORDER]
#define C5 __exp_data.poly[8 - EXP_POLY_ORDER]

/* Handle cases that may overflow or underflow when computing the result that
   is scale*(1+TMP) without intermediate rounding.  The bit representation of
   scale is in SBITS, however it has a computed exponent that may have
   overflown into the sign bit so that needs to be adjusted before using it as
   a double.  (int32_t)KI is the k used in the argument reduction and exponent
   adjustment of scale, positive k here means the result may overflow and
   negative k means the result may underflow.  */
static inline double specialcase(double_t tmp, uint64_t sbits, uint64_t ki)
{
	double_t scale, y;

	if ((ki & 0x80000000) == 0) {
		/* k > 0, the exponent of scale might have overflowed by <= 460.  */
		sbits -= 1009ull << 52;
		scale = asdouble(sbits);
		y = 0x1p1009 * (scale + scale * tmp);
		return eval_as_double(y);
	}
	/* k < 0, need special care in the subnormal range.  */
	sbits += 1022ull << 52;
	scale = asdouble(sbits);
	y = scale + scale * tmp;
	if (y < 1.0) {
		/* Round y to the right precision before scaling it into the subnormal
		 range to avoid double rounding that can cause 0.5+E/2 ulp error where
		 E is the worst-case ulp error outside the subnormal range.  So this
		 is only useful if the goal is better than 1 ulp worst-case error.  */
		double_t hi, lo;
		lo = scale - y + scale * tmp;
		hi = 1.0 + y;
		lo = 1.0 - hi + y + lo;
		y = eval_as_double(hi + lo) - 1.0;
		/* Avoid -0.0 with downward rounding.  */
		if (WANT_ROUNDING && y == 0.0)
			y = 0.0;
		/* The underflow exception needs to be signaled explicitly.  */
		fp_force_eval(fp_barrier(0x1p-1022) * 0x1p-1022);
	}
	y = 0x1p-1022 * y;
	return eval_as_double(y);
}

/* Top 12 bits of a double (sign and exponent bits).  */
static inline uint32_t top12(double x)
{
	return asuint64(x) >> 52;
}

vector double
_ZGVbN2v_exp (vector double x)
{
	v64u abstop, is_special_case, is_special_case2, ki, idx, top, sbits;
	vector double kd, z, r, r2, scale, tail, tmp;
	vector double InvLn2Nv = {InvLn2N, InvLn2N};
	vector double NegLn2hiNv = {NegLn2hiN, NegLn2hiN};
	vector double NegLn2loNv = {NegLn2loN, NegLn2loN};
	vector double Shiftv = {Shift, Shift};
	vector double c2 = {C2, C2};
	vector double c3 = {C3, C3};
	vector double c4 = {C4, C4};
	vector double c5 = {C5, C5};
	u t, t2;
	u res;
	u xu;
	xu.d = x;
	v64u zero = {0, 0};
	u c512;
	c512.d = {512.0, 512.0};
	u normalize;
	normalize.d = {0x1p-54, 0x1p-54};
	c512.l &= 0x7ff0000000000000;
	normalize.l &= 0x7ff0000000000000;
	abstop = xu.l & 0x7ff0000000000000;
	is_special_case = vec_cmpge (xu.l - normalize.l, c512 - normalize.l);
	if (__glibc_unlikely (!vec_all_eq (is_special_case, zero)))
	{
		for (int m=0;m<2;++m)
		{
			double v;
			uint32_t abstops = abstop[m] >> 52;
			if (!is_special_case[m])
				continue;
			v = x[m];
			if (abstops - top12(0x1p-54) >= 0x80000000)
				/* Avoid spurious underflow for tiny x.  */
				/* Note: 0 is common input.  */
				res.d[m] = WANT_ROUNDING ? 1.0 + x : 1.0;
			if (abstops >= top12(1024.0)) {
				if (asuint64(v) == asuint64(-INFINITY))
					res.d[m] = 0.0;
				if (abstops >= top12(INFINITY))
					res.d[m] = 1.0 + v;
				if (asuint64(v) >> 63)
					res.d[m] = __math_uflow(0);
				else
					res.d[m] = __math_oflow(0);
			} else {
	 			/* Large x is special cased below.  */
 				abstop[m] = 0;
				is_special_case[m] = 0;
			}
		}
	}
	/* exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].  */
	/* x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].  */
	z = InvLn2Nv * x;
#if TOINT_INTRINSICS
	kd = roundtoint(z);
	ki = converttoint(z);
#elif EXP_USE_TOINT_NARROW
	/* z - kd is in [-0.5-2^-16, 0.5] in all rounding modes.  */
	kd = z + Shift;
	u kdu;
	kdu.d = kd;
	kdu.l >>= 16;
	kd = (vector double)(vector int)kdu.l;
#else
	/* z - kd is in [-1, 1] in non-nearest rounding modes.  */
	kd = z + Shift;
	u kdu;
	kdu.d = kd;
	kd -= Shift;
#endif
	r = x + kd * NegLn2hiNv + kd * NegLn2loNv;
	/* 2^(k/N) ~= scale * (1 + tail).  */
	idx = 2 * (kdu.l % N);
	top = ki << (52 - EXP_TABLE_BITS);
	t.l = {T[idx, T[idx]};
	t2.l = {T[idx + 1], T[idx + 1]};
	tail = t.d;
	/* This is only a valid scale when -1023*N < k < 1024*N.  */
	sbits = t2.l + top;
	/* exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).  */
	/* Evaluation is optimized assuming superscalar pipelined execution.  */
	r2 = r * r;
	/* Without fma the worst case error is 0.25/N ulp larger.  */
	/* Worst case error is less than 0.5+1.11/N+(abs poly error * 2^53) ulp.  */
	tmp = tail + r + r2 * (c2 + r * c3) + r2 * r2 * (c4 + r * c5);
	is_special_case2 = vec_cmpeq (abstop, zero);
	if (__glibc_unlikely (!vec_all_eq (is_special_case2, zero)))
	{
		u sc2;
		for (int m=0;m<2;++m)
		{
			sc2.d[m] = specialcase(tmp[m], sbits[m], ki[m]);
		}
		res.l = vec_sel (res.l, sc2.l, is_special_case2);
	}
	u sbitsu;
	sbitsu.l = sbits;
	scale = sbitsu.d;
	u res2;
	res2.d = scale + scale * tmp;
	/* Note: tmp == 0 or |tmp| > 2^-200 and scale > 2^-739, so there
	   is no spurious underflow here even without fma.  */
	return vec_sel (res.l, res2.l, ~is_special_case & ~is_special_case2);
}
