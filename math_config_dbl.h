/* Configuration for double precision math routines.
   Copyright (C) 2018-2019 Free Software Foundation, Inc.
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

#ifndef _MATH_CONFIG_H
#define _MATH_CONFIG_H
#define math_opt_barrier(x)                                     \
  ({ __typeof (x) __x = (x); __asm ("" : "+m" (__x)); __x; })
#define math_force_eval(x)                                              \
  ({ __typeof (x) __x = (x); __asm __volatile__ ("" : : "m" (__x)); })
#if FLT_EVAL_METHOD == 0
# define math_narrow_eval(x) (x)
#else
# if FLT_EVAL_METHOD == 1
#  define excess_precision(type) __builtin_types_compatible_p (type, float)
# else
#  define excess_precision(type) (__builtin_types_compatible_p (type, float) \
                                  || __builtin_types_compatible_p (type, \
                                                                   double))
# endif
# define math_narrow_eval(x)                                    \
  ({                                                            \
    __typeof (x) math_narrow_eval_tmp = (x);                    \
    if (excess_precision (__typeof (math_narrow_eval_tmp)))     \
      __asm__ ("" : "+m" (math_narrow_eval_tmp));               \
    math_narrow_eval_tmp;                                       \
   })
#endif


#define attribute_hidden
#include <math.h>
#include <stdint.h>
#include <altivec.h>

#ifndef WANT_ROUNDING
/* Correct special case results in non-nearest rounding modes.  */
# define WANT_ROUNDING 1
#endif
#ifndef WANT_ERRNO
/* Set errno according to ISO C with (math_errhandling & MATH_ERRNO) != 0.  */
# define WANT_ERRNO 1
#endif
#ifndef WANT_ERRNO_UFLOW
/* Set errno to ERANGE if result underflows to 0 (in all rounding modes).  */
# define WANT_ERRNO_UFLOW (WANT_ROUNDING && WANT_ERRNO)
#endif

#ifndef TOINT_INTRINSICS
/* When set, the roundtoint and converttoint functions are provided with
   the semantics documented below.  */
# define TOINT_INTRINSICS 0
#endif

#if TOINT_INTRINSICS
/* Round x to nearest int in all rounding modes, ties have to be rounded
   consistently with converttoint so the results match.  If the result
   would be outside of [-2^31, 2^31-1] then the semantics is unspecified.  */
static inline vector double
roundtoint (vector double x);

/* Convert x to nearest int in all rounding modes, ties have to be rounded
   consistently with roundtoint.  If the result is not representible in an
   int32_t then the semantics is unspecified.  */
static inline vector double
converttoint (double x);
#endif

static inline uint64_t
asuint64 (double f)
{
  union
  {
    double f;
    uint64_t i;
  } u = {f};
  return u.i;
}

static inline double
asdouble (uint64_t i)
{
  union
  {
    uint64_t i;
    double f;
  } u = {i};
  return u.f;
}

static inline int
issignaling_inline (double x)
{
  uint64_t ix = asuint64 (x);
  return 2 * (ix ^ 0x0008000000000000) > 2 * 0x7ff8000000000000ULL;
}

#define NOINLINE __attribute__ ((noinline))

/* NOINLINE reduces code size.  */
NOINLINE static double
xflow (uint32_t sign, double y)
{
  y = math_opt_barrier (sign ? -y : y) * y;
  return y;
}
/*
double
__math_uflow (uint32_t sign)
{
  return xflow (sign, 0x1p-767);
}

double
__math_oflow (uint32_t sign)
{
  return xflow (sign, 0x1p769);
}

*/
/* Error handling tail calls for special cases, with a sign argument.
   The sign of the return value is set if the argument is non-zero.  */

/* The result overflows.  */
//attribute_hidden double __math_oflow (uint32_t);
/* The result underflows to 0 in nearest rounding mode.  */
//attribute_hidden double __math_uflow (uint32_t);
/* The result underflows to 0 in some directed rounding mode only.  */
attribute_hidden double __math_may_uflow (uint32_t);
/* Division by zero.  */
attribute_hidden double __math_divzero (uint32_t);

/* Error handling using input checking.  */

/* Invalid input unless it is a quiet NaN.  */
attribute_hidden double __math_invalid (double);

/* Error handling using output checking, only for errno setting.  */


#define EXP_TABLE_BITS 7
#define EXP_POLY_ORDER 5
#define EXP2_POLY_ORDER 5
extern const struct exp_data
{
  vector double invln2N;
  vector double shift;
  vector double negln2hiN;
  vector double negln2loN;
  vector double poly[4]; /* Last four coefficients.  */
  vector double exp2_shift;
  vector double exp2_poly[EXP2_POLY_ORDER];
  uint64_t tab[2*(1 << EXP_TABLE_BITS)];
} __exp_data attribute_hidden;

#define LOG_TABLE_BITS 7
#define LOG_POLY_ORDER 6
#define LOG_POLY1_ORDER 12
extern const struct log_data
{
  vector double ln2hi;
  vector double ln2lo;
  vector double poly[LOG_POLY_ORDER - 1]; /* First coefficient is 1.  */
  vector double poly1[LOG_POLY1_ORDER - 1];
  /* See e_log_data.c for details.  */
  struct {double invc, logc;} tab[1 << LOG_TABLE_BITS];
  struct {double chi, clo;} tab2[1 << LOG_TABLE_BITS];
} __log_data attribute_hidden;

#define LOG2_TABLE_BITS 6
#define LOG2_POLY_ORDER 7
#define LOG2_POLY1_ORDER 11
extern const struct log2_data
{
  double invln2hi;
  double invln2lo;
  double poly[LOG2_POLY_ORDER - 1];
  double poly1[LOG2_POLY1_ORDER - 1];
  /* See e_log2_data.c for details.  */
  struct {double invc, logc;} tab[1 << LOG2_TABLE_BITS];
#ifndef __FP_FAST_FMA
  struct {double chi, clo;} tab2[1 << LOG2_TABLE_BITS];
#endif
} __log2_data attribute_hidden;

#define POW_LOG_TABLE_BITS 7
#define POW_LOG_POLY_ORDER 8
extern const struct pow_log_data
{
  double ln2hi;
  double ln2lo;
  double poly[POW_LOG_POLY_ORDER - 1]; /* First coefficient is 1.  */
  /* Note: the pad field is unused, but allows slightly faster indexing.  */
  /* See e_pow_log_data.c for details.  */
  struct {double invc, pad, logc, logctail;} tab[1 << POW_LOG_TABLE_BITS];
} __pow_log_data attribute_hidden;

#endif
