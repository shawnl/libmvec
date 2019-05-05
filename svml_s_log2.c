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
// Ported from musl, which is licensed under the MIT license:
// https://git.musl-libc.org/cgit/musl/tree/COPYRIGHT
//
// https://git.musl-libc.org/cgit/musl/tree/src/math/log.c

#define LOG_TABLE_BITS 7
#define LOG_POLY_ORDER 6
#define LOG_POLY1_ORDER 12
const struct log_data {
	double ln2hi;
	double ln2lo;
	double poly[LOG_POLY_ORDER - 1]; /* First coefficient is 1.  */
	double poly1[LOG_POLY1_ORDER - 1];
	struct {
		double invc, logc;
	} tab[1 << LOG_TABLE_BITS];
} __log_data;

#define N (1 << LOG_TABLE_BITS)

const struct log_data __log_data = {
.ln2hi = 0x1.62e42fefa3800p-1,
.ln2lo = 0x1.ef35793c76730p-45,
.poly1 = {
// relative error: 0x1.c04d76cp-63
// in -0x1p-4 0x1.09p-4 (|log(1+x)| > 0x1p-4 outside the interval)
-0x1p-1,
0x1.5555555555577p-2,
-0x1.ffffffffffdcbp-3,
0x1.999999995dd0cp-3,
-0x1.55555556745a7p-3,
0x1.24924a344de3p-3,
-0x1.fffffa4423d65p-4,
0x1.c7184282ad6cap-4,
-0x1.999eb43b068ffp-4,
0x1.78182f7afd085p-4,
-0x1.5521375d145cdp-4,
},
.poly = {
// relative error: 0x1.926199e8p-56
// abs error: 0x1.882ff33p-65
// in -0x1.fp-9 0x1.fp-9
-0x1.0000000000001p-1,
0x1.555555551305bp-2,
-0x1.fffffffeb459p-3,
0x1.999b324f10111p-3,
-0x1.55575e506c89fp-3,
},
/* Algorithm:

	x = 2^k z
	log(x) = k ln2 + log(c) + log(z/c)
	log(z/c) = poly(z/c - 1)

where z is in [1.6p-1; 1.6p0] which is split into N subintervals and z falls
into the ith one, then table entries are computed as

	tab[i].invc = 1/c
	tab[i].logc = (double)log(c)
	tab2[i].chi = (double)c
	tab2[i].clo = (double)(c - (double)c)

where c is near the center of the subinterval and is chosen by trying +-2^29
floating point invc candidates around 1/center and selecting one for which

	1) the rounding error in 0x1.8p9 + logc is 0,
	2) the rounding error in z - chi - clo is < 0x1p-66 and
	3) the rounding error in (double)log(c) is minimized (< 0x1p-66).

Note: 1) ensures that k*ln2hi + logc can be computed without rounding error,
2) ensures that z/c - 1 can be computed as (z - chi - clo)*invc with close to
a single rounding error when there is no fast fma for z*invc - 1, 3) ensures
that logc + poly(z/c - 1) has small error, however near x == 1 when
|log(x)| < 0x1p-4, this is not enough so that is special cased.  */
.tab = {
{0x1.734f0c3e0de9fp+0, -0x1.7cc7f79e69000p-2},
{0x1.713786a2ce91fp+0, -0x1.76feec20d0000p-2},
{0x1.6f26008fab5a0p+0, -0x1.713e31351e000p-2},
{0x1.6d1a61f138c7dp+0, -0x1.6b85b38287800p-2},
{0x1.6b1490bc5b4d1p+0, -0x1.65d5590807800p-2},
{0x1.69147332f0cbap+0, -0x1.602d076180000p-2},
{0x1.6719f18224223p+0, -0x1.5a8ca86909000p-2},
{0x1.6524f99a51ed9p+0, -0x1.54f4356035000p-2},
{0x1.63356aa8f24c4p+0, -0x1.4f637c36b4000p-2},
{0x1.614b36b9ddc14p+0, -0x1.49da7fda85000p-2},
{0x1.5f66452c65c4cp+0, -0x1.445923989a800p-2},
{0x1.5d867b5912c4fp+0, -0x1.3edf439b0b800p-2},
{0x1.5babccb5b90dep+0, -0x1.396ce448f7000p-2},
{0x1.59d61f2d91a78p+0, -0x1.3401e17bda000p-2},
{0x1.5805612465687p+0, -0x1.2e9e2ef468000p-2},
{0x1.56397cee76bd3p+0, -0x1.2941b3830e000p-2},
{0x1.54725e2a77f93p+0, -0x1.23ec58cda8800p-2},
{0x1.52aff42064583p+0, -0x1.1e9e129279000p-2},
{0x1.50f22dbb2bddfp+0, -0x1.1956d2b48f800p-2},
{0x1.4f38f4734ded7p+0, -0x1.141679ab9f800p-2},
{0x1.4d843cfde2840p+0, -0x1.0edd094ef9800p-2},
{0x1.4bd3ec078a3c8p+0, -0x1.09aa518db1000p-2},
{0x1.4a27fc3e0258ap+0, -0x1.047e65263b800p-2},
{0x1.4880524d48434p+0, -0x1.feb224586f000p-3},
{0x1.46dce1b192d0bp+0, -0x1.f474a7517b000p-3},
{0x1.453d9d3391854p+0, -0x1.ea4443d103000p-3},
{0x1.43a2744b4845ap+0, -0x1.e020d44e9b000p-3},
{0x1.420b54115f8fbp+0, -0x1.d60a22977f000p-3},
{0x1.40782da3ef4b1p+0, -0x1.cc00104959000p-3},
{0x1.3ee8f5d57fe8fp+0, -0x1.c202956891000p-3},
{0x1.3d5d9a00b4ce9p+0, -0x1.b81178d811000p-3},
{0x1.3bd60c010c12bp+0, -0x1.ae2c9ccd3d000p-3},
{0x1.3a5242b75dab8p+0, -0x1.a45402e129000p-3},
{0x1.38d22cd9fd002p+0, -0x1.9a877681df000p-3},
{0x1.3755bc5847a1cp+0, -0x1.90c6d69483000p-3},
{0x1.35dce49ad36e2p+0, -0x1.87120a645c000p-3},
{0x1.34679984dd440p+0, -0x1.7d68fb4143000p-3},
{0x1.32f5cceffcb24p+0, -0x1.73cb83c627000p-3},
{0x1.3187775a10d49p+0, -0x1.6a39a9b376000p-3},
{0x1.301c8373e3990p+0, -0x1.60b3154b7a000p-3},
{0x1.2eb4ebb95f841p+0, -0x1.5737d76243000p-3},
{0x1.2d50a0219a9d1p+0, -0x1.4dc7b8fc23000p-3},
{0x1.2bef9a8b7fd2ap+0, -0x1.4462c51d20000p-3},
{0x1.2a91c7a0c1babp+0, -0x1.3b08abc830000p-3},
{0x1.293726014b530p+0, -0x1.31b996b490000p-3},
{0x1.27dfa5757a1f5p+0, -0x1.2875490a44000p-3},
{0x1.268b39b1d3bbfp+0, -0x1.1f3b9f879a000p-3},
{0x1.2539d838ff5bdp+0, -0x1.160c8252ca000p-3},
{0x1.23eb7aac9083bp+0, -0x1.0ce7f57f72000p-3},
{0x1.22a012ba940b6p+0, -0x1.03cdc49fea000p-3},
{0x1.2157996cc4132p+0, -0x1.f57bdbc4b8000p-4},
{0x1.201201dd2fc9bp+0, -0x1.e370896404000p-4},
{0x1.1ecf4494d480bp+0, -0x1.d17983ef94000p-4},
{0x1.1d8f5528f6569p+0, -0x1.bf9674ed8a000p-4},
{0x1.1c52311577e7cp+0, -0x1.adc79202f6000p-4},
{0x1.1b17c74cb26e9p+0, -0x1.9c0c3e7288000p-4},
{0x1.19e010c2c1ab6p+0, -0x1.8a646b372c000p-4},
{0x1.18ab07bb670bdp+0, -0x1.78d01b3ac0000p-4},
{0x1.1778a25efbcb6p+0, -0x1.674f145380000p-4},
{0x1.1648d354c31dap+0, -0x1.55e0e6d878000p-4},
{0x1.151b990275fddp+0, -0x1.4485cdea1e000p-4},
{0x1.13f0ea432d24cp+0, -0x1.333d94d6aa000p-4},
{0x1.12c8b7210f9dap+0, -0x1.22079f8c56000p-4},
{0x1.11a3028ecb531p+0, -0x1.10e4698622000p-4},
{0x1.107fbda8434afp+0, -0x1.ffa6c6ad20000p-5},
{0x1.0f5ee0f4e6bb3p+0, -0x1.dda8d4a774000p-5},
{0x1.0e4065d2a9fcep+0, -0x1.bbcece4850000p-5},
{0x1.0d244632ca521p+0, -0x1.9a1894012c000p-5},
{0x1.0c0a77ce2981ap+0, -0x1.788583302c000p-5},
{0x1.0af2f83c636d1p+0, -0x1.5715e67d68000p-5},
{0x1.09ddb98a01339p+0, -0x1.35c8a49658000p-5},
{0x1.08cabaf52e7dfp+0, -0x1.149e364154000p-5},
{0x1.07b9f2f4e28fbp+0, -0x1.e72c082eb8000p-6},
{0x1.06ab58c358f19p+0, -0x1.a55f152528000p-6},
{0x1.059eea5ecf92cp+0, -0x1.63d62cf818000p-6},
{0x1.04949cdd12c90p+0, -0x1.228fb8caa0000p-6},
{0x1.038c6c6f0ada9p+0, -0x1.c317b20f90000p-7},
{0x1.02865137932a9p+0, -0x1.419355daa0000p-7},
{0x1.0182427ea7348p+0, -0x1.81203c2ec0000p-8},
{0x1.008040614b195p+0, -0x1.0040979240000p-9},
{0x1.fe01ff726fa1ap-1, 0x1.feff384900000p-9},
{0x1.fa11cc261ea74p-1, 0x1.7dc41353d0000p-7},
{0x1.f6310b081992ep-1, 0x1.3cea3c4c28000p-6},
{0x1.f25f63ceeadcdp-1, 0x1.b9fc114890000p-6},
{0x1.ee9c8039113e7p-1, 0x1.1b0d8ce110000p-5},
{0x1.eae8078cbb1abp-1, 0x1.58a5bd001c000p-5},
{0x1.e741aa29d0c9bp-1, 0x1.95c8340d88000p-5},
{0x1.e3a91830a99b5p-1, 0x1.d276aef578000p-5},
{0x1.e01e009609a56p-1, 0x1.07598e598c000p-4},
{0x1.dca01e577bb98p-1, 0x1.253f5e30d2000p-4},
{0x1.d92f20b7c9103p-1, 0x1.42edd8b380000p-4},
{0x1.d5cac66fb5ccep-1, 0x1.606598757c000p-4},
{0x1.d272caa5ede9dp-1, 0x1.7da76356a0000p-4},
{0x1.cf26e3e6b2ccdp-1, 0x1.9ab434e1c6000p-4},
{0x1.cbe6da2a77902p-1, 0x1.b78c7bb0d6000p-4},
{0x1.c8b266d37086dp-1, 0x1.d431332e72000p-4},
{0x1.c5894bd5d5804p-1, 0x1.f0a3171de6000p-4},
{0x1.c26b533bb9f8cp-1, 0x1.067152b914000p-3},
{0x1.bf583eeece73fp-1, 0x1.147858292b000p-3},
{0x1.bc4fd75db96c1p-1, 0x1.2266ecdca3000p-3},
{0x1.b951e0c864a28p-1, 0x1.303d7a6c55000p-3},
{0x1.b65e2c5ef3e2cp-1, 0x1.3dfc33c331000p-3},
{0x1.b374867c9888bp-1, 0x1.4ba366b7a8000p-3},
{0x1.b094b211d304ap-1, 0x1.5933928d1f000p-3},
{0x1.adbe885f2ef7ep-1, 0x1.66acd2418f000p-3},
{0x1.aaf1d31603da2p-1, 0x1.740f8ec669000p-3},
{0x1.a82e63fd358a7p-1, 0x1.815c0f51af000p-3},
{0x1.a5740ef09738bp-1, 0x1.8e92954f68000p-3},
{0x1.a2c2a90ab4b27p-1, 0x1.9bb3602f84000p-3},
{0x1.a01a01393f2d1p-1, 0x1.a8bed1c2c0000p-3},
{0x1.9d79f24db3c1bp-1, 0x1.b5b515c01d000p-3},
{0x1.9ae2505c7b190p-1, 0x1.c2967ccbcc000p-3},
{0x1.9852ef297ce2fp-1, 0x1.cf635d5486000p-3},
{0x1.95cbaeea44b75p-1, 0x1.dc1bd3446c000p-3},
{0x1.934c69de74838p-1, 0x1.e8c01b8cfe000p-3},
{0x1.90d4f2f6752e6p-1, 0x1.f5509c0179000p-3},
{0x1.8e6528effd79dp-1, 0x1.00e6c121fb800p-2},
{0x1.8bfce9fcc007cp-1, 0x1.071b80e93d000p-2},
{0x1.899c0dabec30ep-1, 0x1.0d46b9e867000p-2},
{0x1.87427aa2317fbp-1, 0x1.13687334bd000p-2},
{0x1.84f00acb39a08p-1, 0x1.1980d67234800p-2},
{0x1.82a49e8653e55p-1, 0x1.1f8ffe0cc8000p-2},
{0x1.8060195f40260p-1, 0x1.2595fd7636800p-2},
{0x1.7e22563e0a329p-1, 0x1.2b9300914a800p-2},
{0x1.7beb377dcb5adp-1, 0x1.3187210436000p-2},
{0x1.79baa679725c2p-1, 0x1.377266dec1800p-2},
{0x1.77907f2170657p-1, 0x1.3d54ffbaf3000p-2},
{0x1.756cadbd6130cp-1, 0x1.432eee32fe000p-2},
},
};

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


vector double _ZGV9N2v_log(vector double x) {
  v64u special_cases;
  v64u inf = {INF, INF};
  v64u ninf = inf | (1ULL << 63);
  v64u zero = {0, 0};
  u zerou;
  zerou.l = zero;

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
  v64u is_close_to_one = vec_cmplt(xi - lov, hiv - lov);
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

    /* 1.0 -> 0 */
    u oned;
    vector double one = {1.0, 1.0};
    oned.d = one;
    v64u is_one = vec_cmpeq(xi, oned.l);
    res.l = vec_sel(res.l, zero, is_one);
  } else
    res.l = zero;

  v64u infexp = {0x7ff0000000000000, 0x7ff0000000000000};
  v64u is_special_cases = vec_cmpge(xi - 0x0010000000000000, infexp - 0x0010000000000000);
  if (!vec_all_eq(is_special_cases, zero)) {
    v64u is_zero = vec_cmpeq(xi << 1, zero);
    res.l = vec_sel(res.l, ninf, is_zero);

    v64u is_inf = vec_cmpeq(xi, inf);
    res.l = vec_sel(res.l, inf, is_inf);

    v64u is_neg = vec_cmpne(xi >> 63, zero);
    v64u is_nan = vec_cmpeq(xi & infexp, infexp) & ~is_inf;
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
  v64i k = (v64u)((v64i)tmp >> 52);
  v64u iz = xi - (tmp & 0xfffULL << 52);

  vector double invc = {T[i[0]].invc, T[i[1]].invc};
  vector double logc = {T[i[0]].logc, T[i[1]].logc};

  u z;
  z.l = iz;
  vector double neg1 = {-1.0, -1.0};
  vector double r = vec_madd(z.d, invc, neg1);
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
