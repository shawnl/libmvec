// Ported from musl, which is licensed under the MIT license:
// https://git.musl-libc.org/cgit/musl/tree/COPYRIGHT
//
// https://git.musl-libc.org/cgit/musl/tree/src/math/logf.c
/*
LOGF_TABLE_BITS = 4
LOGF_POLY_ORDER = 4

ULP error: 0.818 (nearest rounding.)
Relative error: 1.957 * 2^-26 (before rounding.)
*/
#define INF 0x7F800000
#define NINF 0xFF800000
#define LOGF_TABLE_BITS 4
#define LOGF_POLY_ORDER 4
const struct logf_data {
        struct {
                double invc, logc;
        } tab[1 << LOGF_TABLE_BITS];
        double ln2;
        double poly[LOGF_POLY_ORDER - 1]; /* First order coefficient is 1.  */
} __logf_data;

const struct logf_data __logf_data = {
  .tab = {
  { 0x1.661ec79f8f3bep+0, -0x1.57bf7808caadep-2 },
  { 0x1.571ed4aaf883dp+0, -0x1.2bef0a7c06ddbp-2 },
  { 0x1.49539f0f010bp+0, -0x1.01eae7f513a67p-2 },
  { 0x1.3c995b0b80385p+0, -0x1.b31d8a68224e9p-3 },
  { 0x1.30d190c8864a5p+0, -0x1.6574f0ac07758p-3 },
  { 0x1.25e227b0b8eap+0, -0x1.1aa2bc79c81p-3 },
  { 0x1.1bb4a4a1a343fp+0, -0x1.a4e76ce8c0e5ep-4 },
  { 0x1.12358f08ae5bap+0, -0x1.1973c5a611cccp-4 },
  { 0x1.0953f419900a7p+0, -0x1.252f438e10c1ep-5 },
  { 0x1p+0, 0x0p+0 },
  { 0x1.e608cfd9a47acp-1, 0x1.aa5aa5df25984p-5 },
  { 0x1.ca4b31f026aap-1, 0x1.c5e53aa362eb4p-4 },
  { 0x1.b2036576afce6p-1, 0x1.526e57720db08p-3 },
  { 0x1.9c2d163a1aa2dp-1, 0x1.bc2860d22477p-3 },
  { 0x1.886e6037841edp-1, 0x1.1058bc8a07ee1p-2 },
  { 0x1.767dcf5534862p-1, 0x1.4043057b6ee09p-2 },
  },
  .ln2 = 0x1.62e42fefa39efp-1,
  .poly = {
  -0x1.00ea348b88334p-2, 0x1.5575b0be00b6ap-2, -0x1.ffffef20a4123p-2,
  }
};

#include <altivec.h>

typedef union {
  vector float f;
  vector unsigned i;
  vector int s;
  vector double d;
  vector long long unsigned l;
} u;

#define T __logf_data.tab
#define A __logf_data.poly
#define Ln2 __logf_data.ln2
#define N (1 << LOGF_TABLE_BITS)
#define OFF 0x3f330000

vector float _ZGV9N4v_logf(vector float x) {
  vector unsigned special_cases;
  vector unsigned inf = {INF, INF, INF, INF};
  vector unsigned ninf = inf | (1 << 31);
  vector unsigned zero = {0, 0, 0, 0};
  u zerou;
  zerou.i = zero;

  u un;
  un.f = x;
  vector unsigned xi = un.i;
  vector unsigned if_special_cases = vec_cmpge(xi - 0x00800000, inf - 0x00800000);
  if (!vec_all_eq(if_special_cases, zero)) {
    // 0 pos -> -inf, 0 neg -> NaN
    vector unsigned is_zero = vec_cmpeq(xi << 1, zero);
    special_cases = is_zero & ninf;

    // inf -> inf
    vector unsigned is_inf = vec_cmpeq(xi, inf);;
    special_cases = vec_sel(special_cases, inf, is_inf);

    // Invalid
    vector unsigned is_negative = vec_cmpne(xi & 0x80000000, zero);
    vector unsigned splat = {0xff000000, 0xff000000, 0xff000000, 0xff000000};
    vector unsigned is_outofrange = vec_cmpge(xi << 1, splat);
    vector unsigned is_invalid = (is_negative | is_outofrange) & ~is_inf & ~is_zero;
    vector unsigned nan = {0x7f800001, 0x7f800001, 0x7f800001, 0x7f800001};
    special_cases = vec_sel(special_cases, nan, is_invalid);

    // normalize subnormals
    vector unsigned is_subnormal = ((if_special_cases & ~is_zero) & ~is_inf) & ~is_invalid;
    vector unsigned subnormals = is_subnormal & xi;
    un.i = subnormals;
    u un2;
    un2.f = un.f * 0x1p23f;
    subnormals = un2.i - 23 << 23;

    // clear subnormals, and merge in normalized ones
    xi = vec_sel(xi, subnormals, is_subnormal);
  }

  // back to the main part
  vector unsigned tmp = xi - OFF;
  vector unsigned i = (tmp >> (23 - LOGF_TABLE_BITS)) % N;
  vector unsigned k = (vector unsigned)((vector int)tmp >> 23);
  vector unsigned iz = xi - (tmp & 0x1ff << 23);

  vector double abinv = {(double)T[i[0]].invc, (double)T[i[1]].invc};
  vector double ablogc = {(double)T[i[0]].logc, (double)T[i[1]].logc};
  vector double cdinv = {(double)T[i[2]].invc, (double)T[i[3]].invc};
  vector double cdlogc = {(double)T[i[2]].logc, (double)T[i[3]].logc};

  u izu;
  izu.i = iz;
  vector double izl = {(double)izu.f[0], (double)izu.f[1]};
  vector double izr = {(double)izu.f[2], (double)izu.f[3]};
  vector double rl = izl * abinv - 1;
  vector double rr = izr * cdinv - 1;
  vector double kl = {(double)k[0], (double)k[1]};
  vector double kr = {(double)k[2], (double)k[3]};
  vector double Ln2v = {Ln2, Ln2};
  vector double y0l = ablogc + kl * Ln2v;
  vector double y0r = cdlogc + kr * Ln2v;

  vector double r2l = rl * rl;
  vector double r2r = rr * rr;
  vector double yl = A[1] * rl + A[2];
  vector double yr = A[1] * rr + A[2];
  yl = yl * r2l + (y0l + rl);
  yr = yr * r2r + (y0r + rr);
  vector float y = {(float)yl[0], (float)yl[1], (float)yr[0], (float)yr[1]};
  un.f = y;
  un.i = vec_sel(un.i, special_cases, if_special_cases);
  return un.f;
}

#include <assert.h>
#include <math.h>
void test(vector float a) {
  vector float b = _ZGV9N4v_logf(a);
for (int i=0;i<4;i++) {
if (isnan(b[i]))
  assert(isnan(logf(a[i])));
else
  assert(logf(a[i]) == b[i]);
}
}

int main() {
 vector float a = {0, 1, 2, 3};
 test(a);
 float r = NAN;
 vector float b = {INFINITY, -INFINITY, -0.0,0.0};
 test(b);

}
