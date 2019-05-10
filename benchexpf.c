vector double _ZGVN2v_log(vector double);
vector float _ZGVN4v_logf(vector float);
vector float _ZGVbN4v_expf(vector float);
typedef union {
float f;
double d;
long long unsigned l;
unsigned u;
} us;
#include <altivec.h>
#include <math.h>
#define EXIT_UNSUPPORTED 77
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#define SIZEM 32
#define SIZE SIZEM * 1024 * 1024
int main() {
  uint64_t ret;
  uint32_t *data = malloc(SIZE * 4);
  float *bench = (float *)((char*)data + SIZE);
  float *bench2 = (float *)((char*)bench + SIZE);
  float *dataconv = (float *)((char*)bench2 + SIZE);
  assert(data && bench && bench2);
  FILE *random = fopen("/dev/urandom", "r");
  if (!random)
    return EXIT_UNSUPPORTED;
  ret = fread(data, 4, SIZE/4, random);
  if (ret != SIZE/4)
    return 1;
  for (int i=0;i<SIZE/4;i+=4) {
    data[i] &= 0xffff;
    data[i+1] &= 0xffff;
    data[i+2] &= 0xffff;
    data[i+3] &= 0xffff;
    vector float in = {(float)data[i], (float)data[i+1], (float)data[i+2], (float)data[i+3]};
    vector float four = _ZGVN4v_expf(in);
    for (int j=0;j<4;j++) {
      four[j] = isnan(four[j]) ? NAN : four[j];
us u;
u.f = four[j];
us u2;
u2.f = expf(data[i+j]);
u2.f = isnan(u2.f) ? NAN : u2.f;
us u3;
u3.u = data[i+j];
dataconv[i+j] = u3.f;
      if (u.u != u2.u) {
        printf("round %u: %.15e -> %.15e and %.15e not equal!\n %.15e\n", i, in[j], four[j], logf(data[i+j]), log(data[i+j]));
        volatile vector float again = _ZGVN4v_expf(in);
        return 1;
      }
    }
  }

  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector float res = _ZGVN4v_expf(in);
    *(vector float*)&bench[i] = res;
  }
  for (int i=0;i<SIZE/4;i+=1) {
    bench2[i] = expf(dataconv[i]);
  }
struct timespec a, b, c;
clock_gettime(CLOCK_MONOTONIC, &a);
  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector float res = _ZGVN4v_expf(in);
    *(vector float*)&bench[i] = res;
  }
clock_gettime(CLOCK_MONOTONIC, &b);
  for (int i=0;i<SIZE/4;i+=1) {
    bench2[i] = expf(dataconv[i]);
  }
clock_gettime(CLOCK_MONOTONIC, &c);

  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector double in1 = vec_unpackh(in);
    vector double in2 = vec_unpackl(in);
    vector double res1 = _ZGVN2v_exp(in1);
    vector double res2 = _ZGVN2v_log(in2);
    vector float res = vec_pack(res1, res2);
    *(vector float*)&bench[i] = res;
  }
//clock_gettime(CLOCK_MONOTONIC, &c);
  double t[3] = {(double)(a.tv_sec * 1000000) + a.tv_nsec / 1000, (double)(b.tv_sec * 1000000) + b.tv_nsec / 1000, (double)(c.tv_sec * 1000000) + c.tv_nsec / 1000};
  printf("non-opt: %f (%f MiB/s)\nopt: %f (%f MiB/s) (%f)\n",
    (t[2] - t[1]) / 1000000, SIZEM / ((t[2] - t[1]) / 1000000),
    (t[1] - t[0]) / 1000000, SIZEM / ((t[1] - t[0]) / 1000000), (t[2] - t[1]) / (t[1] - t[0]));
  return 0;
}
