vector double _ZGV9N2v_log(vector double);
vector float _ZGV9N4v_logf(vector float);
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
#include <sys/time.h>
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
  struct timeval a, b, c, e, f;
  for (int i=0;i<SIZE/4;i+=4) {
    data[i] &= 0xffff;
    data[i+1] &= 0xffff;
    data[i+2] &= 0xffff;
    data[i+3] &= 0xffff;
    vector float in = {(float)data[i], (float)data[i+1], (float)data[i+2], (float)data[i+3]};
    vector float four = _ZGV9N4v_logf(in);
    for (int j=0;j<4;j++) {
      four[j] = isnan(four[j]) ? NAN : four[j];
us u;
u.f = four[j];
us u2;
u2.f = logf(data[i+j]);
u2.f = isnan(u2.f) ? NAN : u2.f;
us u3;
u3.u = data[i+j];
dataconv[i+j] = u3.f;
      if (u.u != u2.u) {
        printf("round %u: %.15e -> %.15e and %.15e not equal!\n %.15e\n", i, in[j], four[j], logf(data[i+j]), log(data[i+j]));
        volatile vector float again = _ZGV9N4v_logf(in);
        return 1;
      }
    }
  }

  gettimeofday(&a, NULL);
  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector float res = _ZGV9N4v_logf(in);
    *(vector float*)&bench[i] = res;
  }
  gettimeofday(&b, NULL);
  for (int i=0;i<SIZE/4;i+=1) {
    bench2[i] = logf(dataconv[i]);
  }
  gettimeofday(&a, NULL);
  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector float res = _ZGV9N4v_logf(in);
    *(vector float*)&bench[i] = res;
  }
  gettimeofday(&b, NULL);
  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector double in1 = vec_unpackh(in);
    vector double in2 = vec_unpackl(in);
    vector double res1 = _ZGV9N2v_log(in1);
    vector double res2 = _ZGV9N2v_log(in2);
    vector float res = vec_pack(res1, res2);
    *(vector float*)&bench[i] = res;
  }
  gettimeofday(&c, NULL);
  double t[3] = {a.tv_sec * 1000000 + a.tv_usec, b.tv_sec * 1000000 + b.tv_usec, c.tv_sec * 1000000 + c.tv_usec};
  printf("non-opt: %f (%f MiB/s)\nopt: %f (%f MiB/s) (%f)\n",
    (t[2] - t[1]) / 1000000, SIZEM / ((t[2] - t[1]) / 1000000),
    (t[1] - t[0]) / 1000000, SIZEM / ((t[1] - t[0]) / 1000000), (t[2] - t[1]) / (t[1] - t[0]));
  return 0;
}
