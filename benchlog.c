vector double _ZGVbN2v_log(vector double);
typedef union {
float f;
double d;
long long unsigned l;
unsigned u;
} us;
#define AVOID_SPECIAL_NUMBERS 1
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
  uint64_t *data = malloc(SIZE * 4);
  double *bench = (double *)((char*)data + SIZE);
  double *bench2 = (double *)((char*)bench + SIZE);
  double *dataconv = (double *)((char*)bench2 + SIZE);
  assert(data && bench && bench2 && dataconv);
  FILE *random = fopen("/dev/urandom", "r");
  if (!random)
    return EXIT_UNSUPPORTED;
  ret = fread(data, 4, SIZE/4, random);
  if (ret != SIZE/4)
    return 1;
  struct timeval a, b, c;
  for (int i=0;i<SIZE/8;i+=2) {
#if AVOID_SPECIAL_NUMBERS
    data[i] &= 0xffffffff;
    data[i+1] &= 0xffffffff;
dataconv[i] = (double)data[i] + 1;
dataconv[i+1] = (double)data[i+1] + 1;
#else
us d;
d.l = data[i];
dataconv[i] = d.d;
d.l = data[i+1];
dataconv[i+1] = d.d;
#endif
    vector double in = vec_ld(0, &dataconv[i]);
    vector double four = _ZGVbN2v_log(in);
    for (int j=0;j<2;j++) {
      four[j] = isnan(four[j]) ? NAN : four[j];
us u;
u.d = four[j];
us u2;
u2.d = log((double)dataconv[i+j]);
u2.d = isnan(u2.d) ? NAN : u2.d;
us u3;
      if (u.l != u2.l) {
        printf("round %u: %.50e -> %.50e and %.50e not equal!\n", i, in[j], four[j], log(dataconv[i+j]));
        volatile vector double again = _ZGVbN2v_log(in);
        return 1;
      }
    }
  }

  for (int i=0;i<SIZE/8;i+=4) {
    vector double in = vec_ld(0, &dataconv[i]);
    vector double res = _ZGVbN2v_log(in);
    *(vector double*)&bench[i] = res;
  }
  for (int i=0;i<SIZE/8;i+=1) {
    bench2[i] = log(dataconv[i]);
  }
struct timespec q, w, e;
clock_gettime(CLOCK_MONOTONIC, &q);
  for (int i=0;i<SIZE/8;i+=4) {
    vector double in = vec_ld(0, &dataconv[i]);
    vector double res = _ZGVbN2v_log(in);
    *(vector double*)&bench[i] = res;
  }
clock_gettime(CLOCK_MONOTONIC, &w);
  for (int i=0;i<SIZE/8;i+=1) {
    bench2[i] = log(dataconv[i]);
  }
clock_gettime(CLOCK_MONOTONIC, &e);
  double t[3] = {(double)(q.tv_sec * 1000000) + q.tv_nsec / 1000, (double)(w.tv_sec * 1000000) + w.tv_nsec / 1000, (double)(e.tv_sec * 1000000) + e.tv_nsec / 1000};
  printf("non-opt: %f (%f MiB/s)\nopt: %f (%f MiB/s) (%f)\n",
    (t[2] - t[1]) / 1000000, SIZEM / ((t[2] - t[1]) / 1000000),
    (t[1] - t[0]) / 1000000, SIZEM / ((t[1] - t[0]) / 1000000),
    (t[2] - t[1]) / (t[1] - t[0]));
  return 0;
}
