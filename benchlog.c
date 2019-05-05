#include "svml_s_log2.c"
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
    data[i] &= 0xffffffff;
    data[i+1] &= 0xffffffff;
dataconv[i] = (double)data[i];
dataconv[i+1] = (double)data[i+1];
    vector double in = vec_ld(0, &dataconv[i]);
    vector double four = _ZGV9N2v_log(in);
    for (int j=0;j<2;j++) {
      four[j] = isnan(four[j]) ? NAN : four[j];
us u;
u.d = four[j];
us u2;
u2.d = log((double)data[i+j]);
u2.d = isnan(u2.d) ? NAN : u2.d;
us u3;
// vector routine more accurate due to use of fused multiply-add
//      if (u.l != u2.l) {
//        printf("round %u: %.25e -> %.25e and %.25e not equal!\n", i, in[j], four[j], log(data[i+j]));
//        volatile vector double again = _ZGV9N2v_log(in);
//        return 1;
//      }
    }
  }

  gettimeofday(&a, NULL);
  for (int i=0;i<SIZE/8;i+=4) {
    vector double in = vec_ld(0, &dataconv[i]);
    vector double res = _ZGV9N2v_log(in);
    *(vector double*)&bench[i] = res;
  }
  gettimeofday(&b, NULL);
  for (int i=0;i<SIZE/8;i+=1) {
    bench2[i] = log(dataconv[i]);
  }
  gettimeofday(&a, NULL);
  for (int i=0;i<SIZE/8;i+=4) {
    vector double in = vec_ld(0, &dataconv[i]);
    vector double res = _ZGV9N2v_log(in);
    *(vector double*)&bench[i] = res;
  }
  gettimeofday(&b, NULL);
  for (int i=0;i<SIZE/8;i+=1) {
    bench2[i] = log(dataconv[i]);
  }
  gettimeofday(&c, NULL);
  double t[3] = {a.tv_sec * 1000000 + a.tv_usec, b.tv_sec * 1000000 + b.tv_usec, c.tv_sec * 1000000 + c.tv_usec};
  printf("non-opt: %f (%f MiB/s)\nopt: %f (%f MiB/s) (%f)\n",
    (t[2] - t[1]) / 1000000, SIZEM / ((t[2] - t[1]) / 1000000),
    (t[1] - t[0]) / 1000000, SIZEM / ((t[1] - t[0]) / 1000000),
    (t[2] - t[1]) / (t[1] - t[0]));
  return 0;
}
