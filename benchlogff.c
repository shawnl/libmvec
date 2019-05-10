vector double _ZGVbN2v_log(vector double);
vector float _ZGVbN4v_logf(vector float);

typedef union {
float f;
double d;
long long unsigned l;
unsigned u;
} us;
#define AVOID_SPECIAL_NUMBERS 0
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
  volatile float *bench = (float *)((char*)data + SIZE);
  volatile float *bench2 = (float *)((char*)bench + SIZE);
  volatile float *dataconv = (float *)((char*)bench2 + SIZE);
  assert(data && bench && bench2 && dataconv);
  FILE *random = fopen("/dev/urandom", "r");
  if (!random)
    return EXIT_UNSUPPORTED;
  ret = fread(data, 4, SIZE/4, random);
  if (ret != SIZE/4)
    return 1;
  struct timeval a, b, c;
  for (int i=0;i<SIZE/4;i+=4) {
#if AVOID_SPECIAL_NUMBERS
    data[i] &= 0xffff;
    data[i+1] &= 0xffff;
dataconv[i] = (float)data[i] + 1;
dataconv[i+1] = (float)data[i+1] + 1;
#else
us d;
d.u = data[i];
dataconv[i] = d.u;
d.u = data[i+1];
dataconv[i+1] = d.u;
#endif
    vector float in = vec_ld(0, &dataconv[i]);
    vector float four = _ZGVbN4v_logf(in);
    for (int j=0;j<4;j++) {
      four[j] = isnan(four[j]) ? NAN : four[j];
us u;
u.f = four[j];
us u2;
u2.f = logf((float)dataconv[i+j]);
u2.f = isnan(u2.f) ? NAN : u2.f;
us u3;
      if (u.u != u2.u) {
        printf("round %u: %.50e -> %.50e and %.50e not equal!\n", i, in[j], four[j], log(dataconv[i+j]));
        volatile vector float again = _ZGVbN4v_logf(in);
        return 1;
      }
    }
  }

  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector float res = _ZGVbN4v_logf(in);
    *(vector float*)&bench[i] = res;
  }
  for (int i=0;i<SIZE/4;i+=1) {
    bench2[i] = logf(dataconv[i]);
  }
uint64_t table[20];
uint64_t tableopt[20];
for (int j=0;j<20;j++) {
struct timespec q, w, e;
clock_gettime(CLOCK_MONOTONIC, &q);
  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, &dataconv[i]);
    vector float res = _ZGVbN4v_logf(in);
    *(vector float*)&bench[i] = res;
  }
clock_gettime(CLOCK_MONOTONIC, &w);
  for (int i=0;i<SIZE/4;i+=1) {
    bench2[i] = logf(dataconv[i]);
  }
clock_gettime(CLOCK_MONOTONIC, &e);
  double t[3] = {(double)(q.tv_sec * 1000000) + q.tv_nsec / 1000, (double)(w.tv_sec * 1000000) + w.tv_nsec / 1000, (double)(e.tv_sec * 1000000) + e.tv_nsec / 1000};
  printf("non-opt: %f (%f MiB/s)\nopt: %f (%f MiB/s) (%f)\n",
    (t[2] - t[1]) / 1000000, SIZEM / ((t[2] - t[1]) / 1000000),
    (t[1] - t[0]) / 1000000, SIZEM / ((t[1] - t[0]) / 1000000),
    (t[2] - t[1]) / (t[1] - t[0]));
table[j] = t[2] - t[1];
tableopt[j] = t[1] - t[0];
}
uint64_t sum = 0;
uint64_t sumopt = 0;
for (int j=0;j<20;j++) {
  sum += table[j];
  sumopt += tableopt[j];
}
double mean = sum / 19;
double meanopt = sumopt / 19;
double sd = 0.0;
double sdopt = 0.0;
for (int j=0;j<20;j++) {
  sd +=( (double)(table[j] - mean) * (double)(table[j] - mean)) / 1000000000000.0;
  sdopt +=( (double)(tableopt[j] - mean) * (double)(tableopt[j] - mean)) / 1000000000000.0;
}
printf("opt: mean %f (sd %f)\nnonopt: mean %f (sd %f)\n", SIZEM / (meanopt / 1000000), sqrt(sdopt), SIZEM / (mean / 1000000), sqrt(sd));
  return 0;
}
