vector double _ZGVbN2v_log(vector double);
vector float _ZGVbN4v_logf(vector float);
vector float _ZGVbN4v_expf(vector float);
vector float _ZGVbN4vv_powf(vector float, vector float);
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
vector float ind = {79, 79, 79, 1.17549435082228750797e-38};
         vector float agdain = _ZGVbN4v_expf(ind);
for (int i=0;i<4;i++) {
assert(agdain[i] == expf(ind[i]));
}
  uint64_t ret;
  uint32_t *data = malloc(SIZE * 6);
  volatile float *bench = (float *)((char*)data + SIZE*2);
  volatile float *bench2 = (float *)((char*)bench + SIZE);
  volatile float *dataconv = (float *)((char*)bench2 + SIZE);
  assert(data && bench && bench2 && dataconv);
  FILE *random = fopen("/dev/urandom", "r");
  if (!random)
    return EXIT_UNSUPPORTED;
  ret = fread(data, 6, SIZE/4, random);
  if (ret != SIZE/4)
    return 1;
  struct timeval a, b, c;
  for (int i=0;i<SIZE/4;i+=4) {
#if AVOID_SPECIAL_NUMBERS
    data[i] &= 0xffff;
    data[i+1] &= 0xffff;
dataconv[i] = (float)data[i] + 1;
dataconv[i+1] = (float)data[i+1] + 1;
    data[i+2] &= 0xffff;
    data[i+3] &= 0xffff;
dataconv[i+2] = (float)data[i+2] + 1;
dataconv[i+3] = (float)data[i+3] + 1;
    data[i+4] &= 0xffff;
    data[i+5] &= 0xffff;
dataconv[i+4] = (float)data[i+4] + 1;
dataconv[i+5] = (float)data[i+5] + 1;
    data[i+6] &= 0xffff;
    data[i+7] &= 0xffff;
dataconv[i+6] = (float)data[i+6] + 1;
dataconv[i+7] = (float)data[i+7] + 1;
#else
us d;
d.u = data[i];
dataconv[i] = d.u;
d.u = data[i+1];
dataconv[i+1] = d.u;
d.u = data[i+2];
dataconv[i+2] = d.u;
d.u = data[i+3];
dataconv[i+3] = d.u;
d.u = data[i+4];
dataconv[i+4] = d.u;
d.u = data[i+5];
dataconv[i+5] = d.u;
d.u = data[i+6];
dataconv[i+6] = d.u;
d.u = data[i+7];
dataconv[i+7] = d.u;
#endif
    vector float in = vec_ld(0, (vector float*)&dataconv[i]);
    vector float in2 = vec_ld(0, (vector float*)&dataconv[i+4]);
    vector float four = _ZGVbN4vv_powf(in, in2);
    for (int j=0;j<4;j++) {
      four[j] = isnan(four[j]) ? NAN : four[j];
us u;
u.f = four[j];
us u2;
u2.f = powf((float)dataconv[i+j], (float)dataconv[i+j+4]);
u2.f = isnan(u2.f) ? NAN : u2.f;
      if (u.u != u2.u) {
        printf("round %u: %.50e -> %.50e and %.50e not equal!\n", i, in[j], four[j], powf(dataconv[i+j]));
        volatile vector float again = _ZGVbN4vv_powf(in);
        return 1;
      }
    }
  }

  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, (vector float*)&dataconv[i]);
    vector float res = _ZGVbN4v_expf(in);
    *(vector float*)&bench[i] = res;
  }
  for (int i=0;i<SIZE/4;i+=1) {
    bench2[i] = expf(dataconv[i]);
  }
uint64_t table[20];
uint64_t tableopt[20];
for (int j=0;j<20;j++) {
struct timespec q, w, e;
clock_gettime(CLOCK_MONOTONIC, &q);
  for (int i=0;i<SIZE/4;i+=4) {
    vector float in = vec_ld(0, (vector float*)&dataconv[i]);
    vector float in2 = vec_ld(0, (vector float*)&dataconv[i+4]);
    vector float res = _ZGVbN4vv_expf(in, in2);
    *(vector float*)&bench[i] = res;
  }
clock_gettime(CLOCK_MONOTONIC, &w);
  for (int i=0;i<SIZE/4;i+=2) {
    bench2[i] = powf(dataconv[i], dataconv[i+4]);
    bench2[i+1] = powf(dataconv[i+1], dataconv[i+5]);
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
