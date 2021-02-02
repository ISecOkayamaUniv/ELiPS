#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"

struct timeval tv_s,tv_e;
void bench_start(){
  gettimeofday(&tv_s,NULL);
}
float bench_end(float time){
  gettimeofday(&tv_e,NULL);
  float tmp=0;
  tmp+=timedifference_msec(tv_s,tv_e);
  return tmp;
}
float bench_result(float time,int cnt){
  return time/cnt;
}
