#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) * 1000.0f + (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_usec - tv_start.tv_usec);
}
int main(){
    struct timeval tv_A,tv_B;
    float time;
    int n=100,j=0;
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++){
        //loop
    }
    gettimeofday(&tv_B,NULL);
    time+=timedifference_msec(tv_A,tv_B)/n;
}