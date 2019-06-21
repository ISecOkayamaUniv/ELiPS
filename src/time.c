#include <ELiPS/time.h>
/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) * 1000.0f + (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_usec - tv_start.tv_usec);
}

float timedifference_nsec(struct timespec tv_start, struct timespec tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_nsec - tv_start.tv_nsec);
}
