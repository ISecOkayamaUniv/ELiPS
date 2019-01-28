#ifndef TIME_H
#define TIME_H

#include <ELiPS/Define.h>

float timedifference_msec(struct timeval tv_start, struct timeval tv_end);
float timedifference_usec(struct timeval tv_start, struct timeval tv_end);
#endif
