#ifndef TIME_H
#define TIME_H

#include <ELiPS/Define.h>

extern float timedifference_msec(struct timeval tv_start, struct timeval tv_end);
extern float timedifference_usec(struct timeval tv_start, struct timeval tv_end);
#endif
