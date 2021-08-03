#ifndef TIME_H
#define TIME_H

#include <ELiPS/define.h>

extern float timedifference_msec(struct timeval tv_start, struct timeval tv_end);
extern float timedifference_usec(struct timeval tv_start, struct timeval tv_end);
//extern float timedifference_nsec(struct timespec tv_start, struct timespec tv_end);
extern void cost_zero();
extern void cost_init(cost *A);
extern void cost_check(cost *A);
extern void cost_addition(cost *A, cost *B);
extern void cost_substruction(cost *ANS, cost *A, cost *B);
extern void cost_printf(char *str, cost *A, int n);
#endif
