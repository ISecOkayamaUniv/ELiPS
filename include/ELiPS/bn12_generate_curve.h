#ifndef BN12_GENERATE_CURVE_H
#define BN12_GENERATE_CURVE_H

#include <ELiPS/mpn.h>

extern void bn12_generate_X();
extern int  bn12_generate_prime();
extern int  bn12_generate_order();
extern void bn12_generate_trace();
extern void bn12_weil();
#endif
