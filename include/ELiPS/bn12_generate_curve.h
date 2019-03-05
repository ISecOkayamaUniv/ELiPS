#ifndef BN12_GENERATE_CURVE_H
#define BN12_GENERATE_CURVE_H

#include <ELiPS/mpn.h>

extern void BN12_generate_X();
extern int  BN12_generate_prime();
extern int  BN12_generate_order();
extern void BN12_generate_trace();
extern void BN12_weil();
#endif
