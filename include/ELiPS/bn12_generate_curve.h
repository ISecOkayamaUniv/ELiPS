#ifndef BN12_GENERATE_CURVE_H
#define BN12_GENERATE_CURVE_H

#include <ELiPS/mpn.h>

void BN12_generate_X();
int  BN12_generate_prime();
int  BN12_generate_order();
void BN12_generate_trace();
void BN12_weil();
#endif
