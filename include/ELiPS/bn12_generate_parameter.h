#ifndef BN12_GENERATE_PARAMETER_H
#define BN12_GENERATE_PARAMETER_H

#include <ELiPS/Fp.h>
#include <ELiPS/Fp2.h>
#include <ELiPS/bn12_generate_curve.h>

extern void BN12_get_epsilon();
extern void BN12_get_Two_inv();
extern void BN12_set_basis();
extern void BN12_set_frobenius_constant();
extern void BN12_set_curve_parameter();

#endif
