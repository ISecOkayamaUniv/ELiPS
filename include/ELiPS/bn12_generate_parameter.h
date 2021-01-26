#ifndef BN12_GENERATE_PARAMETER_H
#define BN12_GENERATE_PARAMETER_H

#include <ELiPS/fp.h>
#include <ELiPS/fp2.h>
#include <ELiPS/bn12_generate_curve.h>

extern void bn12_get_epsilon();
extern void bn12_get_Two_inv();
extern void bn12_set_basis();
extern void bn12_set_frobenius_constant();
extern void bn12_set_curve_parameter();

#endif
