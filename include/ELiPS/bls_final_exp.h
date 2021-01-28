#ifndef INCLUDED_bls_final_exp_h
#define INCLUDED_bls_final_exp_h

#include <ELiPS/fpm2.h>

extern void bls_final_exp_optimal(fpm2_t *ANS, fpm2_t *a);
extern void bls_final_exp(fpm2_t *ANS, fpm2_t *a);
extern int Eulers_totient_function(int n);
extern int gcd(int a, int b);
extern void bls_fpm2_pow_X(fpm2_t *ANS, fpm2_t *a);
extern void bls_fpm2_pow_X2(fpm2_t *ANS, fpm2_t *a);
#endif