#ifndef BN12_FINAL_EXP_H
#define BN12_FINAL_EXP_H

#include <ELiPS/Fp12.h>
#include <ELiPS/test.h>

//final exp
void BN12_Final_exp_plain(Fp12 *ANS,Fp12 *A);
void BN12_Fp12_pow_X(Fp12 *ANS,Fp12 *A);
void BN12_Fp12_pow_X_lazy(Fp12 *ANS,Fp12 *A);
void BN12_Final_exp_optimal(Fp12 *ANS,Fp12 *A);
void BN12_Final_exp_optimal_lazy(Fp12 *ANS,Fp12 *A);

#endif
