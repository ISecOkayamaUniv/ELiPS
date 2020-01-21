#ifndef BN12_FINAL_EXP_H
#define BN12_FINAL_EXP_H

#include <ELiPS/fp12.h>
#include <ELiPS/test.h>

//final exp
extern void bn12_final_exp_plain(fp12_t *ANS,fp12_t *A);
extern void bn12_fp12_pow_X(fp12_t *ANS,fp12_t *A);
extern void bn12_fp12_pow_X_lazy(fp12_t *ANS,fp12_t *A);
extern void bn12_final_exp_optimal(fp12_t *ANS,fp12_t *A);
extern void bn12_final_exp_optimal_lazy(fp12_t *ANS,fp12_t *A);

#endif
