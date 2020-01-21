#ifndef BLS12_FINAL_EXP_H
#define BLS12_FINAL_EXP_H

#include <ELiPS/fp12.h>
#include <ELiPS/time.h>


//final exp

/**
 * @brief Calculation final exponentiation on prime field for ate pairing
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]A --a pointer in fp12_t.
 */


extern void bls12_fp12_pow_X(fp12_t *ANS,fp12_t *A);
extern void bls12_fp12_pow_X2(fp12_t *ANS,fp12_t *A);

extern void bls12_fp12_pow_X_compress(fp12_t *ANS,fp12_t *A);
extern void bls12_fp12_pow_X2_compress(fp12_t *ANS,fp12_t *A);

/**
 * @brief Calculation final exponentiation on prime field for Optimal-ate pairing
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]A --a pointer in fp12_t.
 */
extern void bls12_final_exp_optimal(fp12_t *ANS,fp12_t *A);
/**
 * @brief Calculation final exponentiation on prime field for Optimal-ate pairing
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]A --a pointer in fp12_t.
 */
extern void bls12_final_exp_optimal_compress(fp12_t *ANS,fp12_t *A);
#endif
