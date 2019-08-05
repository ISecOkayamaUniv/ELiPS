#ifndef BLS12_FINAL_EXP_H
#define BLS12_FINAL_EXP_H

#include <ELiPS/Fp12.h>
#include <ELiPS/time.h>


//final exp

/**
 * @brief Calculation final exponentiation on prime field for ate pairing
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]A --a pointer in Fp12.
 */
extern void BLS12_Final_exp_plain(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X_compress(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X_compress2(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X_lazy(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X_compress_lazy(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X_compress_lazy_montgomery(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X2(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X2_compress(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X2_lazy(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X2_compress_lazy(Fp12 *ANS,Fp12 *A);
extern void BLS12_Fp12_pow_X2_compress_lazy_montgomery(Fp12 *ANS,Fp12 *A);

/**
 * @brief Calculation final exponentiation on prime field for Optimal-ate pairing
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]A --a pointer in Fp12.
 */
extern void BLS12_Final_exp_optimal(Fp12 *ANS,Fp12 *A);

/**
 * @brief Calculation final exponentiation on prime field for Optimal-ate pairing
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]A --a pointer in Fp12.
 */
extern void BLS12_Final_exp_optimal_compress(Fp12 *ANS,Fp12 *A);

/**
 * @brief Calculation final exponentiation on prime field for Optimal-ate pairing (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]A --a pointer in Fp12.
 */
extern void BLS12_Final_exp_optimal_lazy(Fp12 *ANS,Fp12 *A);

/**
 * @brief Calculation final exponentiation on prime field for Optimal-ate pairing
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]A --a pointer in Fp12.
 */
extern void BLS12_Final_exp_optimal_compress_lazy(Fp12 *ANS,Fp12 *A);

extern void BLS12_Final_exp_optimal_compress_lazy_montgomery(Fp12 *ANS,Fp12 *A);
#endif
