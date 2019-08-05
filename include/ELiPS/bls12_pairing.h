#ifndef BLS12_PAIRING_H
#define BLS12_PAIRING_H

#include <ELiPS/bls12_miller.h>
#include <ELiPS/bls12_final_exp.h>
#include <ELiPS/time.h>

/**
 * @brief Calculation ate pairing on prime field for BLS12
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
extern void BLS12_Plain_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for BLS12
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
extern void BLS12_Opt_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for BLS12
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
extern void BLS12_Opt_ate_pairing_compress(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for BLS12 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
extern void BLS12_Opt_ate_pairing_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for BLS12
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
extern void BLS12_Opt_ate_pairing_compress_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q);
extern void BLS12_Opt_ate_pairing_compress_lazy_montgomery(Fp12 *ANS,EFp12 *P,EFp12 *Q);


#endif
