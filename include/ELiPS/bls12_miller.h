#ifndef BLs12_MILLER_H
#define BLs12_MILLER_H

#include <ELiPS/sparse.h>
#include <ELiPS/twist.h>

//miller

/**
 * @brief Calculation Miller's Alg on prime field for ate pairing
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
void BLS12_Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * @brief Calculation Miller's Alg on prime field for Optimal-ate pairing
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
void BLS12_Miller_algo_for_opt_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * @brief Calculation Miller's Alg on prime field for Optimal-ate pairing (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
void BLS12_Miller_algo_for_opt_ate_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * @brief Calculation Miller's Alg on prime field for X-ate pairing
 *
 * @param[out]ANS --a pointer of answer in Fp12.
 * @param[in]P --a pointer in EFp12 on G1.
 * @param[in]Q --a pointer in EFp12 on G2.
 */
void BLS12_Miller_algo_for_x_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q);

#endif
