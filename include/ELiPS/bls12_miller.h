#ifndef BLS12_MILLER_H
#define BLS12_MILLER_H

#include <ELiPS/line_calculate.h>
#include <ELiPS/twist.h>

//miller

/**
 * @brief Calculation Miller's Alg on prime field for ate pairing
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_miller_algo_for_plain_ate(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Miller's Alg on prime field for Optimal-ate pairing
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
 
 //TODO: opt_ate->optate
extern void bls12_miller_algo_for_opt_ate(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Miller's Alg on prime field for Optimal-ate pairing (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_miller_algo_for_opt_ate_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bls12_miller_algo_for_opt_ate_projective_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bls12_miller_algo_for_opt_ate_lazy_montgomery(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bls12_miller_algo_for_opt_ate_projective_lazy_montgomery(fp12_t *ANS,efp12_t *P,efp12_t *Q);


#endif
