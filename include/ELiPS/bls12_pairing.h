#ifndef BLS12_PAIRING_H
#define BLS12_PAIRING_H

#include <ELiPS/bls12_miller.h>
#include <ELiPS/bls12_final_exp.h>
#include <ELiPS/time.h>

/**
 * @brief Calculation ate pairing on prime field for bls12
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_ate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for bls12
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_optate_pairing_basic(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for bls12
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_optate_pairing_compress(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for bls12 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_optate_pairing_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for bls12
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_optate_pairing_compress_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bls12_optate_pairing_projective_compress_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bls12_optate_pairing_compress_lazy_montgomery(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bls12_optate_pairing_projective_compress_lazy_montgomery(fp12_t *ANS,efp12_t *P,efp12_t *Q);


#endif
