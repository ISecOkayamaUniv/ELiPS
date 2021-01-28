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
extern void bls12_optate_pairing_basic(fp12_t *ANS,efp12_t *P,efp12_t *Q);
/**
 * @brief Calculation Optimal-ate pairing on prime field for bls12
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_optate_pairing_affine(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for bls12
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_optate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q);

/**
 * @brief Calculation Optimal-ate pairing on prime field for bls12(symmetric)
 *
 * @param[out]ANS --a pointer of answer in fp12_t.
 * @param[in]P --a pointer in efp12_t on G1.
 * @param[in]Q --a pointer in efp12_t on G2.
 */
extern void bls12_symmetric_optate_pairing(fp12_t *ANS,sym_t *A,sym_t *B);

extern void bls12_optate_pairing_basic_cvma(fp12_t *ANS,efp12_t *P,efp12_t *Q);
#endif
