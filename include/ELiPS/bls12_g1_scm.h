#ifndef BLS12_G1_SCM_H
#define BLS12_G1_SCM_H

#include <ELiPS/jsf_naf.h>
#include <ELiPS/time.h>
#include <ELiPS/twist.h>

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 (basic type)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_basic(efp12_t *ANS, efp12_t *P, mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 (2split_5naf_interleaving_mixture)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm(efp12_t *ANS, efp12_t *P, mpz_t scalar);

#endif
