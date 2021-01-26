#ifndef BLS12_G2_SCM_H
#define BLS12_G2_SCM_H

#include <ELiPS/twist.h>
#include <ELiPS/jsf_naf.h>
#include <ELiPS/time.h>


/**
 * @brief Scalar multiplication a efp12_t type struct on G2 (basic type)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_basic(efp12_t *ANS,efp12_t *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G2 (4split_5naf_interleaving_mixture)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm(efp12_t *ANS,efp12_t *Q,mpz_t scalar);
#endif
