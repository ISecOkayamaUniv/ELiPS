#ifndef BN12_G1_SCM_H
#define BN12_G1_SCM_H

#include <ELiPS/twist.h>
#include <ELiPS/jsf_naf.h>
#include <ELiPS/time.h>


/**
 * @brief Scalar multiplication a efp12_t type struct on G1
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g1_scm_plain(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g1_scm_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for BN12 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g1_scm_2split(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for BN12 (GLV-2split + JSF)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g1_scm_2split_JSF(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for BN12 (GLV-2split + JSF +Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g1_scm_2split_JSF_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);
#endif
