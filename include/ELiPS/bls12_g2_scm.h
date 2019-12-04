#ifndef BLS12_G2_SCM_H
#define BLS12_G2_SCM_H

#include <ELiPS/twist.h>
#include <ELiPS/jsf_naf.h>
#include <ELiPS/time.h>


/**
 * @brief Scalar multiplication a efp12_t type struct on G2
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_basic(efp12_t *ANS,efp12_t *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G2 (Lazy Redution)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G2 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_2split(efp12_t *ANS,efp12_t *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G2 (GLV-2split + jsf)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_2split_jsf(efp12_t *ANS,efp12_t *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G2 for bls12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_4split(efp12_t *ANS,efp12_t *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G2 for bls12 (GLV-4split + Lazy)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_4split_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G2 for bls12 (GLV-4split + jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_4split_jacobian_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar);


/**
 * @brief Scalar multiplication a efp12_t type struct on G2 for bls12 (GLV-4split + jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g2_scm_4split_mixture_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar);


extern void bls12_g2_scm_4split_3naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar);
extern void bls12_g2_scm_4split_5naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar);
extern void bls12_g2_scm_4split_5naf_interleaving_mixture_lazy_montgomery(efp12_t *ANS,efp12_t *Q,mpz_t scalar);
#endif
