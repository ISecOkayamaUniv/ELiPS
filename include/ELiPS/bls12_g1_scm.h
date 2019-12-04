#ifndef BLS12_G1_SCM_H
#define BLS12_G1_SCM_H

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
extern void bls12_g1_scm_basic(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + naf)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_2naf(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + naf)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_3naf_shamia(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + naf)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_3naf_interleaving(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + naf)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_5naf_interleaving(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + jsf)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_jsf(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + jsf + jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_jsf_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + jsf + jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_jsf_jacobian_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + jsf + jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_jsf_mixture_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + 5naf_interleaving + jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_5naf_interleaving_mixture(efp12_t *ANS,efp12_t *P,mpz_t scalar);
/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + 5naf_interleaving + jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_5naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);
extern void bls12_g1_scm_2split_5naf_interleaving_mixture_lazy_montgomery(efp12_t *ANS,efp12_t *P,mpz_t scalar);
/**
 * @brief Scalar multiplication a efp12_t type struct on G1 for bls12 (GLV-2split + 5naf_interleaving + jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g1_scm_2split_7naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);

#endif
