#ifndef G2_SCM_H
#define G2_SCM_H

#include <ELiPS/twist.h>
#include <ELiPS/JSF.h>
#include <ELiPS/time.h>

/**
 * @brief Scalar multiplication a EFp12 type struct on G2
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp12_G2_SCM_plain(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 (Lazy Redution)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp12_G2_SCM_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp12_G2_SCM_2split(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 (GLV-2split + JSF)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp12_G2_SCM_2split_JSF(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 for BN12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BN12_EFp12_G2_SCM_4split(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 for BN12 (GLV-4split + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BN12_EFp12_G2_SCM_4split_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 for BLS12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_EFp12_G2_SCM_4split(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 for BLS12 (GLV-4split + Lazy)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_EFp12_G2_SCM_4split_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 for BLS12 (GLV-4split + Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_EFp12_G2_SCM_4split_Jacobian_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G2 for BLS12 (GLV-4split + Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_EFp12_G2_SCM_4split_Jacobian_table(EFp12 *ANS,EFp12 *Q,mpz_t scalar);


/**
 * @brief Scalar multiplication a EFp12 type struct on G2 for BLS12 (GLV-4split + Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_EFp12_G2_SCM_4split_Mixture_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar);


extern void BLS12_EFp12_G2_SCM_4split_3NAF_interleaving_Mixture_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
extern void BLS12_EFp12_G2_SCM_4split_5NAF_interleaving_Mixture_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
extern void BLS12_EFp12_G2_SCM_4split_5NAF_interleaving_Mixture_lazy_montgomery(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
#endif
