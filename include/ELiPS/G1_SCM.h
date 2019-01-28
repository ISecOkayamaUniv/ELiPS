#ifndef G1_SCM_H
#define G1_SCM_H

#include <ELiPS/twist.h>
#include <ELiPS/JSF.h>
#include <ELiPS/test.h>


/**
 * @brief Scalar multiplication a EFp12 type struct on G1
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void EFp12_G1_SCM_plain(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void EFp12_G1_SCM_plain_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 for BN12 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void BN12_EFp12_G1_SCM_2split(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 for BN12 (GLV-2split + JSF)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void BN12_EFp12_G1_SCM_2split_JSF(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 for BN12 (GLV-2split + JSF +Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void BN12_EFp12_G1_SCM_2split_JSF_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 for BLS12 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void BLS12_EFp12_G1_SCM_2split(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 for BLS12 (GLV-2split + JSF)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void BLS12_EFp12_G1_SCM_2split_JSF(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 for BLS12 (GLV-2split + JSF + Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void BLS12_EFp12_G1_SCM_2split_JSF_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct on G1 for BLS12 (GLV-2split + JSF + Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void BLS12_EFp12_G1_SCM_2split_JSF_Jacobian_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar);

#endif
