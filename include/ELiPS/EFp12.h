#ifndef EFP12_H
#define EFP12_H

#include <ELiPS/EFp6.h>

/**
 * @brief Initializes a EFp12 type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
void EFp12_init(EFp12 *P);

/**
 * @brief Print a EFp12 type struct
 *
 * @param[in]P --a pointer to be printed.
 * @param[in]str --a pointer to be printed.
 */
void EFp12_printf(EFp12 *P,char *str);

/**
 * @brief Set a EFp12 type struct to a EFp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void EFp12_set(EFp12 *ANS,EFp12 *A);

/**
 * @brief Set an unsigned int to a EFp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
void EFp12_set_ui(EFp12 *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a EFp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void EFp12_set_mpn(EFp12 *ANS,mp_limb_t *A);

/**
 * @brief Negate EFp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
void EFp12_set_neg(EFp12 *ANS,EFp12 *A);

/**
 * @brief Compare EFp12 type construct and EFp12 type construct
 *
 * @param[in]A --a pointer in EFp12.
 * @param[in]B --a pointer in EFp12.
 * 
 * @return int --(A=B 0 or other 1)
 */
int  EFp12_cmp(EFp12 *A,EFp12 *B);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in EFp12.
 */
void EFp12_rational_point(EFp12 *P);

/**
 * @brief Generate rational point on G1 (BN12).
 *
 * @param[out]P --a pointer in EFp12.
 */
void BN12_EFp12_generate_G1(EFp12 *P);

/**
 * @brief Generate rational point on G1 (BLS12).
 *
 * @param[out]P --a pointer in EFp12.
 */
void BLS12_EFp12_generate_G1(EFp12 *P);

/**
 * @brief Generate rational point on G2.
 *
 * @param[out]P --a pointer in EFp12.
 */
void EFp12_generate_G2(EFp12 *Q);

/**
 * @brief Doubling a EFp12 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 */
void EFp12_ECD(EFp12 *ANS,EFp12 *P);

/**
 * @brief Doubling a EFp12 type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 */
void EFp12_ECD_lazy(EFp12 *ANS,EFp12 *P);

/**
 * @brief Addition a EFp12 type struct and a EFp12 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFp12.
 * @param[in]P2 --a pointer in EFp12.
 */
void EFp12_ECA(EFp12 *ANS,EFp12 *P1,EFp12 *P2);

/**
 * @brief Addition a EFp12 type struct and a EFp12 type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFp12.
 * @param[in]P2 --a pointer in EFp12.
 */
void EFp12_ECA_lazy(EFp12 *ANS,EFp12 *P1,EFp12 *P2);

/**
 * @brief Scalar multiplication a EFp12 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void EFp12_SCM(EFp12 *ANS,EFp12 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp12 type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
void EFp12_SCM_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar);

#endif

