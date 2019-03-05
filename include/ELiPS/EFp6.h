#ifndef EFP6_H
#define EFP6_H

#include <ELiPS/EFp2.h>

/**
 * @brief Initializes a EFp6 type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFp6_init(EFp6 *P);

/**
 * @brief Print a EFp6 type struct
 *
 * @param[in]P --a pointer to be printed.
 * @param[in]str --a pointer to be printed.
 */
extern void EFp6_printf(EFp6 *P,char *str);

/**
 * @brief Set a EFp6 type struct to a EFp6 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp6_set(EFp6 *ANS,EFp6 *A);

/**
 * @brief Set an unsigned int to a EFp6 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void EFp6_set_ui(EFp6 *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a EFp6 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp6_set_mpn(EFp6 *ANS,mp_limb_t *A);

/**
 * @brief Negate EFp6 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void EFp6_set_neg(EFp6 *ANS,EFp6 *A);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in EFp6.
 */
extern void EFp6_rational_point(EFp6 *P);

/**
 * @brief Doubling a EFp6 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp6.
 */
extern void EFp6_ECD(EFp6 *ANS,EFp6 *P);

/**
 * @brief Addition a EFp6 type struct and a EFp6 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFp6.
 * @param[in]P2 --a pointer in EFp6.
 */
extern void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2);

/**
 * @brief Scalar multiplication a EFp6 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp6.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar);

#endif

