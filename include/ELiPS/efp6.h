#ifndef EFP6_H
#define EFP6_H

#include <ELiPS/efp2.h>

/**
 * @brief Initializes a efp6_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp6_init(efp6_t *P);

/**
 * @brief Print a efp6_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp6_printf(char *str, efp6_t *P);

/**
 * @brief Print a efp6_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp6_println(char *str, efp6_t *P);

/**
 * @brief Set a efp6_t type struct to a efp6_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp6_set(efp6_t *ANS, efp6_t *A);

/**
 * @brief Set an unsigned int to a efp6_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void efp6_set_ui(efp6_t *ANS, unsigned long int UI);

/**
 * @brief Set a mpn type struct to a efp6_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp6_set_mpn(efp6_t *ANS, mp_limb_t *A);

/**
 * @brief Negate efp6_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp6_set_neg(efp6_t *ANS, efp6_t *A);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in efp6_t.
 */
extern void efp6_rational_point(efp6_t *P);

/**
 * @brief Doubling a efp6_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp6_t.
 */
extern void efp6_ecd(efp6_t *ANS, efp6_t *P);

/**
 * @brief Addition a efp6_t type struct and a efp6_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp6_t.
 * @param[in]P2 --a pointer in efp6_t.
 */
extern void efp6_eca(efp6_t *ANS, efp6_t *P1, efp6_t *P2);

/**
 * @brief Scalar multiplication a efp6_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp6_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void efp6_scm(efp6_t *ANS, efp6_t *P, mpz_t scalar);

#endif
