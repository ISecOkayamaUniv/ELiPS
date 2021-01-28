#ifndef EFP12_H
#define EFP12_H

#include <ELiPS/efp6.h>

/**
 * @brief Initializes a efp12_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp12_init(efp12_t *P);

/**
 * @brief Initializes a efp12_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void sym_init(sym_t *P);

/**
 * @brief Print a efp12_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp12_printf(char *str,efp12_t *P);

/**
 * @brief Print a efp12_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp12_println(char *str,efp12_t *P);

/**
 * @brief Set a efp12_t type struct to a efp12_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp12_set(efp12_t *ANS,efp12_t *A);

/**
 * @brief Set an unsigned int to a efp12_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]UI1 --an unsigned long int to set.
 * @param[in]UI2 --an unsigned long int to set.
 */
extern void efp12_set_ui(efp12_t *ANS,unsigned long int UI1,unsigned long int UI2);

/**
 * @brief Set a mpn type struct to a efp12_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp12_set_mpn(efp12_t *ANS,mp_limb_t *A);

/**
 * @brief Negate efp12_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp12_set_neg(efp12_t *ANS,efp12_t *A);

/**
 * @brief Compare efp12_t type construct and efp12_t type construct
 *
 * @param[in]A --a pointer in efp12_t.
 * @param[in]B --a pointer in efp12_t.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  efp12_cmp(efp12_t *A,efp12_t *B);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in efp12_t.
 */
extern void efp12_rational_point(efp12_t *P);

/**
 * @brief Generate rational point on G1 (BLS12).
 *
 * @param[out]P --a pointer in efp12_t.
 */
extern void bls12_generate_g1(efp12_t *P);

/**
 * @brief Generate rational point on G2.
 *
 * @param[out]P --a pointer in efp12_t.
 */
extern void bls12_generate_g2(efp12_t *Q);

/**
 * @brief Generate symmetric point..
 *
 * @param[in]A --a pointer to be setted.
 */
extern void bls12_generate_symmetric_point(sym_t *A,mpz_t a);

/**
 * @brief Doubling a efp12_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 */
extern void efp12_ecd(efp12_t *ANS,efp12_t *P);

/**
 * @brief Doubling a efp12_t type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 */
extern void efp12_ecd_lazy(efp12_t *ANS,efp12_t *P);

/**
 * @brief Addition a efp12_t type struct and a efp12_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp12_t.
 * @param[in]P2 --a pointer in efp12_t.
 */
extern void efp12_eca(efp12_t *ANS,efp12_t *P1,efp12_t *P2);

/**
 * @brief Addition a efp12_t type struct and a efp12_t type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp12_t.
 * @param[in]P2 --a pointer in efp12_t.
 */
extern void efp12_eca_lazy(efp12_t *ANS,efp12_t *P1,efp12_t *P2);

/**
 * @brief Scalar multiplication a efp12_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void efp12_scm(efp12_t *ANS,efp12_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp12_t type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void efp12_scm_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar);

#endif

