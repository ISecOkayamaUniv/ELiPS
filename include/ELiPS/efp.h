#ifndef EFP_H
#define EFP_H

#include <ELiPS/fp12.h>
/**
 * @brief Initializes a efp_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp_init(efp_t *P);

/**
 * @brief Initializes a efp_jacobian_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp_jacobian_init(efp_jacobian_t *P);

/**
 * @brief Initializes a efp_jacobian_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp_projective_init(efp_projective_t *P);

/**
 * @brief Print a efp_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp_printf(char *str,efp_t *P);

/**
 * @brief Print a efp_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp_println(char *str,efp_t *P);

/**
 * @brief Print a efp_jacobian_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp_projective_printf(char *str,efp_projective_t *P);

/**
 * @brief Print a efp_jacobian_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp_jacobian_printf(char *str,efp_jacobian_t *P);

/**
 * @brief Set a efp_t type struct to a efp_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_set(efp_t *P,efp_t *A);

/**
 * @brief Set a efp_jacobian_t type struct to a efp_jacobian_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_projective_set(efp_projective_t *P,efp_projective_t *A);

/**
 * @brief Set a efp_jacobian_t type struct to a efp_jacobian_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_jacobian_set(efp_jacobian_t *P,efp_jacobian_t *A);
/**
 * @brief Set Affine to jacobian
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_affine_to_projective(efp_projective_t *ANS,efp_t *A);
/**
 * @brief Set Affine to jacobian
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_affine_to_jacobian(efp_jacobian_t *ANS,efp_t *A);

extern void efp_affine_to_jacobian_montgomery(efp_jacobian_t *ANS,efp_t *A);
/**
 * @brief Set jacobian to Affine
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_projective_to_affine(efp_t *ANS,efp_projective_t *A);

/**
 * @brief Set jacobian to Affine
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_jacobian_to_affine(efp_t *ANS,efp_jacobian_t *A);

extern void efp_jacobian_to_affine_montgomery(efp_t *ANS,efp_jacobian_t *A);

extern void efp_jacobian_to_mixture_noninv(efp_jacobian_t *ANS,efp_jacobian_t *A,fp_t *Zi);
extern void efp_jacobian_to_mixture_noninv_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *A,fp_t *Zi);

extern void efp_mod_montgomery(efp_t *ANS,efp_t *A);

extern void efp_to_montgomery(efp_t *ANS,efp_t *A);

/**
 * @brief Set an unsigned int to a efp_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void efp_set_ui(efp_t *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a efp_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp_set_mpn(efp_t *ANS,mp_limb_t *A);

/**
 * @brief Negate efp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp_set_neg(efp_t *ANS,efp_t *A);

/**
 * @brief Negate efp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp_projective_set_neg(efp_projective_t *ANS,efp_projective_t *A);

/**
 * @brief Negate efp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp_jacobian_set_neg(efp_jacobian_t *ANS,efp_jacobian_t *A);

/**
 * @brief Compare efp_t type construct and efp_t type construct
 *
 * @param[in]A --a pointer in efp_t.
 * @param[in]B --a pointer in efp_t.
 *
 * @return int --(A=B 0 or other 1)
 */
extern int  efp_cmp(efp_t *A,efp_t *B);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in efp_t.
 */
extern void efp_set_random(efp_t *P);

/**
 * @brief Doubling a efp_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp_t.
 */
extern void efp_ecd(efp_t *ANS,efp_t *P);

extern void efp_ecd_jacobian_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P);
/**
 * @brief Addition a efp_t type struct and a efp_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp_t.
 * @param[in]P2 --a pointer in efp_t.
 */
extern void efp_eca(efp_t *ANS,efp_t *P1,efp_t *P2);


extern void efp_eca_jacobian_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2);


extern void efp_eca_mixture_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2);

/**
 * @brief Scalar multiplication a efp_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void efp_scm(efp_t *ANS,efp_t *P,mpz_t scalar);

//skew_frobenius_map
extern void efp_skew_frobenius_map_p2(efp_t *ANS,efp_t *A);
extern void efp_jacobian_skew_frobenius_map_p2(efp_jacobian_t *ANS,efp_jacobian_t *A);
extern void efp_jacobian_skew_frobenius_map_p2_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *A);
extern void efp_ecd_jacobian_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P);
extern void efp_eca_jacobian_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2);
extern void efp_eca_mixture_lazy_montgomery_ignore_inf(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2);
#endif
