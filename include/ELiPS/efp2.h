#ifndef EFP2_H
#define EFP2_H

#include <ELiPS/efp.h>

/**
 * @brief Initializes a efp2_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp2_init(efp2_t *P);

/**
 * @brief Initializes a efp2_jacobian_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp2_jacobian_init(efp2_jacobian_t *P);

extern void efp2_projective_init(efp2_projective_t *P);
/**
 * @brief Print a efp2_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp2_printf(char *str,efp2_t *P);
extern void efp2_projective_printf(char *str,efp2_projective_t *P);
extern void efp2_printf_montgomery(char *str,efp2_t *P);
extern void efp2_projective_printf_affine(char *str,efp2_projective_t *P);
extern void efp2_projective_printf_affine_montgomery(char *str,efp2_projective_t *P);

/**
 * @brief Print a efp2_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp2_println(char *str,efp2_t *P);

/**
 * @brief Print a efp2_jacobian_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp2_jacobian_printf(char *str,efp2_jacobian_t *P);
extern void efp2_jacobian_printf_montgomery(char *str,efp2_jacobian_t *P);
extern void efp2_projective_printf_montgomery(char *str,efp2_projective_t *P);
/**
 * @brief Set a efp2_t type struct to a efp2_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_set(efp2_t *P,efp2_t *A);

/**
 * @brief Set a efp2_jacobian_t type struct to a efp2_jacobian_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_jacobian_set(efp2_jacobian_t *ANS,efp2_jacobian_t *A);

/**
 * @brief Set a efp2_jacobian_t type struct to a efp2_projective_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_projective_set(efp2_projective_t *ANS,efp2_projective_t *A);

/**
 * @brief Set Affine to jacobian
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_affine_to_jacobian(efp2_jacobian_t *ANS,efp2_t *A);
/**
 * @brief Set Affine to jacobian
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_affine_to_projective(efp2_projective_t *ANS,efp2_t *A);
extern void efp2_affine_to_projective_montgomery(efp2_projective_t *ANS,efp2_t *A);

extern void efp2_affine_to_jacobian_montgomery(efp2_jacobian_t *ANS,efp2_t *A);

/**
 * @brief Set jacobian to Affine
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_jacobian_to_affine(efp2_t *ANS,efp2_jacobian_t *A);
/**
 * @brief Set jacobian to Affine
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_projective_to_affine(efp2_t *ANS,efp2_projective_t *A);

extern void efp2_mix(efp2_jacobian_t *ANS,efp2_jacobian_t *A,fp2_t *Zi);

extern void efp2_jacobian_to_affine_montgomery(efp2_t *ANS,efp2_jacobian_t *A);
extern void efp2_projective_to_affine_montgomery(efp2_t *ANS,efp2_projective_t *A);

extern void efp2_mix_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A,fp2_t *Zi);
extern void efp2_to_montgomery(efp2_t *ANS,efp2_t *A);
extern void efp2_projective_to_montgomery(efp2_projective_t *ANS,efp2_projective_t *A);
extern void efp2_mod_montgomery(efp2_t *ANS,efp2_t *A);
extern void efp2_projective_mod_montgomery(efp2_projective_t *ANS,efp2_projective_t *A);
/**
 * @brief Set an unsigned int to a efp2_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void efp2_set_ui(efp2_t *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a efp2_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp2_set_mpn(efp2_t *ANS,mp_limb_t *A);

/**
 * @brief Negate efp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp2_set_neg(efp2_t *ANS,efp2_t *A);

/**
 * @brief Negate efp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp2_jacobian_set_neg(efp2_jacobian_t *ANS,efp2_jacobian_t *A);

/**
 * @brief Compare efp2_t type construct and efp2_t type construct
 *
 * @param[in]A --a pointer in efp2_t.
 * @param[in]B --a pointer in efp2_t.
 *
 * @return int --(A=B 0 or other 1)
 */
extern int  efp2_cmp(efp2_t *A,efp2_t *B);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in efp2_t.
 */
extern void efp2_rational_point(efp2_t *P);

/**
 * @brief Doubling a efp2_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp2_t.
 */
extern void efp2_ecd(efp2_t *ANS,efp2_t *P);
extern void efp2_ecd_lazy_montgomery(efp2_t *ANS,efp2_t *P);

/**
 * @brief Doubling a efp2_jacobian_t type struct(jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp2_jacobian_t.
 */
extern void efp2_ecd_jacobian_lazy_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *P);

/**
 * @brief Addition a efp2_t type struct and a efp2_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp2_t.
 * @param[in]P2 --a pointer in efp2_t.
 */
extern void efp2_eca(efp2_t *ANS,efp2_t *P1,efp2_t *P2);

extern void efp2_eca_lazy_montgomery(efp2_t *ANS,efp2_t *P1,efp2_t *P2);


/**
 * @brief Addition a efp2_t type struct and a efp2_t type struct (jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp2_jacobian_t.
 * @param[in]P2 --a pointer in efp2_jacobian_t.
 */
extern void efp2_eca_jacobian_lazy_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *P1,efp2_jacobian_t *P2);

/**
 * @brief Addition a efp2_t type struct and a efp2_t type struct (jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp2_jacobian_t.
 * @param[in]P2 --a pointer in efp2_jacobian_t.
 */
extern void efp2_eca_mixture_lazy_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *P1,efp2_jacobian_t *P2);

/**
 * @brief Scalar multiplication a efp2_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp2_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void efp2_scm(efp2_t *ANS,efp2_t *P,mpz_t scalar);


//skew_frobenius_map
extern void efp2_skew_frobenius_map_p1(efp2_t *ANS,efp2_t *A);
extern void efp2_skew_frobenius_map_p2(efp2_t *ANS,efp2_t *A);
extern void efp2_skew_frobenius_map_p3(efp2_t *ANS,efp2_t *A);
extern void efp2_jacobian_skew_frobenius_map_p1(efp2_jacobian_t *ANS,efp2_jacobian_t *A);
extern void efp2_jacobian_skew_frobenius_map_p2(efp2_jacobian_t *ANS,efp2_jacobian_t *A);
extern void efp2_jacobian_skew_frobenius_map_p3(efp2_jacobian_t *ANS,efp2_jacobian_t *A);
extern void efp2_skew_frobenius_map_p1_montgomery(efp2_t *ANS,efp2_t *A);
extern void efp2_jacobian_skew_frobenius_map_p1_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A);
extern void efp2_jacobian_skew_frobenius_map_p2_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A);
extern void efp2_jacobian_skew_frobenius_map_p3_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A);
extern void efp2_skew_frobenius_map_p10(efp2_t *ANS,efp2_t *A);

#endif
