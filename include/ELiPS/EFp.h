#ifndef EFP_H
#define EFP_H

#include <ELiPS/Fp12.h>

/**
 * @brief Initializes a EFp type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFp_init(EFp *P);

/**
 * @brief Initializes a EFpZ type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFpZ_init(EFpZ *P);

/**
 * @brief Initializes a EFpZ type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFpZT_init(EFpZT *P);

/**
 * @brief Print a EFp type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFp_printf(char *str,EFp *P);

/**
 * @brief Print a EFpZ type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFpZ_printf(char *str,EFpZ *P);

/**
 * @brief Print a EFpZ type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFpZT_printf(char *str,EFpZT *P);

/**
 * @brief Set a EFp type struct to a EFp type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp_set(EFp *P,EFp *A);

/**
 * @brief Set a EFpZ type struct to a EFpZ type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFpZ_set(EFpZ *P,EFpZ *A);

/**
 * @brief Set a EFpZT type struct to a EFpZT type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFpZT_set(EFpZT *P,EFpZT *A);

/**
 * @brief Set a EFpZT type struct to a EFpZ type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFpZT_to_EFpZ(EFpZ *P,EFpZT *A);

/**
 * @brief Set a EFpZT type struct to a EFpZ type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void EFpZ_to_EFpZT(EFpZT *ANS,EFpZ *A);
/**
 * @brief Set Affine to Jacobian
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp_to_EFpZ(EFpZ *ANS,EFp *A);

/**
 * @brief Set Jacobian to Affine
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp_Jacobian(EFp *ANS,EFpZ *A);

/**
 * @brief Set an unsigned int to a EFp type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void EFp_set_ui(EFp *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a EFp type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp_set_mpn(EFp *ANS,mp_limb_t *A);

/**
 * @brief Negate EFp type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void EFp_set_neg(EFp *ANS,EFp *A);

/**
 * @brief Compare EFp type construct and EFp type construct
 *
 * @param[in]A --a pointer in EFp.
 * @param[in]B --a pointer in EFp.
 * 
 * @return int --(A=B 0 or other 1)
 */
int  EFp_cmp(EFp *A,EFp *B);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in EFp.
 */
extern void EFp_rational_point(EFp *P);

/**
 * @brief Doubling a EFp type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp.
 */
extern void EFp_ECD(EFp *ANS,EFp *P);

/**
 * @brief Doubling a EFpZ type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ.
 */
extern void EFp_ECD_Jacobian(EFpZ *ANS,EFpZ *P);

/**
 * @brief Doubling a EFp type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp.
 */
extern void EFp_ECD_lazy(EFp *ANS,EFp *P);

/**
 * @brief Doubling a EFpZ type struct(Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ.
 */
extern void EFp_ECD_Jacobian_lazy(EFpZ *ANS,EFpZ *P);

/**
 * @brief Addition a EFp type struct and a EFp type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFp.
 * @param[in]P2 --a pointer in EFp.
 */
extern void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);

/**
 * @brief Addition a EFp type struct and a EFp type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ.
 * @param[in]P2 --a pointer in EFpZ.
 */
extern void EFp_ECA_Jacobian(EFpZ *ANS,EFpZ *P1,EFpZ *P2);

/**
 * @brief Addition a EFp type struct and a EFp type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFp.
 * @param[in]P2 --a pointer in EFp.
 */
extern void EFp_ECA_lazy(EFp *ANS,EFp *P1,EFp *P2);

/**
 * @brief Addition a EFp type struct and a EFp type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ.
 * @param[in]P2 --a pointer in EFpZ.
 */
extern void EFp_ECA_Jacobian_lazy(EFpZ *ANS,EFpZ *P1,EFpZ *P2);

/**
 * @brief Addition a EFp type struct and a EFp type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ.
 * @param[in]P2 --a pointer in EFpZ.
 */
extern void EFp_ECA_Jacobian_table(EFpZ *ANS,EFpZ *P1,EFpZT *P2);

/**
 * @brief Scalar multiplication a EFp type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFpZ type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp_SCM_Jacobian(EFp *ANS,EFp *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp_SCM_lazy(EFp *ANS,EFp *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFpZ type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp_SCM_Jacobian_lazy(EFp *ANS,EFp *P,mpz_t scalar);
//skew_frobenius_map
extern void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A);

#endif
