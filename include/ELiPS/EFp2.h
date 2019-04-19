#ifndef EFP2_H
#define EFP2_H

#include <ELiPS/EFp.h>

/**
 * @brief Initializes a EFp2 type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFp2_init(EFp2 *P);

/**
 * @brief Initializes a EFpZ2 type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFpZ2_init(EFpZ2 *P);

/**
 * @brief Initializes a EFp2 type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFpZT2_init(EFpZT2 *P);

/**
 * @brief Print a EFp2 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFp2_printf(char *str,EFp2 *P);

/**
 * @brief Print a EFpZ2 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFpZ2_printf(char *str,EFpZ2 *P);

/**
 * @brief Set a EFp2 type struct to a EFp2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_set(EFp2 *P,EFp2 *A);

/**
 * @brief Set a EFpZ2 type struct to a EFpZ2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFpZ2_set(EFpZ2 *P,EFpZ2 *A);

/**
 * @brief Set a EFpZT2 type struct to a EFpZT2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFpZT2_set(EFpZT2 *P,EFpZT2 *A);

/**
 * @brief Set a EFpZT2 type struct to a EFpZ2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void EFpZT2_to_EFpZ2(EFpZ2 *ANS,EFpZT2 *A);

/**
 * @brief Set a EFpZT2 type struct to a EFpZ2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void EFpZ2_to_EFpZT2(EFpZT2 *ANS,EFpZ2 *A);

/**
 * @brief Set Affine to Jacobian
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_to_EFpZ2(EFpZ2 *ANS,EFp2 *A);

/**
 * @brief Set Jacobian to Affine
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_Jacobian(EFp2 *ANS,EFpZ2 *A);

/**
 * @brief Set an unsigned int to a EFp2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void EFp2_set_ui(EFp2 *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a EFp2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_set_mpn(EFp2 *ANS,mp_limb_t *A);

/**
 * @brief Negate EFp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void EFp2_set_neg(EFp2 *ANS,EFp2 *A);

/**
 * @brief Compare EFp2 type construct and EFp2 type construct
 *
 * @param[in]A --a pointer in EFp2.
 * @param[in]B --a pointer in EFp2.
 * 
 * @return int --(A=B 0 or other 1)
 */
int  EFp2_cmp(EFp2 *A,EFp2 *B);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in EFp2.
 */
extern void EFp2_rational_point(EFp2 *P);

/**
 * @brief Doubling a EFp2 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp2.
 */
extern void EFp2_ECD(EFp2 *ANS,EFp2 *P);

/**
 * @brief Doubling a EFpZ2 type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ2.
 */
extern void EFp2_ECD_Jacobian(EFpZ2 *ANS,EFpZ2 *P);

/**
 * @brief Doubling a EFp2 type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp2.
 */
extern void EFp2_ECD_lazy(EFp2 *ANS,EFp2 *P);

/**
 * @brief Doubling a EFpZ2 type struct(Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ2.
 */
extern void EFp2_ECD_Jacobian_lazy(EFpZ2 *ANS,EFpZ2 *P);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFp2.
 * @param[in]P2 --a pointer in EFp2.
 */
extern void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ2.
 * @param[in]P2 --a pointer in EFpZ2.
 */
extern void EFp2_ECA_Jacobian(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ2.
 * @param[in]P2 --a pointer in EFpZ2.
 *
 * @note Z==1
 */
extern void EFp2_ECA_Jacobian_scm(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFp2.
 * @param[in]P2 --a pointer in EFp2.
 */
extern void EFp2_ECA_lazy(EFp2 *ANS,EFp2 *P1,EFp2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ2.
 * @param[in]P2 --a pointer in EFpZ2.
 */
extern void EFp2_ECA_Jacobian_lazy(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ2.
 * @param[in]P2 --a pointer in EFpZ2.
 */
extern void EFp2_ECA_Jacobian_table(EFpZ2 *ANS,EFpZ2 *P1,EFpZT2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpZ2.
 * @param[in]P2 --a pointer in EFpZ2.
 *
 * @note Z==1
 */
extern void EFp2_ECA_Jacobian_scm_lazy(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2);

/**
 * @brief Scalar multiplication a EFp2 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp2.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFpZ2 type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ2.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp2_SCM_Jacobian(EFp2 *ANS,EFp2 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFp2 type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp2.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp2_SCM_lazy(EFp2 *ANS,EFp2 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFpZ2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpZ2.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp2_SCM_Jacobian_lazy(EFp2 *ANS,EFp2 *P,mpz_t scalar);
//skew_frobenius_map
extern void EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A);

#endif
