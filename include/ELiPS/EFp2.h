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
 * @brief Initializes a EFpJ2 type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFpJ2_init(EFpJ2 *P);

/**
 * @brief Initializes a EFp2 type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void EFpJT2_init(EFpJT2 *P);

/**
 * @brief Print a EFp2 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFp2_printf(char *str,EFp2 *P);

/**
 * @brief Print a EFp2 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFp2_println(char *str,EFp2 *P);

/**
 * @brief Print a EFpJ2 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void EFpJ2_printf(char *str,EFpJ2 *P);

/**
 * @brief Set a EFp2 type struct to a EFp2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_set(EFp2 *P,EFp2 *A);

/**
 * @brief Set a EFpJ2 type struct to a EFpJ2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFpJ2_set(EFpJ2 *P,EFpJ2 *A);

/**
 * @brief Set a EFpJT2 type struct to a EFpJT2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFpJT2_set(EFpJT2 *P,EFpJT2 *A);

/**
 * @brief Set a EFpJT2 type struct to a EFpJ2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void EFpJT2_to_EFpJ2(EFpJ2 *ANS,EFpJT2 *A);

/**
 * @brief Set a EFpJT2 type struct to a EFpJ2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void EFpJ2_to_EFpJT2(EFpJT2 *ANS,EFpJ2 *A);

/**
 * @brief Set Affine to Jacobian
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_to_EFpJ2(EFpJ2 *ANS,EFp2 *A);

/**
 * @brief Set Jacobian to Affine
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void EFp2_Jacobian(EFp2 *ANS,EFpJ2 *A);

extern void EFp2_mix(EFpJ2 *ANS,EFpJ2 *A,Fp2 *Zi);
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
 * @brief Negate EFp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void EFpJ2_set_neg(EFpJ2 *ANS,EFpJ2 *A);

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
 * @brief Doubling a EFpJ2 type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpJ2.
 */
extern void EFp2_ECD_Jacobian(EFpJ2 *ANS,EFpJ2 *P);

/**
 * @brief Doubling a EFp2 type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp2.
 */
extern void EFp2_ECD_lazy(EFp2 *ANS,EFp2 *P);

/**
 * @brief Doubling a EFpJ2 type struct(Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpJ2.
 */
extern void EFp2_ECD_Jacobian_lazy(EFpJ2 *ANS,EFpJ2 *P);

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
 * @param[in]P1 --a pointer in EFpJ2.
 * @param[in]P2 --a pointer in EFpJ2.
 */
extern void EFp2_ECA_Jacobian(EFpJ2 *ANS,EFpJ2 *P1,EFpJ2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpJ2.
 * @param[in]P2 --a pointer in EFpJ2.
 *
 * @note Z==1
 */
extern void EFp2_ECA_Jacobian_scm(EFpJ2 *ANS,EFpJ2 *P1,EFpJ2 *P2);

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
 * @param[in]P1 --a pointer in EFpJ2.
 * @param[in]P2 --a pointer in EFpJ2.
 */
extern void EFp2_ECA_Jacobian_lazy(EFpJ2 *ANS,EFpJ2 *P1,EFpJ2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpJ2.
 * @param[in]P2 --a pointer in EFpJ2.
 */
extern void EFp2_ECA_Mixture_lazy(EFpJ2 *ANS,EFpJ2 *P1,EFpJ2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpJ2.
 * @param[in]P2 --a pointer in EFpJ2.
 */
extern void EFp2_ECA_Jacobian_table(EFpJ2 *ANS,EFpJ2 *P1,EFpJT2 *P2);

/**
 * @brief Addition a EFp2 type struct and a EFp2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in EFpJ2.
 * @param[in]P2 --a pointer in EFpJ2.
 *
 * @note Z==1
 */
extern void EFp2_ECA_Jacobian_scm_lazy(EFpJ2 *ANS,EFpJ2 *P1,EFpJ2 *P2);

/**
 * @brief Scalar multiplication a EFp2 type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp2.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a EFpJ2 type struct (Jacobian)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpJ2.
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
 * @brief Scalar multiplication a EFpJ2 type struct (Jacobian + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFpJ2.
 * @param[in]scalar --a pointer in mpz.
 */
extern void EFp2_SCM_Jacobian_lazy(EFp2 *ANS,EFp2 *P,mpz_t scalar);
//skew_frobenius_map
extern void EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A);
extern void EFpJ2_skew_frobenius_map_p1(EFpJ2 *ANS,EFpJ2 *A);
extern void EFpJ2_skew_frobenius_map_p2(EFpJ2 *ANS,EFpJ2 *A);
extern void EFpJ2_skew_frobenius_map_p3(EFpJ2 *ANS,EFpJ2 *A);
extern void EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A);

#endif
