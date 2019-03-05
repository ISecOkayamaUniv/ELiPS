#ifndef FP2_H
#define FP2_H

#include <ELiPS/Fp.h>

/**
 * @brief Initializes a Fp2 type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void Fp2_init(Fp2 *A);

/**
 * @brief Print a Fp2 type struct
 *
 * @param[in]A --a pointer to be printed.
 * @param[in]str --a pointer to be printed.
 */
extern void Fp2_printf(Fp2 *A,char *str);

/**
 * @brief Set a Fp2 type struct to a Fp2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void Fp2_set(Fp2 *ANS,Fp2 *A);

/**
 * @brief Set an unsigned int to a Fp2 type struct (x0=UI,x1=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void Fp2_set_ui(Fp2 *ANS,unsigned long int UI);

/**
 * @brief Set an unsigned int to a Fp2 type struct (x0=UI,x1=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void Fp2_set_ui_ui(Fp2 *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a Fp2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void Fp2_set_mpn(Fp2 *ANS,mp_limb_t *A);

/**
 * @brief Negate Fp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void Fp2_set_neg(Fp2 *ANS,Fp2 *A);

/**
 * @brief Left Shift Fp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be shiftted.
 * @param[in]UI --an unsigned long int to shift.
 */
extern void Fp2_lshift(Fp2 *ANS,Fp2 *A,unsigned long int UI);

/**
 * @brief Set a random number to a Fp2 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state);

/**
 * @brief Multiplication a Fp2 type struct and a Fp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Multiplication a Fp2 type struct and a Fp2 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_mul_lazy(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Multiplication a Fp2 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --an unsigned long int.
 */
extern void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);

/**
 * @brief Multiplication a Fp2 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp2_mul_mpn(Fp2 *ANS,Fp2 *A,mp_limb_t *B);

/**
 * @brief Multiplication a Fp2 type struct and alpha on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 */
extern void Fp2_mul_basis(Fp2 *ANS,Fp2 *A);

/**
 * @brief Multiplication a Fp2 type struct and alpha on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 */
extern void Fp2_mul_basis_lazy(Fp2 *ANS,Fp2 *A);

extern void Fp2_inv_basis(Fp2 *ANS,Fp2 *A);

/**
 * @brief Squaring a Fp2 type struct and a Fp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_sqr(Fp2 *ANS,Fp2 *A);

/**
 * @brief Squaring a Fp2 type struct and a Fp2 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_sqr_lazy(Fp2 *ANS,Fp2 *A);

/**
 * @brief Addition a Fp2 type struct and a Fp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Addition a Fp2 type struct and a Fp2 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_add_lazy(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Addition a Fp2 type struct and a Fp2 type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_add_final(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Addition a Fp2 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --an unsigned long int.
 */
extern void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);

/**
 * @brief Addition a Fp2 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --an unsigned long int.
 */
extern void Fp2_add_ui_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);

/**
 * @brief Addition a Fp2 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp2_add_mpn(Fp2 *ANS,Fp2 *A,mp_limb_t *B);

/**
 * @brief Subtraction a Fp2 type struct and a Fp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Subtraction a Fp2 type struct and a Fp2 type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_sub_final(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Subtraction a Fp2 type struct and a Fp2 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 */
extern void Fp2_sub_lazy(Fp2 *ANS,Fp2 *A,Fp2 *B);

/**
 * @brief Subtraction a Fp2 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --an unsigned long int.
 */
extern void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);

/**
 * @brief Subtraction a Fp2 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --an unsigned long int.
 */
extern void Fp2_sub_ui_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);

/**
 * @brief Subtraction a Fp2 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp2_sub_mpn(Fp2 *ANS,Fp2 *A,mp_limb_t *B);

/**
 * @brief Invert a Fp2 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void Fp2_inv(Fp2 *ANS,Fp2 *A);

/**
 * @brief Invert a Fp2 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void Fp2_inv_lazy(Fp2 *ANS,Fp2 *A);

/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in Fp2.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
extern int  Fp2_legendre(Fp2 *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in Fp2.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
extern int  Fp2_isCNR(Fp2 *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in Fp2.
 * @param[out]ANS --a pointer of answer.
 */
extern void Fp2_sqrt(Fp2 *ANS,Fp2 *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in Fp2.
 * @param[in]A --a pointer in Fp2.
 * @param[out]ANS --a pointer in Fp2.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar);

/**
 * @brief Compare Fp2 type construct and Fp2 type construct
 *
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in Fp2.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp2_cmp(Fp2 *A,Fp2 *B);

/**
 * @brief Compare Fp2 type construct and mpn type construct
 *
 * @param[in]A --a pointer in Fp2.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI);

/**
 * @brief Compare Fp2 type construct and mpn type construct
 *
 * @param[in]A --a pointer in Fp2.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp2_cmp_mpn(Fp2 *A,mp_limb_t *B);

/**
 * @brief Compare Fp2 type struct and zero
 *
 * @param[in]A --a pointer in Fp2.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  Fp2_cmp_zero(Fp2 *A);

/**
 * @brief Compare Fp2 type struct and one
 *
 * @param[in]A --a pointer in Fp2.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  Fp2_cmp_one(Fp2 *A);

#endif
