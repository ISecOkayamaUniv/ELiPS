#ifndef FP6_H
#define FP6_H

#include <ELiPS/Fp2.h>

/**
 * @brief Initializes a Fp6 type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void Fp6_init(Fp6 *A);

/**
 * @brief Print a Fp6 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void Fp6_printf(char *str,Fp6 *A);

/**
 * @brief Print a Fp6 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void Fp6_println(char *str,Fp6 *A);

extern void Fp6_printf_montgomery(char *str,Fp6 *A);
/**
 * @brief Set a Fp6 type struct to a Fp6 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void Fp6_set(Fp6 *ANS,Fp6 *A);

/**
 * @brief Set an unsigned int to a Fp6 type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void Fp6_set_ui(Fp6 *ANS,unsigned long int UI);


/**
 * @brief Set an unsigned int to a Fp6 type struct (x0=UI,x1=UI,x2=UI)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void Fp6_set_ui_ui(Fp6 *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a Fp6 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void Fp6_set_mpn(Fp6 *ANS,mp_limb_t *A);

/**
 * @brief Negate Fp6 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void Fp6_set_neg(Fp6 *ANS,Fp6 *A);

extern void Fp6_to_montgomery(Fp6 *ANS,Fp6 *A);
extern void Fp6_mod_montgomery(Fp6 *ANS,Fp6 *A);
/**
 * @brief Set a random number to a Fp6 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state);

/**
 * @brief Multiplication a Fp6 type struct and a Fp6 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_mul(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Multiplication a Fp6 type struct and a Fp6 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_mul_lazy(Fp6 *ANS,Fp6 *A,Fp6 *B);
extern void Fp6_mul_lazy_montgomery(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Multiplication a Fp6 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --an unsigned long int.
 */
extern void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);

/**
 * @brief Multiplication a Fp6 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp6_mul_mpn(Fp6 *ANS,Fp6 *A,mp_limb_t *B);

/**
 * @brief Multiplication a Fp6 type struct and beta on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 */
extern void Fp6_mul_basis(Fp6 *ANS,Fp6 *A);

/**
 * @brief Multiplication a Fp6 type struct and beta on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 */
extern void Fp6_mul_basis_lazy(Fp6 *ANS,Fp6 *A);

/**
 * @brief Squaring a Fp6 type struct and a Fp6 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_sqr(Fp6 *ANS,Fp6 *A);

/**
 * @brief Squaring a Fp6 type struct and a Fp6 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_sqr_lazy(Fp6 *ANS,Fp6 *A);
extern void Fp6_sqr_lazy_montgomery(Fp6 *ANS,Fp6 *A);

/**
 * @brief Addition a Fp6 type struct and a Fp6 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Addition a Fp6 type struct and a Fp6 type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_add_final(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Addition a Fp6 type struct and a Fp6 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_add_lazy(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Addition a Fp6 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --an unsigned long int.
 */
extern void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);

/**
 * @brief Addition a Fp6 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --an unsigned long int.
 */
extern void Fp6_add_ui_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);

/**
 * @brief Addition a Fp6 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp6_add_mpn(Fp6 *ANS,Fp6 *A,mp_limb_t *B);

/**
 * @brief Subtraction a Fp6 type struct and a Fp6 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Subtraction a Fp6 type struct and a Fp6 type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_sub_final(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Subtraction a Fp6 type struct and a Fp6 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 */
extern void Fp6_sub_lazy(Fp6 *ANS,Fp6 *A,Fp6 *B);

/**
 * @brief Subtraction a Fp6 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --an unsigned long int.
 */
extern void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);

/**
 * @brief Subtraction a Fp6 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --an unsigned long int.
 */
extern void Fp6_sub_ui_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);

/**
 * @brief Subtraction a Fp6 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp6_sub_mpn(Fp6 *ANS,Fp6 *A,mp_limb_t *B);

/**
 * @brief Invert a Fp6 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void Fp6_inv(Fp6 *ANS,Fp6 *A);

/**
 * @brief Invert a Fp6 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void Fp6_inv_lazy(Fp6 *ANS,Fp6 *A);
extern void Fp6_inv_lazy_montgomery(Fp6 *ANS,Fp6 *A);

/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in Fp6.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
int  Fp6_legendre(Fp6 *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in Fp6.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
int  Fp6_isCNR(Fp6 *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in Fp6.
 * @param[out]ANS --a pointer of answer.
 */
extern void Fp6_sqrt(Fp6 *ANS,Fp6 *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in Fp6.
 * @param[in]A --a pointer in Fp6.
 * @param[out]ANS --a pointer in Fp6.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern void Fp6_pow(Fp6 *ANS,Fp6 *A,mpz_t scalar);

/**
 * @brief Compare Fp6 type construct and Fp6 type construct
 *
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in Fp6.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp6_cmp(Fp6 *A,Fp6 *B);

/**
 * @brief Compare Fp6 type construct and mpn type construct
 *
 * @param[in]A --a pointer in Fp6.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI);

/**
 * @brief Compare Fp6 type construct and mpn type construct
 *
 * @param[in]A --a pointer in Fp6.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp6_cmp_mpn(Fp6 *A,mp_limb_t *B);

/**
 * @brief Compare Fp6 type struct and zero
 *
 * @param[in]A --a pointer in Fp6.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  Fp6_cmp_zero(Fp6 *A);

/**
 * @brief Compare Fp6 type struct and one
 *
 * @param[in]A --a pointer in Fp6.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  Fp6_cmp_one(Fp6 *A);

#endif
