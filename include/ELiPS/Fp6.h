#ifndef FP6_H
#define FP6_H

#include <ELiPS/fp2.h>

/**
 * @brief Initializes a fp6_t type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void fp6_init(fp6_t *A);

/**
 * @brief Print a fp6_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp6_printf(char *str,fp6_t *A);

/**
 * @brief Print a fp6_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp6_println(char *str,fp6_t *A);

extern void fp6_printf_montgomery(char *str,fp6_t *A);
/**
 * @brief Set a fp6_t type struct to a fp6_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp6_set(fp6_t *ANS,fp6_t *A);

/**
 * @brief Set an unsigned int to a fp6_t type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp6_set_ui(fp6_t *ANS,unsigned long int UI);


/**
 * @brief Set an unsigned int to a fp6_t type struct (x0=UI,x1=UI,x2=UI)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp6_set_ui_ui(fp6_t *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a fp6_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp6_set_mpn(fp6_t *ANS,mp_limb_t *A);

/**
 * @brief Negate fp6_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void fp6_set_neg(fp6_t *ANS,fp6_t *A);

extern void fp6_to_montgomery(fp6_t *ANS,fp6_t *A);
extern void fp6_mod_montgomery(fp6_t *ANS,fp6_t *A);
/**
 * @brief Set a random number to a fp6_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void fp6_set_random(fp6_t *ANS,gmp_randstate_t state);

/**
 * @brief Multiplication a fp6_t type struct and a fp6_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_mul(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Multiplication a fp6_t type struct and a fp6_t type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_mul_lazy(fp6_t *ANS,fp6_t *A,fp6_t *B);
extern void fp6_mul_lazy_montgomery(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Multiplication a fp6_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp6_mul_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI);

/**
 * @brief Multiplication a fp6_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp6_mul_mpn(fp6_t *ANS,fp6_t *A,mp_limb_t *B);

/**
 * @brief Multiplication a fp6_t type struct and beta on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 */
extern void fp6_mul_basis(fp6_t *ANS,fp6_t *A);

/**
 * @brief Multiplication a fp6_t type struct and beta on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 */
extern void fp6_mul_basis_lazy(fp6_t *ANS,fp6_t *A);

/**
 * @brief Squaring a fp6_t type struct and a fp6_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_sqr(fp6_t *ANS,fp6_t *A);

/**
 * @brief Squaring a fp6_t type struct and a fp6_t type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_sqr_lazy(fp6_t *ANS,fp6_t *A);
extern void fp6_sqr_lazy_montgomery(fp6_t *ANS,fp6_t *A);

/**
 * @brief Addition a fp6_t type struct and a fp6_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_add(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Addition a fp6_t type struct and a fp6_t type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_add_final(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Addition a fp6_t type struct and a fp6_t type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_add_lazy(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Addition a fp6_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp6_add_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI);

/**
 * @brief Addition a fp6_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp6_add_ui_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI);

/**
 * @brief Addition a fp6_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp6_add_mpn(fp6_t *ANS,fp6_t *A,mp_limb_t *B);

/**
 * @brief Subtraction a fp6_t type struct and a fp6_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_sub(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Subtraction a fp6_t type struct and a fp6_t type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_sub_final(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Subtraction a fp6_t type struct and a fp6_t type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 */
extern void fp6_sub_lazy(fp6_t *ANS,fp6_t *A,fp6_t *B);

/**
 * @brief Subtraction a fp6_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp6_sub_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp6_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp6_sub_ui_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp6_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp6_sub_mpn(fp6_t *ANS,fp6_t *A,mp_limb_t *B);

/**
 * @brief Invert a fp6_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp6_inv(fp6_t *ANS,fp6_t *A);

/**
 * @brief Invert a fp6_t type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp6_inv_lazy(fp6_t *ANS,fp6_t *A);
extern void fp6_inv_lazy_montgomery(fp6_t *ANS,fp6_t *A);

/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in fp6_t.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
int  fp6_legendre(fp6_t *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in fp6_t.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
int  fp6_isCNR(fp6_t *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in fp6_t.
 * @param[out]ANS --a pointer of answer.
 */
extern void fp6_sqrt(fp6_t *ANS,fp6_t *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in fp6_t.
 * @param[in]A --a pointer in fp6_t.
 * @param[out]ANS --a pointer in fp6_t.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern void fp6_pow(fp6_t *ANS,fp6_t *A,mpz_t scalar);

/**
 * @brief Compare fp6_t type construct and fp6_t type construct
 *
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in fp6_t.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp6_cmp(fp6_t *A,fp6_t *B);

/**
 * @brief Compare fp6_t type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp6_t.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  fp6_cmp_ui(fp6_t *A,unsigned long int UI);

/**
 * @brief Compare fp6_t type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp6_t.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp6_cmp_mpn(fp6_t *A,mp_limb_t *B);

/**
 * @brief Compare fp6_t type struct and zero
 *
 * @param[in]A --a pointer in fp6_t.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  fp6_cmp_zero(fp6_t *A);

/**
 * @brief Compare fp6_t type struct and one
 *
 * @param[in]A --a pointer in fp6_t.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  fp6_cmp_one(fp6_t *A);

#endif
