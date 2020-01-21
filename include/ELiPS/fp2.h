#ifndef FP2_H
#define FP2_H

#include <ELiPS/fp.h>

/**
 * @brief Initializes a fp2_t type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void fp2_init(fp2_t *A);

/**
 * @brief Print a fp2_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp2_printf(char *str,fp2_t *A);
extern void fp2_printf_montgomery(char *str,fp2_t *A);
/**
 * @brief Print a fp2_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp2_println(char *str,fp2_t *A);
extern void fpd2_println(char *str,fpd2_t *A);

/**
 * @brief Set a fp2_t type struct to a fp2_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp2_set(fp2_t *ANS,fp2_t *A);
extern void fpd2_set(fpd2_t *ANS,fpd2_t *A);

/**
 * @brief Set an unsigned int to a fp2_t type struct (x0=UI,x1=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp2_set_ui(fp2_t *ANS,unsigned long int UI);

/**
 * @brief Set an unsigned int to a fp2_t type struct (x0=UI,x1=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp2_set_ui_ui(fp2_t *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a fp2_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp2_set_mpn(fp2_t *ANS,mp_limb_t *A);

/**
 * @brief Negate fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void fp2_set_neg(fp2_t *ANS,fp2_t *A);

extern void fp2_to_montgomery(fp2_t *ANS,fp2_t *A);
extern void fp2_mod_montgomery(fp2_t *ANS,fp2_t *A);
extern void fp2_mod_montgomery_double(fp2_t *ANS,fpd2_t *A);
/**
 * @brief Left Shift fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be shiftted.
 * @param[in]UI --an unsigned long int to shift.
 */
extern void fp2_lshift(fp2_t *ANS,fp2_t *A,unsigned long int UI);

/**
 * @brief Left Shift 2bit fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be shiftted.
 * @param[in]UI --an unsigned long int to shift.
 */
extern void fp2_lshift2(fp2_t *ANS,fp2_t *A);
extern void fp2_div2(fp2_t *ANS,fp2_t *A);
extern void fp2_dbl(fp2_t *ANS, fp2_t *A);

/**
 * @brief Set a random number to a fp2_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void fp2_set_random(fp2_t *ANS,gmp_randstate_t state);

/**
 * @brief Multiplication a fp2_t type struct and a fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_mul(fp2_t *ANS,fp2_t *A,fp2_t *B);

/**
 * @brief Multiplication a fp2_t type struct and a fp2_t type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_mul_lazy(fp2_t *ANS,fp2_t *A,fp2_t *B);
extern void fp2_mul_lazy_montgomery(fp2_t *ANS,fp2_t *A,fp2_t *B);
extern void fp2_mul_nonmod_montgomery(fpd2_t *ANS,fp2_t *A,fp2_t *B);
/**
 * @brief Multiplication a fp2_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp2_mul_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI);

/**
 * @brief Multiplication a fp2_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp2_mul_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B);
extern void fp2_mul_mpn_montgomery(fp2_t *ANS,fp2_t *A,mp_limb_t *B);
/**
 * @brief Multiplication a fp2_t type struct and alpha on prime field A*beta^2
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 */
extern void fp2_mul_basis(fp2_t *ANS,fp2_t *A);

/**
 * @brief Multiplication a fp2_t type struct and alpha on prime field A+B*beta^2
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_add_basis(fp2_t *ANS,fp2_t *A,fp2_t *B);

/**
 * @brief Multiplication a fp2_t type struct and alpha on prime field A-B*beta^2
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_sub_basis(fp2_t *ANS,fp2_t *A,fp2_t *B);

/**
 * @brief Multiplication a fp2_t type struct and alpha on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 */
extern void fp2_mul_basis_lazy(fp2_t *ANS,fp2_t *A);
extern void fp2_mul_basis_lazy_double(fpd2_t *ANS,fpd2_t *A);

extern void fp2_inv_basis(fp2_t *ANS,fp2_t *A);

/**
 * @brief Squaring a fp2_t type struct and a fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_sqr(fp2_t *ANS,fp2_t *A);

/**
 * @brief Squaring a fp2_t type struct and a fp2_t type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_sqr_lazy(fp2_t *ANS,fp2_t *A);
extern void fp2_sqr_lazy_montgomery(fp2_t *ANS,fp2_t *A);

/**
 * @brief Addition a fp2_t type struct and a fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_add(fp2_t *ANS,fp2_t *A,fp2_t *B);

/**
 * @brief Addition a fp2_t type struct and a fp2_t type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_add_nonmod_single(fp2_t *ANS,fp2_t *A,fp2_t *B);
extern void fp2_add_nonmod_double(fpd2_t *ANS,fpd2_t *A,fpd2_t *B);

/**
 * @brief Addition a fp2_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp2_add_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI);

/**
 * @brief Addition a fp2_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp2_add_ui_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI);

/**
 * @brief Addition a fp2_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp2_add_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B);

/**
 * @brief Subtraction a fp2_t type struct and a fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_sub(fp2_t *ANS,fp2_t *A,fp2_t *B);


/**
 * @brief Subtraction a fp2_t type struct and a fp2_t type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 */
extern void fp2_sub_nonmod_single(fp2_t *ANS,fp2_t *A,fp2_t *B);
extern void fp2_sub_nonmod_double(fpd2_t *ANS,fpd2_t *A,fpd2_t *B);

/**
 * @brief Subtraction a fp2_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp2_sub_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp2_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp2_sub_ui_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp2_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp2_sub_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B);

/**
 * @brief Invert a fp2_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp2_inv(fp2_t *ANS,fp2_t *A);

/**
 * @brief Invert a fp2_t type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp2_inv_lazy(fp2_t *ANS,fp2_t *A);

extern void fp2_inv_lazy_montgomery(fp2_t *ANS,fp2_t *A);
/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in fp2_t.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
extern int  fp2_legendre(fp2_t *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in fp2_t.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
extern int  fp2_isCNR(fp2_t *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in fp2_t.
 * @param[out]ANS --a pointer of answer.
 */
extern void fp2_sqrt(fp2_t *ANS,fp2_t *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in fp2_t.
 * @param[in]A --a pointer in fp2_t.
 * @param[out]ANS --a pointer in fp2_t.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern void fp2_pow(fp2_t *ANS,fp2_t *A,mpz_t scalar);

/**
 * @brief Compare fp2_t type construct and fp2_t type construct
 *
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in fp2_t.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp2_cmp(fp2_t *A,fp2_t *B);

/**
 * @brief Compare fp2_t type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp2_t.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  fp2_cmp_ui(fp2_t *A,unsigned long int UI);

/**
 * @brief Compare fp2_t type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp2_t.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp2_cmp_mpn(fp2_t *A,mp_limb_t *B);

/**
 * @brief Compare fp2_t type struct and zero
 *
 * @param[in]A --a pointer in fp2_t.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  fp2_cmp_zero(fp2_t *A);

/**
 * @brief Compare fp2_t type struct and one
 *
 * @param[in]A --a pointer in fp2_t.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  fp2_cmp_one(fp2_t *A);

extern int fp2_montgomery_trick(fp2_t *A_inv,fp2_t *A,int n);
extern int fp2_montgomery_trick_montgomery(fp2_t *A_inv,fp2_t *A,int n);
#endif
