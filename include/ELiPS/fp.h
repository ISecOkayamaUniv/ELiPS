#ifndef FP_H
#define FP_H

#include <ELiPS/mpn.h>

/**
 * @brief Initializes a fp_t type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void fp_init(fp_t *A);

/**
 * @brief Print a fp_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp_printf(char *str,fp_t *A);
extern void fpd_printf(char *str,fpd_t *A);

/**
 * @brief Print a fp_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp_println(char *str,fp_t *A);

extern void fp_printf_montgomery(char *str,fp_t *A);
/**
 * @brief Set a fp_t type struct to a fp_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp_set(fp_t *ANS,fp_t *A);
extern void fpd_set(fpd_t *ANS,fpd_t *A);

/**
 * @brief Set an unsigned int to a fp_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp_set_ui(fp_t *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a fp_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp_set_mpn(fp_t *ANS,mp_limb_t *A);

/**
 * @brief Negate fp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void fp_set_neg(fp_t *ANS,fp_t *A);

/**
 * @brief Left Shift fp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be shiftted.
 * @param[in]UI --an unsigned long int to shift.
 */
extern void fp_lshift(fp_t *ANS,fp_t *A,unsigned long int UI);

/**
 * @brief Left Shift 1bit fp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be shiftted.
 * @param[in]UI --an unsigned long int to shift.
 */
extern void fp_l1shift(fp_t *ANS,fp_t *A);

extern void fp_r1shift(fp_t *ANS,fp_t *A);

/**
 * @brief Set a random number to a fp_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void fp_set_random(fp_t *ANS,gmp_randstate_t state);

extern void fp_set_random_montgomery(fp_t *ANS,gmp_randstate_t state);


/**
 * @brief Modulo a fp_t type struct on prime field (Montgomery Reduction)
 *
 * @param[out]ANS --a pointer of answer in fp_t.
 * @param[in]A --a pointer in fp_t.
 */
extern void pre_montgomery();
extern void fp_mulmod_montgomery(fp_t *ANS,fp_t *A,fp_t *B);
extern void fp_sqrmod_montgomery(fp_t *ANS,fp_t *A);
extern void fp_mod_montgomery(fp_t *ANS,fp_t *A);
extern void fp_to_montgomery(fp_t *ANS, fp_t *A);


/**
 * @brief Modulo a mpn type struct on prime field
 *
 * @param[out]ans --a pointer of answer in fp_t.
 * @param[in]a --a pointer to be mod.
 * @param[in]a_size --a size of a.
 */
extern void fp_mod(fp_t *ans,mp_limb_t *a,mp_size_t size_a);

/**
 * @brief Modulo a mpn type struct
 *
 * @param[out]ans --a pointer of answer in fp_t.
 * @param[in]a --a pointer to be mod.
 * @param[in]a_size --a size of a.
 * @param[in]ui --a unsigned long int to mod.
 */
extern void fp_mod_ui(fp_t *ans,mp_limb_t *a,mp_size_t size_a,unsigned long int UI);

/**
 * @brief Multiplication a fp_t type struct and a fp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in fp_t.
 */
extern void fp_mul(fp_t *ANS,fp_t *A,fp_t *B);

extern void fp_mul_nonmod(fpd_t *ANS,fp_t *A,fp_t *B);

/**
 * @brief Multiplication a fp_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp_mul_ui(fp_t *ANS,fp_t *A,unsigned long int UI);

/**
 * @brief Multiplication a fp_t type struct and an unsigned long int
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp_mul_ui_nonmod_single(fp_t *ANS,fp_t *A,unsigned long int UI);

/**
 * @brief Multiplication a fp_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp_mul_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B);

/**
 * @brief Squraring a mpn type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer in fp_t.
 * @param[in]A --a pointer in fp_t.
 */
extern void fp_sqr(fp_t *ANS,fp_t *A);

extern void fp_sqr_nonmod(fpd_t *ANS,fp_t *A);

/**
 * @brief Addition a fp_t type struct and a fp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in fp_t.
 */
extern void fp_add(fp_t *ANS,fp_t *A,fp_t *B);

/**
 * @brief Addition a mpn type struct and a mpn type struct
 *
 * @param[out]ANS --a pointer of answer in mpn.
 * @param[out]ANS_size --a size of answer.
 * @param[in]A --a pointer in mpn.
 * @param[in]A_size --a size of A.
 * @param[in]B --a pointer in mpn.
 * @param[in]B_size --a size of B.
 */
extern  void fp_add_nonmod_single(fp_t *ANS,fp_t *A,fp_t *B);
extern void fp_add_nonmod_double(fpd_t *ANS,fpd_t *A,fpd_t *B);

/**
 * @brief Addition a fp_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp_add_ui(fp_t *ANS,fp_t *A,unsigned long int UI);

/**
 * @brief Addition a fp_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp_add_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B);

/**
 * @brief Subtraction a fp_t type struct and a fp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in fp_t.
 */
extern void fp_sub(fp_t *ANS,fp_t *A,fp_t *B);
/**
 * @brief Substruction a mpn type struct and a mpn type struct (A_size >= B_size)
 *
 * @param[out]ANS --a pointer of answer in mpn.
 * @param[out]ANS_size --a size of answer.
 * @param[in]A --a pointer in mpn.
 * @param[in]A_size --a size of A.
 * @param[in]B --a pointer in mpn.
 * @param[in]B_size --a size of B.
 */
extern void fp_sub_nonmod_single(fp_t *ANS,fp_t *A,fp_t *B);
extern void fp_sub_nonmod_double(fpd_t *ANS,fpd_t *A,fpd_t *B);


/**
 * @brief Subtraction a fp_t type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --an unsigned long int.
 */
extern void fp_sub_ui(fp_t *ANS,fp_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp_t type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in mpn.
 */
extern void fp_sub_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B);

/**
 * @brief Invert a fp_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp_inv(fp_t *ANS,fp_t *A);

extern void fp_inv_montgomery(fp_t *ANS,fp_t *A);


/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in fp_t.
 *
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
extern int  fp_legendre(fp_t *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in fp_t.
 *
 * @return int --a CNR (0 or 1 or -1)
 */
extern int  fp_isCNR(fp_t *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in fp_t.
 * @param[out]ANS --a pointer of answer.
 */
extern void fp_sqrt(fp_t *ANS,fp_t *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in fp_t.
 * @param[in]A --a pointer in fp_t.
 * @param[out]ANS --a pointer in fp_t.
 *
 * @return int --(A=UI 0 or other 1)
 */
extern void fp_pow(fp_t *ANS,fp_t *A,mpz_t scalar);

/**
 * @brief Power a by mpn type struct
 *
 * @param[in]a --a pointer in fp_t.
 * @param[in]r --a pointer in mpn.
 * @param[in]n --a size of r.
 * @param[out]ans --a pointer answer.
 *
 * @return int --(A=B 0 or other 1)
 */
extern void fp_pow_mpn(fp_t *ans,fp_t *a,mp_limb_t *r,mp_size_t n);

/**
 * @brief Compare fp_t type construct and fp_t type construct
 *
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in fp_t.
 *
 * @return int --(A=B 0 or other 1)
 */
extern int  fp_cmp(fp_t *A,fp_t *B);

/**
 * @brief Compare fp_t type construct and unsigned long int
 *
 * @param[in]A --a pointer in fp_t.
 * @param[in]UI --an unsigned long int.
 *
 * @return int --(A=UI 0 or other 1)
 */
extern int  fp_cmp_ui(fp_t *A,unsigned long int UI);

/**
 * @brief Compare fp_t type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp_t.
 * @param[in]B --a pointer in mpn.
 *
 * @return int --(A=B 0 or other 1)
 */
extern int  fp_cmp_mpn(fp_t *A,mp_limb_t *B);

/**
 * @brief Compare fp_t type struct and zero
 *
 * @param[in]A --a pointer in fp_t.
 *
 * @return int --(one 0 or other 1)
 */
extern int  fp_cmp_zero(fp_t *A);

/**
 * @brief Compare fp_t type struct and one
 *
 * @param[in]A --a pointer in fp_t.
 *
 * @return int --(zero 0 or other 1)
 */
extern int  fp_cmp_one(fp_t *A);

extern int fp_montgomery_trick(fp_t *A_inv,fp_t *A,int n);
extern int fp_montgomery_trick_montgomery(fp_t *A_inv,fp_t *A,int n);
extern void fp_lshift_ui_nonmod_single(fp_t *ANS,fp_t *A,int s);
extern void fp_lshift_ui_nonmod_double(fpd_t *ANS,fpd_t *A,int s);
extern void fpd_mod_montgomery(fp_t *ANS,fpd_t *A);
extern int fp_legendre_sqrt(fp_t *ANS,fp_t *A);
#endif
