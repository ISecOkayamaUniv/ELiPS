#ifndef FP12_H
#define FP12_H

#include <ELiPS/fp6.h>


/**
 * @brief Initializes a fp12 type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void fp12_init(fp12_t *A);

/**
 * @brief Print a fp12 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp12_printf(char *str,fp12_t *A);

/**
 * @brief Print a fp12 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp12_println(char *str,fp12_t *A);
extern void fp12_printf_montgomery(char *str,fp12_t *A);
/**
 * @brief Set a fp12 type struct to a fp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp12_set(fp12_t *ANS,fp12_t *A);

/**
 * @brief Set an unsigned int to a fp12 type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp12_set_ui(fp12_t *ANS,unsigned long int UI);

/**
 * @brief Set an unsigned int to a fp12 type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp12_set_ui_ui(fp12_t *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a fp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp12_set_mpn(fp12_t *ANS,mp_limb_t *A);

/**
 * @brief Negate fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void fp12_set_neg(fp12_t *ANS,fp12_t *A);

extern void fp12_to_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp12_mod_montgomery(fp12_t *ANS,fp12_t *A);
/**
 * @brief Set a random number to a fp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void fp12_set_random(fp12_t *ANS,gmp_randstate_t state);

/**
 * @brief Multiplication a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_mul(fp12_t *ANS,fp12_t *A,fp12_t *B);

/**
 * @brief Multiplication a fp12 type struct and a fp12 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_mul_lazy(fp12_t *ANS,fp12_t *A,fp12_t *B);
extern void fp12_mul_lazy_montgomery(fp12_t *ANS,fp12_t *A,fp12_t *B);

/**
 * @brief Multiplication a fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --an unsigned long int.
 */
extern void fp12_mul_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI);

/**
 * @brief Multiplication a fp12 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in mpn.
 */
extern void fp12_mul_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B);

/**
 * @brief Squaring a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_sqr(fp12_t *ANS,fp12_t *A);

/**
 * @brief Squaring a fp12 type struct and a fp12 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_sqr_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_lazy_montgomery(fp12_t *ANS,fp12_t *A);

/**
 * @brief Squaring a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_sqr_karat(fp12_t *ANS,fp12_t *A);

/**
 * @brief Cyclotomic Squaring a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_sqr_cyclotomic(fp12_t *ANS,fp12_t *A);

/**
 * @brief Cyclotomic Squaring a fp12 type struct and a fp12 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_sqr_cyclotomic_lazy(fp12_t *ANS,fp12_t *A);

extern void fp12_sqr_Karabina(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_compressed(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_compressed_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_compressed_lazy_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_recover_g1(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_recover_g1_noninv(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_recover_g1_montrick(fp12_t *A,fp12_t *B,fp12_t *C,fp12_t *D);
extern void fp12_sqr_recover_g1_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_recover_g1_lazy_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_recover_g0(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_recover_g0_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_recover_g0_lazy_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_GS(fp12_t *ANS,fp12_t *A);
extern void fp4_sqr(fp2_t *t0,fp2_t *t1,fp2_t *g0,fp2_t *g1);
extern void fp12_sqr_GS_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_sqr_GS_lazy_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp4_sqr_lazy(fp2_t *t0,fp2_t *t1,fp2_t *g0,fp2_t *g1);
extern void fp4_sqr_lazy_montgomery(fp2_t *t0,fp2_t *t1,fp2_t *g0,fp2_t *g1);
/**
 * @brief Addition a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_add(fp12_t *ANS,fp12_t *A,fp12_t *B);
/**
 * @brief Addition a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_add_nonmod_single(fp12_t *ANS,fp12_t *A,fp12_t *B);

/**
 * @brief Addition a fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --an unsigned long int.
 */
extern void fp12_add_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI);

/**
 * @brief Addition a fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --an unsigned long int.
 */
extern void fp12_add_ui_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI);

/**
 * @brief Addition a fp12 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in mpn.
 */
extern void fp12_add_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B);

/**
 * @brief Subtraction a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_sub(fp12_t *ANS,fp12_t *A,fp12_t *B);

/**
 * @brief Subtraction a fp12 type struct and a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 */
extern void fp12_sub_nonmod_single(fp12_t *ANS,fp12_t *A,fp12_t *B);

/**
 * @brief Subtraction a fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --an unsigned long int.
 */
extern void fp12_sub_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --an unsigned long int.
 */
extern void fp12_sub_ui_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp12 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in mpn.
 */
extern void fp12_sub_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B);

/**
 * @brief Invert a fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp12_inv(fp12_t *ANS,fp12_t *A);

/**
 * @brief Invert a fp12 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp12_inv_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_inv_lazy_montgomery(fp12_t *ANS,fp12_t *A);

/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in fp12.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
extern int  fp12_legendre(fp12_t *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in fp12.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
extern int  fp12_isCNR(fp12_t *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in fp12.
 * @param[out]ANS --a pointer of answer.
 */
extern void fp12_sqrt(fp12_t *ANS,fp12_t *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in fp12.
 * @param[in]A --a pointer in fp12.
 * @param[out]ANS --a pointer in fp12.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern void fp12_pow(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Compare fp12 type construct and fp12 type construct
 *
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in fp12.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp12_cmp(fp12_t *A,fp12_t *B);

/**
 * @brief Compare fp12 type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp12.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  fp12_cmp_ui(fp12_t *A,unsigned long int UI);

/**
 * @brief Compare fp12 type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp12.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp12_cmp_mpn(fp12_t *A,mp_limb_t *B);

/**
 * @brief Compare fp12 type struct and zero
 *
 * @param[in]A --a pointer in fp12.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  fp12_cmp_zero(fp12_t *A);

/**
 * @brief Compare fp12 type struct and one
 *
 * @param[in]A --a pointer in fp12.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  fp12_cmp_one(fp12_t *A);

extern void fp12_frobenius_map_p1(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p1_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p2(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p3(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p3_lazy(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p4(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p6(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p8(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p10(fp12_t *ANS,fp12_t *A);

extern void fp12_frobenius_map_p1_lazy_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p2_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p3_lazy_montgomery(fp12_t *ANS,fp12_t *A);
extern void fp12_frobenius_map_p6_montgomery(fp12_t *ANS,fp12_t *A);

extern int fp12_montgomery_trick(fp12_t *A_inv,fp12_t *A,int n);
#endif
