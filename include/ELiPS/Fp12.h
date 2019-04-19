#ifndef FP12_H
#define FP12_H

#include <ELiPS/Fp6.h>


/**
 * @brief Initializes a Fp12 type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void Fp12_init(Fp12 *A);

/**
 * @brief Print a Fp12 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void Fp12_printf(char *str,Fp12 *A);

/**
 * @brief Set a Fp12 type struct to a Fp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void Fp12_set(Fp12 *ANS,Fp12 *A);

/**
 * @brief Set an unsigned int to a Fp12 type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void Fp12_set_ui(Fp12 *ANS,unsigned long int UI);

/**
 * @brief Set an unsigned int to a Fp12 type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void Fp12_set_ui_ui(Fp12 *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a Fp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void Fp12_set_mpn(Fp12 *ANS,mp_limb_t *A);

/**
 * @brief Negate Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void Fp12_set_neg(Fp12 *ANS,Fp12 *A);

/**
 * @brief Set a random number to a Fp12 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state);

/**
 * @brief Multiplication a Fp12 type struct and a Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);

/**
 * @brief Multiplication a Fp12 type struct and a Fp12 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_mul_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B);

/**
 * @brief Multiplication a Fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --an unsigned long int.
 */
extern void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);

/**
 * @brief Multiplication a Fp12 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp12_mul_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B);

/**
 * @brief Squaring a Fp12 type struct and a Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_sqr(Fp12 *ANS,Fp12 *A);

/**
 * @brief Squaring a Fp12 type struct and a Fp12 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_sqr_lazy(Fp12 *ANS,Fp12 *A);

/**
 * @brief Squaring a Fp12 type struct and a Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_sqr_karat(Fp12 *ANS,Fp12 *A);

/**
 * @brief Cyclotomic Squaring a Fp12 type struct and a Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_sqr_cyclotomic(Fp12 *ANS,Fp12 *A);

/**
 * @brief Cyclotomic Squaring a Fp12 type struct and a Fp12 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_sqr_cyclotomic_lazy(Fp12 *ANS,Fp12 *A);

extern void Fp12_sqr_Karabina(Fp12 *ANS,Fp12 *A);
extern void Fp12_sqr_compressed(Fp12 *ANS,Fp12 *A);
extern void Fp12_sqr_compressed_lazy(Fp12 *ANS,Fp12 *A);
extern void Fp12_sqr_recover_g1(Fp12 *ANS,Fp12 *A);
extern void Fp12_sqr_recover_g1_lazy(Fp12 *ANS,Fp12 *A);
extern void Fp12_sqr_recover_g0(Fp12 *ANS,Fp12 *A);
extern void Fp12_sqr_recover_g0_lazy(Fp12 *ANS,Fp12 *A);
extern void Fp12_sqr_GS(Fp12 *ANS,Fp12 *A);
extern void Fp4_sqr(Fp2 *t0,Fp2 *t1,Fp2 *g0,Fp2 *g1);
extern void Fp12_sqr_GS_lazy(Fp12 *ANS,Fp12 *A);
extern void Fp4_sqr_lazy(Fp2 *t0,Fp2 *t1,Fp2 *g0,Fp2 *g1);
/**
 * @brief Addition a Fp12 type struct and a Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B);

/**
 * @brief Addition a Fp12 type struct and a Fp12 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_add_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B);

/**
 * @brief Addition a Fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --an unsigned long int.
 */
extern void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);

/**
 * @brief Addition a Fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --an unsigned long int.
 */
extern void Fp12_add_ui_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);

/**
 * @brief Addition a Fp12 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp12_add_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B);

/**
 * @brief Subtraction a Fp12 type struct and a Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B);

/**
 * @brief Subtraction a Fp12 type struct and a Fp12 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 */
extern void Fp12_sub_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B);

/**
 * @brief Subtraction a Fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --an unsigned long int.
 */
extern void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);

/**
 * @brief Subtraction a Fp12 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --an unsigned long int.
 */
extern void Fp12_sub_ui_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);

/**
 * @brief Subtraction a Fp12 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in mpn.
 */
extern void Fp12_sub_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B);

/**
 * @brief Invert a Fp12 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void Fp12_inv(Fp12 *ANS,Fp12 *A);

/**
 * @brief Invert a Fp12 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void Fp12_inv_lazy(Fp12 *ANS,Fp12 *A);

/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in Fp12.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
extern int  Fp12_legendre(Fp12 *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in Fp12.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
extern int  Fp12_isCNR(Fp12 *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in Fp12.
 * @param[out]ANS --a pointer of answer.
 */
extern void Fp12_sqrt(Fp12 *ANS,Fp12 *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in Fp12.
 * @param[in]A --a pointer in Fp12.
 * @param[out]ANS --a pointer in Fp12.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern void Fp12_pow(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Compare Fp12 type construct and Fp12 type construct
 *
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in Fp12.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp12_cmp(Fp12 *A,Fp12 *B);

/**
 * @brief Compare Fp12 type construct and mpn type construct
 *
 * @param[in]A --a pointer in Fp12.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI);

/**
 * @brief Compare Fp12 type construct and mpn type construct
 *
 * @param[in]A --a pointer in Fp12.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp12_cmp_mpn(Fp12 *A,mp_limb_t *B);

/**
 * @brief Compare Fp12 type struct and zero
 *
 * @param[in]A --a pointer in Fp12.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  Fp12_cmp_zero(Fp12 *A);

/**
 * @brief Compare Fp12 type struct and one
 *
 * @param[in]A --a pointer in Fp12.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  Fp12_cmp_one(Fp12 *A);

extern void Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p1_lazy(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p3_lazy(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A);

#endif
