#ifndef FP_H
#define FP_H


#include <ELiPS/bn12_generate_curve.h>
#include <ELiPS/bls12_generate_curve.h>

/**
 * @brief Initializes a Fp type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
void Fp_init(Fp *A);

/**
 * @brief Print a Fp type struct
 *
 * @param[in]A --a pointer to be printed.
 * @param[in]str --a pointer to be printed.
 */
void Fp_printf(Fp *A,char *str);

/**
 * @brief Set a Fp type struct to a Fp type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void Fp_set(Fp *ANS,Fp *A);

/**
 * @brief Set an unsigned int to a Fp type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
void Fp_set_ui(Fp *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a Fp type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
void Fp_set_mpn(Fp *ANS,mp_limb_t *A);

/**
 * @brief Negate Fp type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
void Fp_set_neg(Fp *ANS,Fp *A);

/**
 * @brief Left Shift Fp type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be shiftted.
 * @param[in]UI --an unsigned long int to shift.
 */
void Fp_lshift(Fp *ANS,Fp *A,unsigned long int UI);

/**
 * @brief Set a random number to a Fp type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
void Fp_set_random(Fp *ANS,gmp_randstate_t state);

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
void Lazy_add(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size);

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
void Lazy_sub(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size);

/**
 * @brief Substruction a mpn type struct and a mpn type struct and Modulo ANS(A_size == B_size == FPLIMB2)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in mpn.
 * @param[in]B --a pointer in mpn.
 */
void Lazy_sub_mod(Fp *ANS,mp_limb_t *A,mp_limb_t *B);

/**
 * @brief Multiplication a mpn type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer in mpn.
 * @param[in]A --a pointer in mpn.
 * @param[in]B --a pointer in mpn.
 */
void Lazy_mul(mp_limb_t *ANS,mp_limb_t *A,mp_limb_t *B);

/**
 * @brief Modulo a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer in mpn.
 * @param[in]A --a pointer in mpn.
 * @param[in]A_size --a size of A.
 */
void Lazy_mod(mp_limb_t *ans,mp_limb_t *a,mp_size_t size_a);

/**
 * @brief Multiplication a Fp type struct and a Fp type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in Fp.
 */
void Fp_mul(Fp *ANS,Fp *A,Fp *B);

/**
 * @brief Multiplication a Fp type struct and a Fp type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in Fp.
 */
void Fp_mul_final(Fp *ANS,Fp *A,Fp *B);

/**
 * @brief Multiplication a Fp type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --an unsigned long int.
 */
void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI);

/**
 * @brief Multiplication a Fp type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in mpn.
 */
void Fp_mul_mpn(Fp *ANS,Fp *A,mp_limb_t *B);

/**
 * @brief Addition a Fp type struct and a Fp type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in Fp.
 */
void Fp_add(Fp *ANS,Fp *A,Fp *B);

/**
 * @brief Addition a Fp type struct and a Fp type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in Fp.
 */
void Fp_add_final(Fp *ANS,Fp *A,Fp *B);

/**
 * @brief Addition a Fp type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --an unsigned long int.
 */
void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI);

/**
 * @brief Addition a Fp type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in mpn.
 */
void Fp_add_mpn(Fp *ANS,Fp *A,mp_limb_t *B);

/**
 * @brief Subtraction a Fp type struct and a Fp type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in Fp.
 */
void Fp_sub(Fp *ANS,Fp *A,Fp *B);

/**
 * @brief Subtraction a Fp type struct and a Fp type struct on prime field (Always mod)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in Fp.
 */
void Fp_sub_final(Fp *ANS,Fp *A,Fp *B);

/**
 * @brief Subtraction a Fp type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --an unsigned long int.
 */
void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI);

/**
 * @brief Subtraction a Fp type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in mpn.
 */
void Fp_sub_mpn(Fp *ANS,Fp *A,mp_limb_t *B);

/**
 * @brief Invert a Fp type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
void Fp_inv(Fp *ANS,Fp *A);

/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in Fp.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
extern int  Fp_legendre(Fp *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in Fp.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
extern int  Fp_isCNR(Fp *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in Fp.
 * @param[out]ANS --a pointer of answer.
 */
void Fp_sqrt(Fp *ANS,Fp *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in Fp.
 * @param[in]A --a pointer in Fp.
 * @param[out]ANS --a pointer in Fp.
 * 
 * @return int --(A=UI 0 or other 1)
 */
void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar);

/**
 * @brief Power a by mpn type struct
 *
 * @param[in]a --a pointer in Fp.
 * @param[in]r --a pointer in mpn.
 * @param[in]n --a size of r.
 * @param[out]ans --a pointer answer.
 * 
 * @return int --(A=B 0 or other 1)
 */
void Fp_pow_mpn(Fp *ans,Fp *a,mp_limb_t *r,mp_size_t n);

/**
 * @brief Compare Fp type construct and Fp type construct
 *
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in Fp.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp_cmp(Fp *A,Fp *B);

/**
 * @brief Compare Fp type construct and unsigned long int
 *
 * @param[in]A --a pointer in Fp.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  Fp_cmp_ui(Fp *A,unsigned long int UI);

/**
 * @brief Compare Fp type construct and mpn type construct
 *
 * @param[in]A --a pointer in Fp.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  Fp_cmp_mpn(Fp *A,mp_limb_t *B);

/**
 * @brief Compare Fp type struct and zero
 *
 * @param[in]A --a pointer in Fp.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  Fp_cmp_zero(Fp *A);

/**
 * @brief Compare Fp type struct and one
 *
 * @param[in]A --a pointer in Fp.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  Fp_cmp_one(Fp *A);

#endif
