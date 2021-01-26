#ifndef MPN_H
#define MPN_H

#include <ELiPS/define.h>


/**
 * @brief Initializes a mpn type struct
 *
 * @param[in]a --a pointer to be initialized.
 * @param[in]size --a size of A.
 */
extern void mpn_init(mp_limb_t *a,mp_size_t size);

/**
 * @brief Set a charcter to a mpn type struct
 *
 * @param[out]ans --a pointer to be setted.
 * @param[in]mp_size --a size of ans.
 * @param[in]str --a pointer to set.
 */
extern void mpn_set_char(mp_limb_t *ans,mp_size_t mp_size,char *str);

/**
 * @brief Set an unsigned long int to a mpn type struct
 *
 * @param[out]ans --a pointer to be setted.
 * @param[in]mp_size --a size of ans.
 * @param[in]ui --a unsigned long int to set.
 */
extern void mpn_set_ui(mp_limb_t *ans,mp_size_t size,unsigned long int ui);

/**
 * @brief Set an unsigned long int to a mpn type struct
 *
 * @param[out]ans --a pointer to be setted.
 * @param[in]mp_size --a size of ans.
 * @param[in]a --a pointer to set.
 */
extern void mpn_set_mpz(mp_limb_t *ans,mpz_t a);

/**
 * @brief Modulo a mpn type struct on prime field
 *
 * @param[out]ans --a pointer of answer in mp_limb_t.
 * @param[in]a --a pointer to be mod.
 * @param[in]size_a --a size of a.
 */
extern void mpn_mod(mp_limb_t *ans,mp_limb_t *a,mp_size_t size_a);

/**
 * @brief Compare mpn type construct and char
 *
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]str --a pointer of char.
 * 
 * @return int --(A==str 0 or other 1)
 */
extern int mpn_cmp_char(mp_limb_t *a,char *str);

/**
 * @brief Compare mpn type construct and unsigned long int
 *
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]ui --an unsigned long int.
 * 
 * @return int --(A==UI 0 or other 1)
 */
extern int mpn_cmp_ui(mp_limb_t *a,mp_size_t size,unsigned long int ui);

/**
 * @brief Make sure sizeA or sizeB.
 *
 * @param[in]a --a pointer in mpn.
 * @param[in]sizeA --a size of patternA.
 * @param[in]sizeB --a size of patternB.
 * 
 * @return int --(sizeB 0 or sizeA 1)
 */
extern int mpn_chk_limb(mp_limb_t *a,mp_size_t sizeA,mp_size_t sizeB);

/**
 * @brief Left Shift mpn type struct on prime field.
 *
 * @param[out]ans --a pointer of answer.
 * @param[in]a --a pointer to be shiftted.
 * @param[in]size --a size of a.
 * @param[in]sizeB --a long int to shift.
 *
 * @note L>64 is OK
 */
extern void mpn_lshift_ext(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,long int L);

/**
 * @brief righ Shift mpn type struct on prime field.
 *
 * @param[out]ans --a pointer of answer.
 * @param[in]a --a pointer to be shiftted.
 * @param[in]size --a size of a.
 * @param[in]sizeB --a long int to shift.
 *
 * @note L>64 is OK
 */
extern void mpn_rshift_ext(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,long int L);

/**
 * @brief double mpn type struct on prime field.
 *
 * @param[out]ans --a pointer of answer.
 * @param[in]a --a pointer to be doubled.
 * @param[in]size --a size of a.
 */
extern void mpn_dbl(mp_limb_t *ans,mp_limb_t *a,mp_size_t size);

/**
 * @brief Addition a mpn type struct and char type
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]str --a char type.
 */
extern void mpn_add_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,char *str);

/**
 * @brief Substruction a mpn type struct and char type
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]str --a char type.
 */
extern void mpn_sub_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,char *str);

/**
 * @brief multiplication a mpn type struct and char type
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]str --a char type.
 */
extern void mpn_mul_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,char *str);

/**
 * @brief Addition a mpn type struct and unsigned long int.
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]ui --an unsigned long int.
 */
extern void mpn_add_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,unsigned long int ui);

/**
 * @brief Substruction a mpn type struct and unsigned long int.
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]ui --an unsigned long int.
 */
extern void mpn_sub_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,unsigned long int ui);

/**
 * @brief Multiplication a mpn type struct and unsigned long int.
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer in mpn.
 * @param[in]size --a size of a.
 * @param[in]ui --an unsigned long int.
 */
extern void mpn_mul_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,unsigned long int ui);

/**
 * @brief Exponentiation a mpn type struct.
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[out]ans_size --a size of aans.
 * @param[in]a --a pointer in mpn.
 * @param[in]a_size --a size of a.
 * @param[in]r --a pointer in mpn.
 * @param[in]n --a size of r.
 */
extern void mpn_pow(mp_limb_t *ans,mp_size_t ans_size,mp_limb_t *a,mp_size_t a_size,mp_limb_t *r,mp_size_t n);

/**
 * @brief Exponentiation a mpn type struct.
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[out]ans_size --a size of aans.
 * @param[in]a --a pointer in mpn.
 * @param[in]a_size --a size of a.
 * @param[in]str --a pointer of char.
 */
extern void mpn_pow_ui(mp_limb_t *ans,mp_size_t ans_size,mp_limb_t *a,mp_size_t a_size,char *str);

/**
 * @brief Division a mpn type struct on prime field
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer to be division.
 * @param[in]a_size --a size of a.
 * @param[in]str --a pointer of char.
 */
extern void mpn_tdiv_q_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size_a,char *str);

/**
 * @brief Division a mpn type struct on prime field
 *
 * @param[out]ans --a pointer of answer in mpn.
 * @param[in]a --a pointer to be division.
 * @param[in]a_size --a size of a.
 * @param[in]UI --an unsigned long int.
 */
extern void mpn_tdiv_q_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size_a,unsigned long int UI);

/**
 * @brief Inversion a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer in mpn.
 * @param[in]A --a pointer to be invert.
 * @param[in]p --a pointer of prime.
 */
extern void mpn_invert(mp_limb_t *ANS,mp_limb_t *A,mp_limb_t *p);

#endif
