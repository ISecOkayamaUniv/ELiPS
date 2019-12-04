#ifndef BN12_G3_EXP_H
#define BN12_G3_EXP_H

#include <ELiPS/fp12.h>
#include <ELiPS/jsf_naf.h>
#include <ELiPS/time.h>


/**
 * @brief Exponentiation a efp12_t type struct on G3
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g3_exp_plain(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12_t type struct on G3 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g3_exp_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12_t type struct on G3 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g3_exp_2split(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12_t type struct on G3 (GLV-2split + JSF)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g3_exp_2split_JSF(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12_t type struct on G3 for BN12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g3_exp_4split(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12_t type struct on G3 for BN12 (GLV-4split + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bn12_g3_exp_4split_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar);

#endif
