#ifndef BLS12_G3_EXP_H
#define BLS12_G3_EXP_H

#include <ELiPS/fp12.h>
#include <ELiPS/jsf_naf.h>
#include <ELiPS/time.h>


/**
 * @brief Exponentiation a efp12 type struct on G3 (basic type)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_basic(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 (4split_5naf_interleaving_GS)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp(fp12_t *ANS,fp12_t *A,mpz_t scalar);
#endif
