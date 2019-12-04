#ifndef BLS12_G3_EXP_H
#define BLS12_G3_EXP_H

#include <ELiPS/fp12.h>
#include <ELiPS/jsf_naf.h>
#include <ELiPS/time.h>


/**
 * @brief Exponentiation a efp12 type struct on G3
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_basic(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_2split(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 (GLV-2split + jsf)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_2split_jsf(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 for bls12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_4split(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 for bls12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_4split_GS(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 for bls12 (GLV-4split + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_4split_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar);

/**
 * @brief Exponentiation a efp12 type struct on G3 for bls12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_g3_exp_4split_GS_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar);


extern void bls12_g3_exp_4split_5naf_interleaving_GS_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar);

extern void bls12_g3_exp_4split_5naf_interleaving_GS_lazy_montgomery(fp12_t *ANS,fp12_t *A,mpz_t scalar);
#endif
