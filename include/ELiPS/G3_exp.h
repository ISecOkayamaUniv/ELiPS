#ifndef G3_EXP_H
#define G3_EXP_H

#include <ELiPS/Fp12.h>
#include <ELiPS/JSF.h>
#include <ELiPS/time.h>


/**
 * @brief Exponentiation a EFp12 type struct on G3
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void Fp12_G3_EXP_plain(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void Fp12_G3_EXP_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 (GLV-2split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void Fp12_G3_EXP_2split(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 (GLV-2split + JSF)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void Fp12_G3_EXP_2split_JSF(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 for BN12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BN12_Fp12_G3_EXP_4split(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 for BN12 (GLV-4split + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BN12_Fp12_G3_EXP_4split_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 for BLS12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_Fp12_G3_EXP_4split(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 for BLS12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_Fp12_G3_EXP_4split_GS(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 for BLS12 (GLV-4split + Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_Fp12_G3_EXP_4split_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar);

/**
 * @brief Exponentiation a EFp12 type struct on G3 for BLS12 (GLV-4split)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in EFp12.
 * @param[in]scalar --a pointer in mpz.
 */
extern void BLS12_Fp12_G3_EXP_4split_GS_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar);


extern void BLS12_Fp12_G3_EXP_4split_5NAF_interleaving_GS_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar);
#endif
