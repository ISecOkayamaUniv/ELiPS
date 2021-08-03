#ifndef BLS12_SYM_H
#define BLS12_SYM_H

#include <ELiPS/bls12_g1_scm.h>
#include <ELiPS/bls12_g2_scm.h>
#include <ELiPS/time.h>

/**
 * @brief Scalar multiplication a efp12_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp12_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void bls12_sym_scm(sym_t *ANS, sym_t *A, mpz_t scalar);

#endif
