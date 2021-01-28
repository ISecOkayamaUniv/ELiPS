#include <ELiPS/bls12_sym.h>

void bls12_sym_scm(sym_t *ANS,sym_t *A,mpz_t scalar){
    bls12_g1_scm(&ANS->p,&A->p,scalar);
    bls12_g2_scm(&ANS->q,&A->q,scalar);
}