#ifndef SPARSE_H
#define SPARSE_H

#include <ELiPS/EFp2.h>

//Pseudo 8-sparse
void Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L);
void Pseudo_8_sparse_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);
void Pseudo_8_sparse_mul_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B);
void ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L);
void ff_ltt_lazy(Fp12 *f,EFp2 *T,EFp *P,Fp *L);
void f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);
void f_ltq_lazy(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);

#endif
