#ifndef LINE_CALCULATE_H
#define LINE_CALCULATE_H

#include <ELiPS/efp2.h>

//Pseudo 8-sparse
extern void Pseudo_8_sparse_mapping(efp_t *P,efp2_t *Q,fp_t *L);
extern void Pseudo_8_sparse_mapping_montgomery(efp_t *P,efp2_t *Q,fp_t *L);
extern void Pseudo_8_sparse_mul(fp12_t *ANS,fp12_t *A,fp12_t *B);
extern void Pseudo_8_sparse_mul_lazy(fp12_t *ANS,fp12_t *A,fp12_t *B);
extern void Pseudo_8_sparse_mul_lazy_montgomery(fp12_t *ANS,fp12_t *A,fp12_t *B);
extern void fp12_6_sparse_mul_lazy(fp12_t *C, fp12_t *A, fp12_t *B);
extern void fp12_6_sparse_mul_lazy_montgomery(fp12_t *C, fp12_t *A, fp12_t *B);
extern void ff_ltt(fp12_t *f,efp2_t *T,efp_t *P,fp_t *L);
extern void ff_ltt_lazy(fp12_t *f,efp2_t *T,efp_t *P,fp_t *L);
extern void ff_ltt_projective_lazy(fp12_t *f,efp2_projective_t *T,efp_t *P);
extern void ff_ltt_lazy_montgomery(fp12_t *f,efp2_t *T,efp_t *P,fp_t *L);
extern void ff_ltt_projective_lazy_montgomery(fp12_t *f,efp2_projective_t *T,efp_t *P);
extern void f_ltq(fp12_t *f,efp2_t *T,efp2_t *Q,efp_t *P,fp_t *L);
extern void f_ltq_lazy(fp12_t *f,efp2_t *T,efp2_t *Q,efp_t *P,fp_t *L);
extern void f_ltq_projective_lazy(fp12_t *f,efp2_projective_t *T,efp2_projective_t *Q,efp_t *P);
extern void f_ltq_lazy_montgomery(fp12_t *f,efp2_t *T,efp2_t *Q,efp_t *P,fp_t *L);
extern void f_ltq_projective_lazy_montgomery(fp12_t *f,efp2_projective_t *T,efp2_projective_t *Q,efp_t *P);

#endif
