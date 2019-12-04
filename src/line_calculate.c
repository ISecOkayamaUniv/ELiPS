#include <ELiPS/line_calculate.h>
//Pseudo 8-sparse
void Pseudo_8_sparse_mapping(efp_t *P,efp2_t *Q,fp_t *L){
    efp2_t Tmp_Q;
	efp2_init(&Tmp_Q);
	efp_t Tmp_P;
	efp_init(&Tmp_P);
	fp_t A,B,C,D,c;
	fp_init(&A);
	fp_init(&B);
	fp_init(&C);
	fp_init(&D);
	fp_init(&c);
	
	efp_set(&Tmp_P,P);
	efp2_set(&Tmp_Q,Q);
	
	fp_mul(&A,&Tmp_P.x,&Tmp_P.y);
	fp_inv(&A,&A);
	fp_mul(&B,&Tmp_P.x,&Tmp_P.x);
	fp_mul(&B,&B,&A);
	fp_mul(&C,&Tmp_P.y,&A);
	fp_mul(&D,&B,&B);
	
	fp2_mul_mpn(&Q->x,&Tmp_Q.x,D.x0);
	fp_mul(&c,&B,&D);
	fp2_mul_mpn(&Q->y,&Tmp_Q.y,c.x0);
	
	fp_mul(&P->x,&D,&Tmp_P.x);
	fp_set(&P->y,&P->x);
	
	fp_mul(L,&C,&Tmp_P.y);
	fp_mul(L,L,L);
	fp_mul(L,L,&C);
}
void Pseudo_8_sparse_mapping_montgomery(efp_t *P,efp2_t *Q,fp_t *L){
    efp2_t Tmp_Q;
	efp2_init(&Tmp_Q);
	efp_t Tmp_P;
	efp_init(&Tmp_P);
	fp_t A,B,C,D,c;
	fp_init(&A);
	fp_init(&B);
	fp_init(&C);
	fp_init(&D);
	fp_init(&c);
	
	efp_set(&Tmp_P,P);
	efp2_set(&Tmp_Q,Q);
	
	fp_mulmod_montgomery(&A,&Tmp_P.x,&Tmp_P.y);
	fp_inv_montgomery(&A,&A);
	fp_mulmod_montgomery(&B,&Tmp_P.x,&Tmp_P.x);
	fp_mulmod_montgomery(&B,&B,&A);
	fp_mulmod_montgomery(&C,&Tmp_P.y,&A);
	fp_mulmod_montgomery(&D,&B,&B);
	
	fp2_mul_mpn_montgomery(&Q->x,&Tmp_Q.x,D.x0);
	fp_mulmod_montgomery(&c,&B,&D);
	fp2_mul_mpn_montgomery(&Q->y,&Tmp_Q.y,c.x0);
	
	fp_mulmod_montgomery(&P->x,&D,&Tmp_P.x);
	fp_set(&P->y,&P->x);
	
	fp_mulmod_montgomery(L,&C,&Tmp_P.y);
	fp_mulmod_montgomery(L,L,L);
	fp_mulmod_montgomery(L,L,&C);
}

void Pseudo_8_sparse_mul(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2;
    static fp12_t tmp3_fp12;
	fp2_mul(&tmp1_fp2,&A->x0.x0,&B->x1.x0);		//tmp1_fp2=b3*f0
	fp2_mul(&tmp2_fp2,&A->x0.x1,&B->x1.x1);		//tmp2_fp2=b4*f1
	fp2_add(&tmp3_fp2,&A->x0.x0,&A->x0.x1);		//tmp3_fp2=f0+f1
	fp2_add(&tmp4_fp2,&B->x1.x0,&B->x1.x1);		//tmp4_fp2=b3+b4

	fp2_mul(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);			//tmp3_fp2=tmp3_fp2*tmp4_fp2
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp1_fp2);			//tmp3_fp2=tmp3_fp2-tmp1_fp2
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp2_fp2);			//tmp3_fp2=tmp3_fp2-tmp2_fp2

	fp2_add(&tmp3_fp12.x1.x1,&tmp3_fp2,&A->x1.x1);	//ans[γ^3]=tmp3_fp2+f4
	fp2_add(&tmp1_fp2,&tmp1_fp2,&A->x1.x0);		//tmp1_fp2=tmp1_fp2+f3
	fp2_mul(&tmp3_fp2,&A->x0.x2,&B->x1.x1);		//tmp3_fp2=b4*f2
	fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);

	fp2_add(&tmp3_fp12.x1.x0,&tmp1_fp2,&tmp3_fp2);		//ans[γ]=tmp1_fp2+tmp3_fp2
	fp2_add(&tmp1_fp2,&tmp2_fp2,&A->x1.x2);		//tmp1_fp2=tmp2_fp2+f5
	fp2_mul(&tmp2_fp2,&A->x0.x2,&B->x1.x0);		//tmp2_fp2=b3*f2	
	fp2_add(&tmp3_fp12.x1.x2,&tmp1_fp2,&tmp2_fp2);		//ans[γ^5]=tmp1_fp2+tmp2_fp2
	fp2_mul(&tmp1_fp2,&A->x1.x0,&B->x1.x0);		//tmp1_fp2=b3*f3
	fp2_mul(&tmp2_fp2,&A->x1.x1,&B->x1.x1);		//tmp2_fp2=b4*f4
	fp2_add(&tmp3_fp2,&A->x1.x0,&A->x1.x1);		//tmp3_fp2=f3+f4
	fp2_mul(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);			//tmp3_fp2=tmp3_fp2*tmp4_fp2

	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp1_fp2);			//tmp3_fp2=tmp3_fp2-tmp1_fp2
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp2_fp2);			//tmp3_fp2=tmp3_fp2-tmp2_fp2
	fp2_add(&tmp3_fp12.x0.x2,&tmp3_fp2,&A->x0.x2);	//ans[γ^4]=tmp3_fp2+f4
	fp2_add(&tmp1_fp2,&tmp1_fp2,&A->x0.x1);		//tmp1_fp2=tmp1_fp2+f1
	fp2_mul(&tmp3_fp2,&A->x1.x2,&B->x1.x1);		//tmp3_fp2=b4*f5
	fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);			//tmp3_fp2=tmp3_fp2*(α+1)
	fp2_add(&tmp3_fp12.x0.x1,&tmp1_fp2,&tmp3_fp2);		//ans[γ^2]=tmp1_fp2+tmp3_fp2
	fp2_mul(&tmp1_fp2,&A->x1.x2,&B->x1.x0);		//tmp1_fp2=b3*f5	
	fp2_add(&tmp1_fp2,&tmp1_fp2,&tmp2_fp2);			//tmp1_fp2=tmp1_fp2+tmp2_fp2
	fp2_mul_basis(&tmp1_fp2,&tmp1_fp2);			//tmp1_fp2=tmp1_fp2*(α+1)
	fp2_add(&tmp3_fp12.x0.x0,&tmp1_fp2,&A->x0.x0);	//ans[1]=tmp1_fp2+f0	
	fp12_set(ANS,&tmp3_fp12);	
}
void Pseudo_8_sparse_mul_lazy(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2;
    static fp12_t tmp3_fp12;
	fp2_mul_lazy(&tmp1_fp2,&A->x0.x0,&B->x1.x0);
	fp2_mul_lazy(&tmp2_fp2,&A->x0.x1,&B->x1.x1);
	fp2_add_lazy(&tmp3_fp2,&A->x0.x0,&A->x0.x1);
	fp2_add_lazy(&tmp4_fp2,&B->x1.x0,&B->x1.x1);

	fp2_mul_lazy(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp1_fp2);		
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp2_fp2);
	fp2_add(&tmp3_fp12.x1.x1,&tmp3_fp2,&A->x1.x1);

	fp2_add(&tmp1_fp2,&tmp1_fp2,&A->x1.x0);
	fp2_mul_lazy(&tmp3_fp2,&A->x0.x2,&B->x1.x1);
	fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);
	fp2_add(&tmp3_fp12.x1.x0,&tmp1_fp2,&tmp3_fp2);

	fp2_add(&tmp1_fp2,&tmp2_fp2,&A->x1.x2);
	fp2_mul_lazy(&tmp2_fp2,&A->x0.x2,&B->x1.x0);
	fp2_add(&tmp3_fp12.x1.x2,&tmp1_fp2,&tmp2_fp2);

	fp2_mul_lazy(&tmp1_fp2,&A->x1.x0,&B->x1.x0);//tmp1
	fp2_mul_lazy(&tmp2_fp2,&A->x1.x1,&B->x1.x1);//tmp2

	fp2_add_lazy(&tmp3_fp2,&A->x1.x0,&A->x1.x1);
	fp2_mul_lazy(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp1_fp2);
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp2_fp2);
	fp2_add(&tmp3_fp12.x0.x2,&tmp3_fp2,&A->x0.x2);

	fp2_add(&tmp1_fp2,&tmp1_fp2,&A->x0.x1);
	fp2_mul_lazy(&tmp3_fp2,&A->x1.x2,&B->x1.x1);
	fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);
	fp2_add(&tmp3_fp12.x0.x1,&tmp1_fp2,&tmp3_fp2);

	fp2_mul_lazy(&tmp1_fp2,&A->x1.x2,&B->x1.x0);
	fp2_add(&tmp1_fp2,&tmp1_fp2,&tmp2_fp2);
	fp2_mul_basis(&tmp1_fp2,&tmp1_fp2);
	fp2_add(&tmp3_fp12.x0.x0,&tmp1_fp2,&A->x0.x0);

	fp12_set(ANS,&tmp3_fp12);	
}
void Pseudo_8_sparse_mul_lazy_montgomery(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2;
    static fp12_t tmp3_fp12;
	fp2_mul_lazy_montgomery(&tmp1_fp2,&A->x0.x0,&B->x1.x0);
	fp2_mul_lazy_montgomery(&tmp2_fp2,&A->x0.x1,&B->x1.x1);
	fp2_add_lazy(&tmp3_fp2,&A->x0.x0,&A->x0.x1);
	fp2_add_lazy(&tmp4_fp2,&B->x1.x0,&B->x1.x1);

	fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp1_fp2);		
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp2_fp2);
	fp2_add(&tmp3_fp12.x1.x1,&tmp3_fp2,&A->x1.x1);

	fp2_add(&tmp1_fp2,&tmp1_fp2,&A->x1.x0);
	fp2_mul_lazy_montgomery(&tmp3_fp2,&A->x0.x2,&B->x1.x1);
	fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);
	fp2_add(&tmp3_fp12.x1.x0,&tmp1_fp2,&tmp3_fp2);

	fp2_add(&tmp1_fp2,&tmp2_fp2,&A->x1.x2);
	fp2_mul_lazy_montgomery(&tmp2_fp2,&A->x0.x2,&B->x1.x0);
	fp2_add(&tmp3_fp12.x1.x2,&tmp1_fp2,&tmp2_fp2);

	fp2_mul_lazy_montgomery(&tmp1_fp2,&A->x1.x0,&B->x1.x0);//tmp1
	fp2_mul_lazy_montgomery(&tmp2_fp2,&A->x1.x1,&B->x1.x1);//tmp2

	fp2_add_lazy(&tmp3_fp2,&A->x1.x0,&A->x1.x1);
	fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp1_fp2);
	fp2_sub(&tmp3_fp2,&tmp3_fp2,&tmp2_fp2);
	fp2_add(&tmp3_fp12.x0.x2,&tmp3_fp2,&A->x0.x2);

	fp2_add(&tmp1_fp2,&tmp1_fp2,&A->x0.x1);
	fp2_mul_lazy_montgomery(&tmp3_fp2,&A->x1.x2,&B->x1.x1);
	fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);
	fp2_add(&tmp3_fp12.x0.x1,&tmp1_fp2,&tmp3_fp2);

	fp2_mul_lazy_montgomery(&tmp1_fp2,&A->x1.x2,&B->x1.x0);
	fp2_add(&tmp1_fp2,&tmp1_fp2,&tmp2_fp2);
	fp2_mul_basis(&tmp1_fp2,&tmp1_fp2);
	fp2_add(&tmp3_fp12.x0.x0,&tmp1_fp2,&A->x0.x0);

	fp12_set(ANS,&tmp3_fp12);	
}
/*
void fp12_6_sparse_mul_lazy(fp12_t *C, fp12_t *A, fp12_t *B)
{
	fp6_t tmp0_fp6, tmp1_fp6, tmp2_fp6;
	fp6_t a_fp6, b_fp6, c_fp6;
	fp2_t v0, v1, tmp0_fp2, tmp1_fp2, tmp2_fp2;

	fp2_mul_lazy(&tmp0_fp6.x0, &A->x0.x0, &B->x0.x0);
	fp2_mul_lazy(&tmp0_fp6.x1, &A->x0.x1, &B->x0.x0);
	fp2_mul_lazy(&tmp0_fp6.x2, &A->x0.x2, &B->x0.x0);
	fp2_add(&tmp2_fp6.x0, &B->x0.x0, &B->x1.x0);
	fp2_set(&tmp2_fp6.x1, &B->x1.x1);
	fp6_mul_dxs(&tmp1_fp6, &A->x1, &B->x1);

	fp6_set(&a_fp6,&A->x1);
	fp6_set(&b_fp6,&B->x1);
	////////
	fp2_mul_lazy(&v0, &a_fp6.x0, &b_fp6.x0);

	fp2_mul_lazy(&v1, &a_fp6.x1, &b_fp6.x1);

	fp2_add(&tmp0_fp2, &a_fp6.x1, &a_fp6.x2);
	fp2_mul_lazy(&tmp0_fp2, &tmp0_fp2, &b_fp6.x1);
	fp2_sub(&tmp0_fp2, &tmp0_fp2, &v1);
	fp2_mul_basis(&tmp2_fp2, &tmp0_fp2);
	fp2_add(&tmp2_fp2, &tmp2_fp2, &v0);

	fp2_add(&tmp0_fp2, &a_fp6.x0, &a_fp6.x1);
	fp2_add(&tmp1_fp2, &b_fp6.x0, &b_fp6.x1);
	fp2_mul_lazy(&c_fp6.x1, &tmp0_fp2, &tmp1_fp2);
	fp2_sub(&c_fp6.x1, &c_fp6.x1, &v0);
	fp2_sub(&c_fp6.x1, &c_fp6.x1, &v1);

	fp2_add(&tmp0_fp2, &a_fp6.x0, &a_fp6.x2);
	fp2_mul_lazy(&c_fp6.x2, &tmp0_fp2, &b_fp6.x0);
	fp2_sub(&c_fp6.x2, &c_fp6.x2, &v0);
	fp2_add(&c_fp6.x2, &c_fp6.x2, &v1);

	fp2_set(&c_fp6.x0, &tmp2_fp2);
	////////

	fp6_set(&tmp1_fp6,&c_fp6);
	
	fp6_add(&C->x1, &A->x0, &A->x1);
	fp6_mul_dxs(&C->x1, &C->x1, &tmp2_fp6);
	
	fp6_set(&a_fp6,&C->x1);
	fp6_set(&b_fp6,&tmp2_fp6);
	///////
	fp2_mul_lazy(&v0, &a_fp6.x0, &b_fp6.x0);

	fp2_mul_lazy(&v1, &a_fp6.x1, &b_fp6.x1);

	fp2_add(&tmp0_fp2, &a_fp6.x1, &a_fp6.x2);
	fp2_mul_lazy(&tmp0_fp2, &tmp0_fp2, &b_fp6.x1);
	fp2_sub(&tmp0_fp2, &tmp0_fp2, &v1);
	fp2_mul_basis(&tmp2_fp2, &tmp0_fp2);
	fp2_add(&tmp2_fp2, &tmp2_fp2, &v0);

	fp2_add(&tmp0_fp2, &a_fp6.x0, &a_fp6.x1);
	fp2_add(&tmp1_fp2, &b_fp6.x0, &b_fp6.x1);
	fp2_mul_lazy(&c_fp6.x1, &tmp0_fp2, &tmp1_fp2);
	fp2_sub(&c_fp6.x1, &c_fp6.x1, &v0);
	fp2_sub(&c_fp6.x1, &c_fp6.x1, &v1);

	fp2_add(&tmp0_fp2, &a_fp6.x0, &a_fp6.x2);
	fp2_mul_lazy(&c_fp6.x2, &tmp0_fp2, &b_fp6.x0);
	fp2_sub(&c_fp6.x2, &c_fp6.x2, &v0);
	fp2_add(&c_fp6.x2, &c_fp6.x2, &v1);

	fp2_set(&c_fp6.x0, &tmp2_fp2);
	////////
	fp6_set(&C->x1,&c_fp6);
	
	fp6_sub(&C->x1, &C->x1, &tmp0_fp6);
	fp6_sub(&C->x1, &C->x1, &tmp1_fp6);
	fp6_mul_basis(&tmp1_fp6, &tmp1_fp6);
	fp6_add(&C->x0, &tmp0_fp6, &tmp1_fp6);
}
*/

void fp12_6_sparse_mul_lazy(fp12_t *C, fp12_t *A, fp12_t *B)
{
	fp6_t tmp0_fp6, tmp1_fp6, tmp2_fp6;
	fp2_t v0, v1, v2, v3, v4;

	fp2_mul_lazy(&tmp0_fp6.x0, &A->x0.x0, &B->x0.x0);
	fp2_mul_lazy(&tmp0_fp6.x1, &A->x0.x1, &B->x0.x0);
	fp2_mul_lazy(&tmp0_fp6.x2, &A->x0.x2, &B->x0.x0);
	fp2_add(&tmp2_fp6.x0, &B->x0.x0, &B->x1.x0);
	fp2_set(&tmp2_fp6.x1, &B->x1.x1);

	//fp6_6_sparse_mul(&tmp1_fp6, &A->x1, &B->x1);
	fp2_mul_lazy(&v0, &A->x1.x0, &B->x1.x0);
	fp2_mul_lazy(&v1, &A->x1.x1, &B->x1.x1);
	fp2_add(&v2, &A->x1.x1, &A->x1.x2);
	fp2_mul_lazy(&v2, &v2, &B->x1.x1);
	fp2_sub(&v2, &v2, &v1);
	fp2_mul_basis(&v4, &v2);
	fp2_add(&v4, &v4, &v0);
	fp2_add(&v2, &A->x1.x0, &A->x1.x1);
	fp2_add(&v3, &B->x1.x0, &B->x1.x1);
	fp2_mul_lazy(&tmp1_fp6.x1, &v2, &v3);
	fp2_sub(&tmp1_fp6.x1, &tmp1_fp6.x1, &v0);
	fp2_sub(&tmp1_fp6.x1, &tmp1_fp6.x1, &v1);
	fp2_add(&v2, &A->x1.x0, &A->x1.x2);
	fp2_mul_lazy(&tmp1_fp6.x2, &v2, &B->x1.x0);
	fp2_sub(&tmp1_fp6.x2, &tmp1_fp6.x2, &v0);
	fp2_add(&tmp1_fp6.x2, &tmp1_fp6.x2, &v1);
	fp2_set(&tmp1_fp6.x0, &v4);
	//fp6_sparse_mul end

	fp6_add(&C->x1, &A->x0, &A->x1);

	//fp6_6_sparse_mul(&C->x1, &C->x1, &tmp2_fp6);
	fp2_mul_lazy(&v0, &C->x1.x0, &tmp2_fp6.x0);
	fp2_mul_lazy(&v1, &C->x1.x1, &tmp2_fp6.x1);
	fp2_add(&v2, &C->x1.x1, &C->x1.x2);
	fp2_mul_lazy(&v2, &v2, &tmp2_fp6.x1);
	fp2_sub(&v2, &v2, &v1);
	fp2_mul_basis(&v4, &v2);
	fp2_add(&v4, &v4, &v0);
	fp2_add(&v2, &C->x1.x0, &C->x1.x1);
	fp2_add(&v3, &tmp2_fp6.x0, &tmp2_fp6.x1);
	fp2_mul_lazy(&C->x1.x1, &v2, &v3);
	fp2_sub(&C->x1.x1, &C->x1.x1, &v0);
	fp2_sub(&C->x1.x1, &C->x1.x1, &v1);
	fp2_add(&v2, &C->x1.x0, &C->x1.x2);
	fp2_mul_lazy(&C->x1.x2, &v2, &tmp2_fp6.x0);
	fp2_sub(&C->x1.x2, &C->x1.x2, &v0);
	fp2_add(&C->x1.x2, &C->x1.x2, &v1);
	fp2_set(&C->x1.x0, &v4);
  	//fp6_sparse_mul end

	fp6_sub(&C->x1, &C->x1, &tmp0_fp6);
	fp6_sub(&C->x1, &C->x1, &tmp1_fp6);
	fp6_mul_basis(&tmp1_fp6, &tmp1_fp6);
	fp6_add(&C->x0, &tmp0_fp6, &tmp1_fp6);
}

void fp12_6_sparse_mul_lazy_montgomery(fp12_t *C, fp12_t *A, fp12_t *B)
{
	fp6_t tmp0_fp6, tmp1_fp6, tmp2_fp6;
	fp2_t v0, v1, v2, v3, v4;

	fp2_mul_lazy_montgomery(&tmp0_fp6.x0, &A->x0.x0, &B->x0.x0);
	fp2_mul_lazy_montgomery(&tmp0_fp6.x1, &A->x0.x1, &B->x0.x0);
	fp2_mul_lazy_montgomery(&tmp0_fp6.x2, &A->x0.x2, &B->x0.x0);
	fp2_add(&tmp2_fp6.x0, &B->x0.x0, &B->x1.x0);
	fp2_set(&tmp2_fp6.x1, &B->x1.x1);

	//fp6_6_sparse_mul(&tmp1_fp6, &A->x1, &B->x1);
	fp2_mul_lazy_montgomery(&v0, &A->x1.x0, &B->x1.x0);
	fp2_mul_lazy_montgomery(&v1, &A->x1.x1, &B->x1.x1);
	fp2_add(&v2, &A->x1.x1, &A->x1.x2);
	fp2_mul_lazy_montgomery(&v2, &v2, &B->x1.x1);
	fp2_sub(&v2, &v2, &v1);
	fp2_mul_basis(&v4, &v2);
	fp2_add(&v4, &v4, &v0);
	fp2_add(&v2, &A->x1.x0, &A->x1.x1);
	fp2_add(&v3, &B->x1.x0, &B->x1.x1);
	fp2_mul_lazy_montgomery(&tmp1_fp6.x1, &v2, &v3);
	fp2_sub(&tmp1_fp6.x1, &tmp1_fp6.x1, &v0);
	fp2_sub(&tmp1_fp6.x1, &tmp1_fp6.x1, &v1);
	fp2_add(&v2, &A->x1.x0, &A->x1.x2);
	fp2_mul_lazy_montgomery(&tmp1_fp6.x2, &v2, &B->x1.x0);
	fp2_sub(&tmp1_fp6.x2, &tmp1_fp6.x2, &v0);
	fp2_add(&tmp1_fp6.x2, &tmp1_fp6.x2, &v1);
	fp2_set(&tmp1_fp6.x0, &v4);
	//fp6_sparse_mul end

	fp6_add(&C->x1, &A->x0, &A->x1);

	//fp6_6_sparse_mul(&C->x1, &C->x1, &tmp2_fp6);
	fp2_mul_lazy_montgomery(&v0, &C->x1.x0, &tmp2_fp6.x0);
	fp2_mul_lazy_montgomery(&v1, &C->x1.x1, &tmp2_fp6.x1);
	fp2_add(&v2, &C->x1.x1, &C->x1.x2);
	fp2_mul_lazy_montgomery(&v2, &v2, &tmp2_fp6.x1);
	fp2_sub(&v2, &v2, &v1);
	fp2_mul_basis(&v4, &v2);
	fp2_add(&v4, &v4, &v0);
	fp2_add(&v2, &C->x1.x0, &C->x1.x1);
	fp2_add(&v3, &tmp2_fp6.x0, &tmp2_fp6.x1);
	fp2_mul_lazy_montgomery(&C->x1.x1, &v2, &v3);
	fp2_sub(&C->x1.x1, &C->x1.x1, &v0);
	fp2_sub(&C->x1.x1, &C->x1.x1, &v1);
	fp2_add(&v2, &C->x1.x0, &C->x1.x2);
	fp2_mul_lazy_montgomery(&C->x1.x2, &v2, &tmp2_fp6.x0);
	fp2_sub(&C->x1.x2, &C->x1.x2, &v0);
	fp2_add(&C->x1.x2, &C->x1.x2, &v1);
	fp2_set(&C->x1.x0, &v4);
  	//fp6_sparse_mul end

	fp6_sub(&C->x1, &C->x1, &tmp0_fp6);
	fp6_sub(&C->x1, &C->x1, &tmp1_fp6);
	fp6_mul_basis(&tmp1_fp6, &tmp1_fp6);
	fp6_add(&C->x0, &tmp0_fp6, &tmp1_fp6);
}
void ff_ltt(fp12_t *f,efp2_t *T,efp_t *P,fp_t *L){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;
    efp2_t Tmp_T;
	efp2_init(&Tmp_T);
	efp2_set(&Tmp_T,T);
	
	fp12_sqr(&tmp1_fp12,f);
	
	//ltt
	fp2_add(&tmp1_fp2,&Tmp_T.y,&Tmp_T.y);		//tmp1_fp2=1/(2*T.y)
	fp2_inv(&tmp1_fp2,&tmp1_fp2);
	fp2_sqr(&tmp2_fp2,&Tmp_T.x);			//tmp2_fp2=3(T.x)^2
	fp2_add(&tmp3_fp2,&tmp2_fp2,&tmp2_fp2);
	fp2_add(&tmp2_fp2,&tmp3_fp2,&tmp2_fp2);
	fp2_mul(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);				//tmp3_fp2=tmp1_fp2*tmp2_fp2
	
	fp2_add(&tmp4_fp2,&Tmp_T.x,&Tmp_T.x);		//tmp4_fp2=2T.x
	fp2_sqr(&T->x,&tmp3_fp2);				//next_T.x=tmp3_fp2^2-tmp4_fp2
	fp2_sub(&T->x,&T->x,&tmp4_fp2);
	fp2_mul(&tmp5_fp2,&tmp3_fp2,&Tmp_T.x);			//tmp5_fp2=tmp3_fp2*T.x-T.y
	fp2_sub(&tmp5_fp2,&tmp5_fp2,&Tmp_T.y);
	fp2_mul(&T->y,&tmp3_fp2,&T->x);			//next_T.y=tmp5_fp2-tmp3_fp2*next_T.x
	fp2_sub(&T->y,&tmp5_fp2,&T->y);
	
	//set ltt
	fp_set_ui(&tmp2_fp12.x0.x0.x0,1);
	fp2_set_neg(&tmp2_fp12.x1.x0,&tmp3_fp2);
	fp2_mul_mpn(&tmp2_fp12.x1.x1,&tmp5_fp2,L->x0);
	
	Pseudo_8_sparse_mul(f,&tmp1_fp12,&tmp2_fp12);
}
void ff_ltt_lazy(fp12_t *f,efp2_t *T,efp_t *P,fp_t *L){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;
    efp2_t Tmp_T;
	efp2_init(&Tmp_T);
	efp2_set(&Tmp_T,T);
	
	fp12_sqr_lazy(&tmp1_fp12,f);
	
	//ltt
	fp2_add(&tmp1_fp2,&Tmp_T.y,&Tmp_T.y);		//tmp1_fp2=1/(2*T.y)
	fp2_inv_lazy(&tmp1_fp2,&tmp1_fp2);
	fp2_sqr_lazy(&tmp2_fp2,&Tmp_T.x);			//tmp2_fp2=3(T.x)^2
	fp2_add_lazy(&tmp3_fp2,&tmp2_fp2,&tmp2_fp2);
	fp2_add_lazy(&tmp2_fp2,&tmp3_fp2,&tmp2_fp2);
	fp2_mul_lazy(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);		//tmp3_fp2=tmp1_fp2*tmp2_fp2
	
	fp2_add(&tmp4_fp2,&Tmp_T.x,&Tmp_T.x);		//tmp4_fp2=2T.x
	fp2_sqr_lazy(&T->x,&tmp3_fp2);				//next_T.x=tmp3_fp2^2-tmp4_fp2
	fp2_sub(&T->x,&T->x,&tmp4_fp2);
	fp2_mul_lazy(&tmp5_fp2,&tmp3_fp2,&Tmp_T.x);		//tmp5_fp2=tmp3_fp2*T.x-T.y
	fp2_sub(&tmp5_fp2,&tmp5_fp2,&Tmp_T.y);
	fp2_mul_lazy(&T->y,&tmp3_fp2,&T->x);			//next_T.y=tmp5_fp2-tmp3_fp2*next_T.x
	fp2_sub(&T->y,&tmp5_fp2,&T->y);
	
	//set ltt
	fp_set_ui(&tmp2_fp12.x0.x0.x0,1);
	fp2_set_neg(&tmp2_fp12.x1.x0,&tmp3_fp2);
	fp2_mul_mpn(&tmp2_fp12.x1.x1,&tmp5_fp2,L->x0);

	Pseudo_8_sparse_mul_lazy(f,&tmp1_fp12,&tmp2_fp12);

}
void ff_ltt_projective_lazy(fp12_t *f,efp2_projective_t *T,efp_t *P){
    static fp2_t tmp0_fp2,tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;
    efp2_projective_t Tmp_T;
	efp2_projective_init(&Tmp_T);
	efp2_projective_set(&Tmp_T,T);
	
	static fp2_t u0, u1;

	fp12_sqr_lazy(&tmp1_fp12, f);

	fp2_sqr_lazy(&tmp0_fp2, &T->z);
	fp2_sqr_lazy(&tmp1_fp2, &T->y);
	fp2_add(&tmp5_fp2, &tmp0_fp2, &tmp1_fp2);
	fp2_mul_ui(&tmp0_fp2, &tmp0_fp2,24);
	fp_add(&tmp2_fp2.x0,&tmp0_fp2.x0,&tmp0_fp2.x1);
	fp_sub(&tmp2_fp2.x1,&tmp0_fp2.x1,&tmp0_fp2.x0);

	fp2_sqr_lazy(&tmp0_fp2, &T->x);
	fp2_mul_lazy(&tmp4_fp2, &T->x, &T->y);
	fp_hlv(&tmp4_fp2.x0, &tmp4_fp2.x0);
	fp_hlv(&tmp4_fp2.x1, &tmp4_fp2.x1);
	fp2_dbl(&tmp3_fp2, &tmp2_fp2);
	fp2_add(&tmp3_fp2, &tmp3_fp2, &tmp2_fp2);
	fp2_sub(&T->x, &tmp1_fp2, &tmp3_fp2);
	fp2_mul_lazy(&T->x, &T->x, &tmp4_fp2);

	fp2_add(&tmp3_fp2, &tmp1_fp2, &tmp3_fp2);
	fp_hlv(&tmp3_fp2.x0, &tmp3_fp2.x0);
	fp_hlv(&tmp3_fp2.x1, &tmp3_fp2.x1);

	fp2_sqr_lazy(&u0, &tmp2_fp2);
	fp2_add(&u1, &u0, &u0);
	fp2_add(&u1, &u1, &u0);
	fp2_sqr_lazy(&u0, &tmp3_fp2);
	fp2_sub(&u0, &u0, &u1);


	fp2_add(&tmp3_fp2, &T->y, &T->z);
	fp2_sqr_lazy(&tmp3_fp2, &tmp3_fp2);
	fp2_sub(&tmp3_fp2, &tmp3_fp2, &tmp5_fp2);

	fp2_set(&T->y, &u0);

	fp2_mul_lazy(&T->z, &tmp1_fp2, &tmp3_fp2);

	fp2_sub(&tmp2_fp12.x1.x1, &tmp2_fp2, &tmp1_fp2);

	fp_mul(&tmp2_fp12.x1.x0.x0, &P->x, &tmp0_fp2.x0);
	fp_mul(&tmp2_fp12.x1.x0.x1, &P->x, &tmp0_fp2.x1);
	fp2_mul_ui(&tmp2_fp12.x1.x0, &tmp2_fp12.x1.x0, 3);

	fp_mul(&tmp2_fp12.x0.x0.x0, &tmp3_fp2.x0, &P->y);
	fp_mul(&tmp2_fp12.x0.x0.x1, &tmp3_fp2.x1, &P->y);
	fp2_set_neg(&tmp2_fp12.x0.x0, &tmp2_fp12.x0.x0);


	//fp12_mul_lazy(f, &tmp1_fp12, &tmp2_fp12);
	fp12_6_sparse_mul_lazy(f, &tmp1_fp12, &tmp2_fp12);

}
void ff_ltt_projective_lazy_montgomery(fp12_t *f,efp2_projective_t *T,efp_t *P){
    static fp2_t tmp0_fp2,tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,buf_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;
    efp2_projective_t Tmp_T;
	efp2_projective_init(&Tmp_T);
	efp2_projective_set(&Tmp_T,T);
	
	static fp2_t u0, u1;

	fp12_sqr_lazy_montgomery(&tmp1_fp12, f);

	fp2_sqr_lazy_montgomery(&tmp0_fp2, &T->z);
	fp2_sqr_lazy_montgomery(&tmp1_fp2, &T->y);
	fp2_add(&tmp5_fp2, &tmp0_fp2, &tmp1_fp2);
	
	fp2_mul_ui(&tmp0_fp2, &tmp0_fp2,24);
	
	fp_add(&tmp2_fp2.x0,&tmp0_fp2.x0,&tmp0_fp2.x1);
	fp_sub(&tmp2_fp2.x1,&tmp0_fp2.x1,&tmp0_fp2.x0);

	fp2_sqr_lazy_montgomery(&tmp0_fp2, &T->x);

	fp2_mul_lazy_montgomery(&tmp4_fp2, &T->x, &T->y);
	fp_hlv(&tmp4_fp2.x0, &tmp4_fp2.x0);
	fp_hlv(&tmp4_fp2.x1, &tmp4_fp2.x1);
	fp2_dbl(&tmp3_fp2, &tmp2_fp2);
	fp2_add(&tmp3_fp2, &tmp3_fp2, &tmp2_fp2);
	fp2_sub(&T->x, &tmp1_fp2, &tmp3_fp2);
	fp2_mul_lazy_montgomery(&T->x, &T->x, &tmp4_fp2);

	fp2_add(&tmp3_fp2, &tmp1_fp2, &tmp3_fp2);
	fp_hlv(&tmp3_fp2.x0, &tmp3_fp2.x0);
	fp_hlv(&tmp3_fp2.x1, &tmp3_fp2.x1);

	fp2_sqr_lazy_montgomery(&u0, &tmp2_fp2);
	fp2_add(&u1, &u0, &u0);
	fp2_add(&u1, &u1, &u0);
	fp2_sqr_lazy_montgomery(&u0, &tmp3_fp2);
	fp2_sub(&u0, &u0, &u1);


	fp2_add(&tmp3_fp2, &T->y, &T->z);
	fp2_sqr_lazy_montgomery(&tmp3_fp2, &tmp3_fp2);
	fp2_sub(&tmp3_fp2, &tmp3_fp2, &tmp5_fp2);

	fp2_set(&T->y, &u0);

	fp2_mul_lazy_montgomery(&T->z, &tmp1_fp2, &tmp3_fp2);

	fp2_sub(&tmp2_fp12.x1.x1, &tmp2_fp2, &tmp1_fp2);

	fp_mulmod_montgomery(&tmp2_fp12.x1.x0.x0, &P->x, &tmp0_fp2.x0);
	fp_mulmod_montgomery(&tmp2_fp12.x1.x0.x1, &P->x, &tmp0_fp2.x1);
	//fp2_mul_ui(&tmp2_fp12.x1.x0, &tmp2_fp12.x1.x0, 3);
	fp2_add(&buf_fp2,&tmp2_fp12.x1.x0, &tmp2_fp12.x1.x0);
	fp2_add(&tmp2_fp12.x1.x0,&buf_fp2, &tmp2_fp12.x1.x0);

	fp_mulmod_montgomery(&tmp2_fp12.x0.x0.x0, &tmp3_fp2.x0, &P->y);
	fp_mulmod_montgomery(&tmp2_fp12.x0.x0.x1, &tmp3_fp2.x1, &P->y);
	fp2_set_neg(&tmp2_fp12.x0.x0, &tmp2_fp12.x0.x0);


	//fp12_mul_lazy_montgomery(f, &tmp1_fp12, &tmp2_fp12);
	fp12_6_sparse_mul_lazy_montgomery(f, &tmp1_fp12, &tmp2_fp12);

}
void ff_ltt_lazy_montgomery(fp12_t *f,efp2_t *T,efp_t *P,fp_t *L){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;
    efp2_t Tmp_T;
	efp2_init(&Tmp_T);
	efp2_set(&Tmp_T,T);
	//TODO:lazy->monty,f seted 1 or RmodP
	fp12_sqr_lazy(&tmp1_fp12,f);
	
	//ltt
	fp2_add(&tmp1_fp2,&Tmp_T.y,&Tmp_T.y);		//tmp1_fp2=1/(2*T.y)
	fp2_inv_lazy_montgomery(&tmp1_fp2,&tmp1_fp2);
	fp2_sqr_lazy_montgomery(&tmp2_fp2,&Tmp_T.x);			//tmp2_fp2=3(T.x)^2
	fp2_add_lazy(&tmp3_fp2,&tmp2_fp2,&tmp2_fp2);
	fp2_add_lazy(&tmp2_fp2,&tmp3_fp2,&tmp2_fp2);
	fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);		//tmp3_fp2=tmp1_fp2*tmp2_fp2
	
	fp2_add(&tmp4_fp2,&Tmp_T.x,&Tmp_T.x);		//tmp4_fp2=2T.x
	fp2_sqr_lazy_montgomery(&T->x,&tmp3_fp2);				//next_T.x=tmp3_fp2^2-tmp4_fp2
	fp2_sub(&T->x,&T->x,&tmp4_fp2);
	fp2_mul_lazy_montgomery(&tmp5_fp2,&tmp3_fp2,&Tmp_T.x);		//tmp5_fp2=tmp3_fp2*T.x-T.y
	fp2_sub(&tmp5_fp2,&tmp5_fp2,&Tmp_T.y);
	fp2_mul_lazy_montgomery(&T->y,&tmp3_fp2,&T->x);			//next_T.y=tmp5_fp2-tmp3_fp2*next_T.x
	fp2_sub(&T->y,&tmp5_fp2,&T->y);
	
	//set ltt
	//TODO:1->RmodP?
	fp_set_ui(&tmp2_fp12.x0.x0.x0,1);
	//fp_set_mpn(&tmp2_fp12.x0.x0.x0,RmodP);
	fp2_set_neg(&tmp2_fp12.x1.x0,&tmp3_fp2);
	fp2_mul_mpn_montgomery(&tmp2_fp12.x1.x1,&tmp5_fp2,L->x0);

	Pseudo_8_sparse_mul_lazy_montgomery(f,&tmp1_fp12,&tmp2_fp12);
}
void f_ltq(fp12_t *f,efp2_t *T,efp2_t *Q,efp_t *P,fp_t *L){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12;
    efp2_t Tmp_T;
	efp2_init(&Tmp_T);
	efp2_set(&Tmp_T,T);
		
	//ltq
	fp2_sub(&tmp1_fp2,&Q->x,&Tmp_T.x);		//tmp1_fp2=(Q->x-T.x)^-1
	fp2_inv(&tmp1_fp2,&tmp1_fp2);
	fp2_sub(&tmp2_fp2,&Q->y,&Tmp_T.y);		//tmp2_fp2=(Q->y-T.y)
	fp2_mul(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);			//tmp3_fp2=tmp1_fp2*tmp2_fp2
	fp2_add(&tmp4_fp2,&Tmp_T.x,&Q->x);		//tmp4_fp2=Q->x+T.x
	fp2_sqr(&T->x,&tmp3_fp2);			//next_T.x=tmp3_fp2^2-tmp4_fp2
	fp2_sub(&T->x,&T->x,&tmp4_fp2);
	fp2_mul(&tmp5_fp2,&tmp3_fp2,&Tmp_T.x);		//tmp5_fp2=tmp3_fp2*T.x-T.y
	fp2_sub(&tmp5_fp2,&tmp5_fp2,&Tmp_T.y);
	fp2_mul(&T->y,&tmp3_fp2,&T->x);		//next_T.y=tmp5_fp2-tmp3_fp2*next_T.x
	fp2_sub(&T->y,&tmp5_fp2,&T->y);
	
	//set ltq
	fp_set_ui(&tmp1_fp12.x0.x0.x0,1);
	fp2_set_neg(&tmp1_fp12.x1.x0,&tmp3_fp2);
	fp2_mul_mpn(&tmp1_fp12.x1.x1,&tmp5_fp2,L->x0);
	
	Pseudo_8_sparse_mul(f,f,&tmp1_fp12);
}
void f_ltq_lazy(fp12_t *f,efp2_t *T,efp2_t *Q,efp_t *P,fp_t *L){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12;
    efp2_t Tmp_T;
	efp2_init(&Tmp_T);
	efp2_set(&Tmp_T,T);
	
	//ltq
	fp2_sub(&tmp1_fp2,&Q->x,&Tmp_T.x);		//tmp1_fp2=(Q->x-T.x)^-1
	fp2_inv_lazy(&tmp1_fp2,&tmp1_fp2);
	fp2_sub_lazy(&tmp2_fp2,&Q->y,&Tmp_T.y);		//tmp2_fp2=(Q->y-T.y)
	fp2_mul_lazy(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);			//tmp3_fp2=tmp1_fp2*tmp2_fp2
	fp2_add(&tmp4_fp2,&Tmp_T.x,&Q->x);		//tmp4_fp2=Q->x+T.x
	fp2_sqr_lazy(&T->x,&tmp3_fp2);			//next_T.x=tmp3_fp2^2-tmp4_fp2
	fp2_sub(&T->x,&T->x,&tmp4_fp2);
	fp2_mul_lazy(&tmp5_fp2,&tmp3_fp2,&Tmp_T.x);		//tmp5_fp2=tmp3_fp2*T.x-T.y
	fp2_sub(&tmp5_fp2,&tmp5_fp2,&Tmp_T.y);
	fp2_mul_lazy(&T->y,&tmp3_fp2,&T->x);		//next_T.y=tmp5_fp2-tmp3_fp2*next_T.x
	fp2_sub(&T->y,&tmp5_fp2,&T->y);
	
	//set ltq
	fp_set_ui(&tmp1_fp12.x0.x0.x0,1);
	fp2_set_neg(&tmp1_fp12.x1.x0,&tmp3_fp2);
	fp2_mul_mpn(&tmp1_fp12.x1.x1,&tmp5_fp2,L->x0);
	
	Pseudo_8_sparse_mul_lazy(f,f,&tmp1_fp12);

}
void f_ltq_projective_lazy(fp12_t *f,efp2_projective_t *T,efp2_projective_t *Q,efp_t *P){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12;
	
	fp2_mul_lazy(&tmp1_fp2, &T->z, &Q->x);
	fp2_sub(&tmp1_fp2, &T->x, &tmp1_fp2);
	fp2_mul_lazy(&tmp2_fp2, &T->z, &Q->y);
	fp2_sub(&tmp2_fp2, &T->y, &tmp2_fp2);

	fp2_sqr_lazy(&tmp3_fp2, &tmp1_fp2);
	fp2_mul_lazy(&T->x, &T->x, &tmp3_fp2);
	fp2_mul_lazy(&tmp3_fp2, &tmp3_fp2, &tmp1_fp2);
	fp2_sqr_lazy(&tmp4_fp2, &tmp2_fp2);
	fp2_mul_lazy(&tmp4_fp2, &tmp4_fp2, &T->z);
	fp2_add(&tmp4_fp2, &tmp3_fp2, &tmp4_fp2);

	fp_mul(&tmp1_fp12.x1.x0.x0, &tmp2_fp2.x0, &P->x);
	fp_mul(&tmp1_fp12.x1.x0.x1, &tmp2_fp2.x1, &P->x);
	fp2_set_neg(&tmp1_fp12.x1.x0, &tmp1_fp12.x1.x0);

	fp2_mul_lazy(&tmp5_fp2, &Q->x, &tmp2_fp2);

	fp2_sub(&tmp4_fp2, &tmp4_fp2, &T->x);
	fp2_sub(&tmp4_fp2, &tmp4_fp2, &T->x);
	fp2_sub(&T->x, &T->x, &tmp4_fp2);
	fp2_mul_lazy(&tmp2_fp2, &tmp2_fp2, &T->x);
	fp2_mul_lazy(&T->y, &tmp3_fp2, &T->y);
	fp2_sub(&T->y, &tmp2_fp2, &T->y);
	fp2_mul_lazy(&T->x, &tmp1_fp2, &tmp4_fp2);
	fp2_mul_lazy(&T->z, &T->z, &tmp3_fp2);

	fp2_mul_lazy(&tmp3_fp2, &Q->y, &tmp1_fp2);
	fp2_sub(&tmp1_fp12.x1.x1, &tmp5_fp2, &tmp3_fp2);

	fp_mul(&tmp1_fp12.x0.x0.x0, &tmp1_fp2.x0, &P->y);
	fp_mul(&tmp1_fp12.x0.x0.x1, &tmp1_fp2.x1, &P->y);
	
	fp12_mul_lazy(f,&tmp1_fp12,f);
	//fp12_6_sparse_mul_lazy(f,&tmp1_fp12,f);

}
void f_ltq_lazy_montgomery(fp12_t *f,efp2_t *T,efp2_t *Q,efp_t *P,fp_t *L){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12;
    efp2_t Tmp_T;
	efp2_init(&Tmp_T);
	efp2_set(&Tmp_T,T);
	
	//ltq
	fp2_sub(&tmp1_fp2,&Q->x,&Tmp_T.x);		//tmp1_fp2=(Q->x-T.x)^-1
	fp2_inv_lazy_montgomery(&tmp1_fp2,&tmp1_fp2);
	fp2_sub_lazy(&tmp2_fp2,&Q->y,&Tmp_T.y);		//tmp2_fp2=(Q->y-T.y)
	fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);			//tmp3_fp2=tmp1_fp2*tmp2_fp2
	fp2_add(&tmp4_fp2,&Tmp_T.x,&Q->x);		//tmp4_fp2=Q->x+T.x
	fp2_sqr_lazy_montgomery(&T->x,&tmp3_fp2);			//next_T.x=tmp3_fp2^2-tmp4_fp2
	fp2_sub(&T->x,&T->x,&tmp4_fp2);
	fp2_mul_lazy_montgomery(&tmp5_fp2,&tmp3_fp2,&Tmp_T.x);		//tmp5_fp2=tmp3_fp2*T.x-T.y
	fp2_sub(&tmp5_fp2,&tmp5_fp2,&Tmp_T.y);
	fp2_mul_lazy_montgomery(&T->y,&tmp3_fp2,&T->x);		//next_T.y=tmp5_fp2-tmp3_fp2*next_T.x
	fp2_sub(&T->y,&tmp5_fp2,&T->y);
	
	//set ltq
	fp_set_ui(&tmp1_fp12.x0.x0.x0,1);
	fp2_set_neg(&tmp1_fp12.x1.x0,&tmp3_fp2);
	fp2_mul_mpn_montgomery(&tmp1_fp12.x1.x1,&tmp5_fp2,L->x0);
	
	Pseudo_8_sparse_mul_lazy_montgomery(f,f,&tmp1_fp12);

}
void f_ltq_projective_lazy_montgomery(fp12_t *f,efp2_projective_t *T,efp2_projective_t *Q,efp_t *P){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12;
	
	fp2_mul_lazy_montgomery(&tmp1_fp2, &T->z, &Q->x);
	fp2_sub(&tmp1_fp2, &T->x, &tmp1_fp2);
	fp2_mul_lazy_montgomery(&tmp2_fp2, &T->z, &Q->y);
	fp2_sub(&tmp2_fp2, &T->y, &tmp2_fp2);

	fp2_sqr_lazy_montgomery(&tmp3_fp2, &tmp1_fp2);
	fp2_mul_lazy_montgomery(&T->x, &T->x, &tmp3_fp2);
	fp2_mul_lazy_montgomery(&tmp3_fp2, &tmp3_fp2, &tmp1_fp2);
	fp2_sqr_lazy_montgomery(&tmp4_fp2, &tmp2_fp2);
	fp2_mul_lazy_montgomery(&tmp4_fp2, &tmp4_fp2, &T->z);
	fp2_add(&tmp4_fp2, &tmp3_fp2, &tmp4_fp2);

	fp_mulmod_montgomery(&tmp1_fp12.x1.x0.x0, &tmp2_fp2.x0, &P->x);
	fp_mulmod_montgomery(&tmp1_fp12.x1.x0.x1, &tmp2_fp2.x1, &P->x);
	fp2_set_neg(&tmp1_fp12.x1.x0, &tmp1_fp12.x1.x0);

	fp2_mul_lazy_montgomery(&tmp5_fp2, &Q->x, &tmp2_fp2);

	fp2_sub(&tmp4_fp2, &tmp4_fp2, &T->x);
	fp2_sub(&tmp4_fp2, &tmp4_fp2, &T->x);
	fp2_sub(&T->x, &T->x, &tmp4_fp2);
	fp2_mul_lazy_montgomery(&tmp2_fp2, &tmp2_fp2, &T->x);
	fp2_mul_lazy_montgomery(&T->y, &tmp3_fp2, &T->y);
	fp2_sub(&T->y, &tmp2_fp2, &T->y);
	fp2_mul_lazy_montgomery(&T->x, &tmp1_fp2, &tmp4_fp2);
	fp2_mul_lazy_montgomery(&T->z, &T->z, &tmp3_fp2);

	fp2_mul_lazy_montgomery(&tmp3_fp2, &Q->y, &tmp1_fp2);
	fp2_sub(&tmp1_fp12.x1.x1, &tmp5_fp2, &tmp3_fp2);

	fp_mulmod_montgomery(&tmp1_fp12.x0.x0.x0, &tmp1_fp2.x0, &P->y);
	fp_mulmod_montgomery(&tmp1_fp12.x0.x0.x1, &tmp1_fp2.x1, &P->y);
	
	fp12_mul_lazy_montgomery(f,&tmp1_fp12,f);
	//fp12_6_sparse_mul_lazy_montgomery(f,&tmp1_fp12,f);
}