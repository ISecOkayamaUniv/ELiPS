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
void ff_ltt_lazy_montgomery(fp12_t *f,efp2_t *T,efp_t *P,fp_t *L){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;
    efp2_t Tmp_T;
	efp2_init(&Tmp_T);
	efp2_set(&Tmp_T,T);
	
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
	fp_set_ui(&tmp2_fp12.x0.x0.x0,1);
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
