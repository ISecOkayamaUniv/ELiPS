#include <ELiPS/sparse.h>
//Pseudo 8-sparse
void Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L){
    EFp2 Tmp_Q;
	EFp2_init(&Tmp_Q);
	EFp Tmp_P;
	EFp_init(&Tmp_P);
	Fp A,B,C,D,c;
	Fp_init(&A);
	Fp_init(&B);
	Fp_init(&C);
	Fp_init(&D);
	Fp_init(&c);
	
	EFp_set(&Tmp_P,P);
	EFp2_set(&Tmp_Q,Q);
	
	Fp_mul(&A,&Tmp_P.x,&Tmp_P.y);
	Fp_inv(&A,&A);
	Fp_mul(&B,&Tmp_P.x,&Tmp_P.x);
	Fp_mul(&B,&B,&A);
	Fp_mul(&C,&Tmp_P.y,&A);
	Fp_mul(&D,&B,&B);
	
	Fp2_mul_mpn(&Q->x,&Tmp_Q.x,D.x0);
	Fp_mul(&c,&B,&D);
	Fp2_mul_mpn(&Q->y,&Tmp_Q.y,c.x0);
	
	Fp_mul(&P->x,&D,&Tmp_P.x);
	Fp_set(&P->y,&P->x);
	
	Fp_mul(L,&C,&Tmp_P.y);
	Fp_mul(L,L,L);
	Fp_mul(L,L,&C);
}
void Pseudo_8_sparse_mapping_montgomery(EFp *P,EFp2 *Q,Fp *L){
    EFp2 Tmp_Q;
	EFp2_init(&Tmp_Q);
	EFp Tmp_P;
	EFp_init(&Tmp_P);
	Fp A,B,C,D,c;
	Fp_init(&A);
	Fp_init(&B);
	Fp_init(&C);
	Fp_init(&D);
	Fp_init(&c);
	
	EFp_set(&Tmp_P,P);
	EFp2_set(&Tmp_Q,Q);
	
	Fp_mulmod_montgomery(&A,&Tmp_P.x,&Tmp_P.y);
	Fp_inv_montgomery(&A,&A);
	Fp_mulmod_montgomery(&B,&Tmp_P.x,&Tmp_P.x);
	Fp_mulmod_montgomery(&B,&B,&A);
	Fp_mulmod_montgomery(&C,&Tmp_P.y,&A);
	Fp_mulmod_montgomery(&D,&B,&B);
	
	Fp2_mul_mpn_montgomery(&Q->x,&Tmp_Q.x,D.x0);
	Fp_mulmod_montgomery(&c,&B,&D);
	Fp2_mul_mpn_montgomery(&Q->y,&Tmp_Q.y,c.x0);
	
	Fp_mulmod_montgomery(&P->x,&D,&Tmp_P.x);
	Fp_set(&P->y,&P->x);
	
	Fp_mulmod_montgomery(L,&C,&Tmp_P.y);
	Fp_mulmod_montgomery(L,L,L);
	Fp_mulmod_montgomery(L,L,&C);
}

void Pseudo_8_sparse_mul(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2;
    static Fp12 tmp3_Fp12;
	Fp2_mul(&tmp1_Fp2,&A->x0.x0,&B->x1.x0);		//tmp1_Fp2=b3*f0
	Fp2_mul(&tmp2_Fp2,&A->x0.x1,&B->x1.x1);		//tmp2_Fp2=b4*f1
	Fp2_add(&tmp3_Fp2,&A->x0.x0,&A->x0.x1);		//tmp3_Fp2=f0+f1
	Fp2_add(&tmp4_Fp2,&B->x1.x0,&B->x1.x1);		//tmp4_Fp2=b3+b4

	Fp2_mul(&tmp3_Fp2,&tmp3_Fp2,&tmp4_Fp2);			//tmp3_Fp2=tmp3_Fp2*tmp4_Fp2
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp1_Fp2);			//tmp3_Fp2=tmp3_Fp2-tmp1_Fp2
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp2_Fp2);			//tmp3_Fp2=tmp3_Fp2-tmp2_Fp2

	Fp2_add(&tmp3_Fp12.x1.x1,&tmp3_Fp2,&A->x1.x1);	//ans[γ^3]=tmp3_Fp2+f4
	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&A->x1.x0);		//tmp1_Fp2=tmp1_Fp2+f3
	Fp2_mul(&tmp3_Fp2,&A->x0.x2,&B->x1.x1);		//tmp3_Fp2=b4*f2
	Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);

	Fp2_add(&tmp3_Fp12.x1.x0,&tmp1_Fp2,&tmp3_Fp2);		//ans[γ]=tmp1_Fp2+tmp3_Fp2
	Fp2_add(&tmp1_Fp2,&tmp2_Fp2,&A->x1.x2);		//tmp1_Fp2=tmp2_Fp2+f5
	Fp2_mul(&tmp2_Fp2,&A->x0.x2,&B->x1.x0);		//tmp2_Fp2=b3*f2	
	Fp2_add(&tmp3_Fp12.x1.x2,&tmp1_Fp2,&tmp2_Fp2);		//ans[γ^5]=tmp1_Fp2+tmp2_Fp2
	Fp2_mul(&tmp1_Fp2,&A->x1.x0,&B->x1.x0);		//tmp1_Fp2=b3*f3
	Fp2_mul(&tmp2_Fp2,&A->x1.x1,&B->x1.x1);		//tmp2_Fp2=b4*f4
	Fp2_add(&tmp3_Fp2,&A->x1.x0,&A->x1.x1);		//tmp3_Fp2=f3+f4
	Fp2_mul(&tmp3_Fp2,&tmp3_Fp2,&tmp4_Fp2);			//tmp3_Fp2=tmp3_Fp2*tmp4_Fp2

	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp1_Fp2);			//tmp3_Fp2=tmp3_Fp2-tmp1_Fp2
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp2_Fp2);			//tmp3_Fp2=tmp3_Fp2-tmp2_Fp2
	Fp2_add(&tmp3_Fp12.x0.x2,&tmp3_Fp2,&A->x0.x2);	//ans[γ^4]=tmp3_Fp2+f4
	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&A->x0.x1);		//tmp1_Fp2=tmp1_Fp2+f1
	Fp2_mul(&tmp3_Fp2,&A->x1.x2,&B->x1.x1);		//tmp3_Fp2=b4*f5
	Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);			//tmp3_Fp2=tmp3_Fp2*(α+1)
	Fp2_add(&tmp3_Fp12.x0.x1,&tmp1_Fp2,&tmp3_Fp2);		//ans[γ^2]=tmp1_Fp2+tmp3_Fp2
	Fp2_mul(&tmp1_Fp2,&A->x1.x2,&B->x1.x0);		//tmp1_Fp2=b3*f5	
	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&tmp2_Fp2);			//tmp1_Fp2=tmp1_Fp2+tmp2_Fp2
	Fp2_mul_basis(&tmp1_Fp2,&tmp1_Fp2);			//tmp1_Fp2=tmp1_Fp2*(α+1)
	Fp2_add(&tmp3_Fp12.x0.x0,&tmp1_Fp2,&A->x0.x0);	//ans[1]=tmp1_Fp2+f0	
	Fp12_set(ANS,&tmp3_Fp12);	
}
void Pseudo_8_sparse_mul_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2;
    static Fp12 tmp3_Fp12;
	Fp2_mul_lazy(&tmp1_Fp2,&A->x0.x0,&B->x1.x0);
	Fp2_mul_lazy(&tmp2_Fp2,&A->x0.x1,&B->x1.x1);
	Fp2_add_lazy(&tmp3_Fp2,&A->x0.x0,&A->x0.x1);
	Fp2_add_lazy(&tmp4_Fp2,&B->x1.x0,&B->x1.x1);

	Fp2_mul_lazy(&tmp3_Fp2,&tmp3_Fp2,&tmp4_Fp2);
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp1_Fp2);		
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp2_Fp2);
	Fp2_add(&tmp3_Fp12.x1.x1,&tmp3_Fp2,&A->x1.x1);

	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&A->x1.x0);
	Fp2_mul_lazy(&tmp3_Fp2,&A->x0.x2,&B->x1.x1);
	Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);
	Fp2_add(&tmp3_Fp12.x1.x0,&tmp1_Fp2,&tmp3_Fp2);

	Fp2_add(&tmp1_Fp2,&tmp2_Fp2,&A->x1.x2);
	Fp2_mul_lazy(&tmp2_Fp2,&A->x0.x2,&B->x1.x0);
	Fp2_add(&tmp3_Fp12.x1.x2,&tmp1_Fp2,&tmp2_Fp2);

	Fp2_mul_lazy(&tmp1_Fp2,&A->x1.x0,&B->x1.x0);//tmp1
	Fp2_mul_lazy(&tmp2_Fp2,&A->x1.x1,&B->x1.x1);//tmp2

	Fp2_add_lazy(&tmp3_Fp2,&A->x1.x0,&A->x1.x1);
	Fp2_mul_lazy(&tmp3_Fp2,&tmp3_Fp2,&tmp4_Fp2);
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp1_Fp2);
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp2_Fp2);
	Fp2_add(&tmp3_Fp12.x0.x2,&tmp3_Fp2,&A->x0.x2);

	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&A->x0.x1);
	Fp2_mul_lazy(&tmp3_Fp2,&A->x1.x2,&B->x1.x1);
	Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);
	Fp2_add(&tmp3_Fp12.x0.x1,&tmp1_Fp2,&tmp3_Fp2);

	Fp2_mul_lazy(&tmp1_Fp2,&A->x1.x2,&B->x1.x0);
	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&tmp2_Fp2);
	Fp2_mul_basis(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_add(&tmp3_Fp12.x0.x0,&tmp1_Fp2,&A->x0.x0);

	Fp12_set(ANS,&tmp3_Fp12);	
}
void Pseudo_8_sparse_mul_lazy_montgomery(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2;
    static Fp12 tmp3_Fp12;
	Fp2_mul_lazy_montgomery(&tmp1_Fp2,&A->x0.x0,&B->x1.x0);
	Fp2_mul_lazy_montgomery(&tmp2_Fp2,&A->x0.x1,&B->x1.x1);
	Fp2_add_lazy(&tmp3_Fp2,&A->x0.x0,&A->x0.x1);
	Fp2_add_lazy(&tmp4_Fp2,&B->x1.x0,&B->x1.x1);

	Fp2_mul_lazy_montgomery(&tmp3_Fp2,&tmp3_Fp2,&tmp4_Fp2);
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp1_Fp2);		
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp2_Fp2);
	Fp2_add(&tmp3_Fp12.x1.x1,&tmp3_Fp2,&A->x1.x1);

	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&A->x1.x0);
	Fp2_mul_lazy_montgomery(&tmp3_Fp2,&A->x0.x2,&B->x1.x1);
	Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);
	Fp2_add(&tmp3_Fp12.x1.x0,&tmp1_Fp2,&tmp3_Fp2);

	Fp2_add(&tmp1_Fp2,&tmp2_Fp2,&A->x1.x2);
	Fp2_mul_lazy_montgomery(&tmp2_Fp2,&A->x0.x2,&B->x1.x0);
	Fp2_add(&tmp3_Fp12.x1.x2,&tmp1_Fp2,&tmp2_Fp2);

	Fp2_mul_lazy_montgomery(&tmp1_Fp2,&A->x1.x0,&B->x1.x0);//tmp1
	Fp2_mul_lazy_montgomery(&tmp2_Fp2,&A->x1.x1,&B->x1.x1);//tmp2

	Fp2_add_lazy(&tmp3_Fp2,&A->x1.x0,&A->x1.x1);
	Fp2_mul_lazy_montgomery(&tmp3_Fp2,&tmp3_Fp2,&tmp4_Fp2);
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp1_Fp2);
	Fp2_sub(&tmp3_Fp2,&tmp3_Fp2,&tmp2_Fp2);
	Fp2_add(&tmp3_Fp12.x0.x2,&tmp3_Fp2,&A->x0.x2);

	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&A->x0.x1);
	Fp2_mul_lazy_montgomery(&tmp3_Fp2,&A->x1.x2,&B->x1.x1);
	Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);
	Fp2_add(&tmp3_Fp12.x0.x1,&tmp1_Fp2,&tmp3_Fp2);

	Fp2_mul_lazy_montgomery(&tmp1_Fp2,&A->x1.x2,&B->x1.x0);
	Fp2_add(&tmp1_Fp2,&tmp1_Fp2,&tmp2_Fp2);
	Fp2_mul_basis(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_add(&tmp3_Fp12.x0.x0,&tmp1_Fp2,&A->x0.x0);

	Fp12_set(ANS,&tmp3_Fp12);	
}

void ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    static Fp12 tmp1_Fp12,tmp2_Fp12;
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	EFp2_set(&Tmp_T,T);
	
	Fp12_sqr(&tmp1_Fp12,f);
	
	//ltt
	Fp2_add(&tmp1_Fp2,&Tmp_T.y,&Tmp_T.y);		//tmp1_Fp2=1/(2*T.y)
	Fp2_inv(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_sqr(&tmp2_Fp2,&Tmp_T.x);			//tmp2_Fp2=3(T.x)^2
	Fp2_add(&tmp3_Fp2,&tmp2_Fp2,&tmp2_Fp2);
	Fp2_add(&tmp2_Fp2,&tmp3_Fp2,&tmp2_Fp2);
	Fp2_mul(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);				//tmp3_Fp2=tmp1_Fp2*tmp2_Fp2
	
	Fp2_add(&tmp4_Fp2,&Tmp_T.x,&Tmp_T.x);		//tmp4_Fp2=2T.x
	Fp2_sqr(&T->x,&tmp3_Fp2);				//next_T.x=tmp3_Fp2^2-tmp4_Fp2
	Fp2_sub(&T->x,&T->x,&tmp4_Fp2);
	Fp2_mul(&tmp5_Fp2,&tmp3_Fp2,&Tmp_T.x);			//tmp5_Fp2=tmp3_Fp2*T.x-T.y
	Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&Tmp_T.y);
	Fp2_mul(&T->y,&tmp3_Fp2,&T->x);			//next_T.y=tmp5_Fp2-tmp3_Fp2*next_T.x
	Fp2_sub(&T->y,&tmp5_Fp2,&T->y);
	
	//set ltt
	Fp_set_ui(&tmp2_Fp12.x0.x0.x0,1);
	Fp2_set_neg(&tmp2_Fp12.x1.x0,&tmp3_Fp2);
	Fp2_mul_mpn(&tmp2_Fp12.x1.x1,&tmp5_Fp2,L->x0);
	
	Pseudo_8_sparse_mul(f,&tmp1_Fp12,&tmp2_Fp12);
}
void ff_ltt_lazy(Fp12 *f,EFp2 *T,EFp *P,Fp *L){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    static Fp12 tmp1_Fp12,tmp2_Fp12;
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	EFp2_set(&Tmp_T,T);
	
	Fp12_sqr_lazy(&tmp1_Fp12,f);
	
	//ltt
	Fp2_add(&tmp1_Fp2,&Tmp_T.y,&Tmp_T.y);		//tmp1_Fp2=1/(2*T.y)
	Fp2_inv_lazy(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_sqr_lazy(&tmp2_Fp2,&Tmp_T.x);			//tmp2_Fp2=3(T.x)^2
	Fp2_add_lazy(&tmp3_Fp2,&tmp2_Fp2,&tmp2_Fp2);
	Fp2_add_lazy(&tmp2_Fp2,&tmp3_Fp2,&tmp2_Fp2);
	Fp2_mul_lazy(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);		//tmp3_Fp2=tmp1_Fp2*tmp2_Fp2
	
	Fp2_add(&tmp4_Fp2,&Tmp_T.x,&Tmp_T.x);		//tmp4_Fp2=2T.x
	Fp2_sqr_lazy(&T->x,&tmp3_Fp2);				//next_T.x=tmp3_Fp2^2-tmp4_Fp2
	Fp2_sub(&T->x,&T->x,&tmp4_Fp2);
	Fp2_mul_lazy(&tmp5_Fp2,&tmp3_Fp2,&Tmp_T.x);		//tmp5_Fp2=tmp3_Fp2*T.x-T.y
	Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&Tmp_T.y);
	Fp2_mul_lazy(&T->y,&tmp3_Fp2,&T->x);			//next_T.y=tmp5_Fp2-tmp3_Fp2*next_T.x
	Fp2_sub(&T->y,&tmp5_Fp2,&T->y);
	
	//set ltt
	Fp_set_ui(&tmp2_Fp12.x0.x0.x0,1);
	Fp2_set_neg(&tmp2_Fp12.x1.x0,&tmp3_Fp2);
	Fp2_mul_mpn(&tmp2_Fp12.x1.x1,&tmp5_Fp2,L->x0);

	Pseudo_8_sparse_mul_lazy(f,&tmp1_Fp12,&tmp2_Fp12);

}
void ff_ltt_lazy_montgomery(Fp12 *f,EFp2 *T,EFp *P,Fp *L){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    static Fp12 tmp1_Fp12,tmp2_Fp12;
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	EFp2_set(&Tmp_T,T);
	
	Fp12_sqr_lazy(&tmp1_Fp12,f);
	
	//ltt
	Fp2_add(&tmp1_Fp2,&Tmp_T.y,&Tmp_T.y);		//tmp1_Fp2=1/(2*T.y)
	Fp2_inv_lazy_montgomery(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_sqr_lazy_montgomery(&tmp2_Fp2,&Tmp_T.x);			//tmp2_Fp2=3(T.x)^2
	Fp2_add_lazy(&tmp3_Fp2,&tmp2_Fp2,&tmp2_Fp2);
	Fp2_add_lazy(&tmp2_Fp2,&tmp3_Fp2,&tmp2_Fp2);
	Fp2_mul_lazy_montgomery(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);		//tmp3_Fp2=tmp1_Fp2*tmp2_Fp2
	
	Fp2_add(&tmp4_Fp2,&Tmp_T.x,&Tmp_T.x);		//tmp4_Fp2=2T.x
	Fp2_sqr_lazy_montgomery(&T->x,&tmp3_Fp2);				//next_T.x=tmp3_Fp2^2-tmp4_Fp2
	Fp2_sub(&T->x,&T->x,&tmp4_Fp2);
	Fp2_mul_lazy_montgomery(&tmp5_Fp2,&tmp3_Fp2,&Tmp_T.x);		//tmp5_Fp2=tmp3_Fp2*T.x-T.y
	Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&Tmp_T.y);
	Fp2_mul_lazy_montgomery(&T->y,&tmp3_Fp2,&T->x);			//next_T.y=tmp5_Fp2-tmp3_Fp2*next_T.x
	Fp2_sub(&T->y,&tmp5_Fp2,&T->y);
	
	//set ltt
	Fp_set_ui(&tmp2_Fp12.x0.x0.x0,1);
	Fp2_set_neg(&tmp2_Fp12.x1.x0,&tmp3_Fp2);
	Fp2_mul_mpn_montgomery(&tmp2_Fp12.x1.x1,&tmp5_Fp2,L->x0);

	Pseudo_8_sparse_mul_lazy_montgomery(f,&tmp1_Fp12,&tmp2_Fp12);

}
void f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    static Fp12 tmp1_Fp12;
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	EFp2_set(&Tmp_T,T);
		
	//ltq
	Fp2_sub(&tmp1_Fp2,&Q->x,&Tmp_T.x);		//tmp1_Fp2=(Q->x-T.x)^-1
	Fp2_inv(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_sub(&tmp2_Fp2,&Q->y,&Tmp_T.y);		//tmp2_Fp2=(Q->y-T.y)
	Fp2_mul(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);			//tmp3_Fp2=tmp1_Fp2*tmp2_Fp2
	Fp2_add(&tmp4_Fp2,&Tmp_T.x,&Q->x);		//tmp4_Fp2=Q->x+T.x
	Fp2_sqr(&T->x,&tmp3_Fp2);			//next_T.x=tmp3_Fp2^2-tmp4_Fp2
	Fp2_sub(&T->x,&T->x,&tmp4_Fp2);
	Fp2_mul(&tmp5_Fp2,&tmp3_Fp2,&Tmp_T.x);		//tmp5_Fp2=tmp3_Fp2*T.x-T.y
	Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&Tmp_T.y);
	Fp2_mul(&T->y,&tmp3_Fp2,&T->x);		//next_T.y=tmp5_Fp2-tmp3_Fp2*next_T.x
	Fp2_sub(&T->y,&tmp5_Fp2,&T->y);
	
	//set ltq
	Fp_set_ui(&tmp1_Fp12.x0.x0.x0,1);
	Fp2_set_neg(&tmp1_Fp12.x1.x0,&tmp3_Fp2);
	Fp2_mul_mpn(&tmp1_Fp12.x1.x1,&tmp5_Fp2,L->x0);
	
	Pseudo_8_sparse_mul(f,f,&tmp1_Fp12);
}
void f_ltq_lazy(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    static Fp12 tmp1_Fp12;
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	EFp2_set(&Tmp_T,T);
	
	//ltq
	Fp2_sub(&tmp1_Fp2,&Q->x,&Tmp_T.x);		//tmp1_Fp2=(Q->x-T.x)^-1
	Fp2_inv_lazy(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_sub_lazy(&tmp2_Fp2,&Q->y,&Tmp_T.y);		//tmp2_Fp2=(Q->y-T.y)
	Fp2_mul_lazy(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);			//tmp3_Fp2=tmp1_Fp2*tmp2_Fp2
	Fp2_add(&tmp4_Fp2,&Tmp_T.x,&Q->x);		//tmp4_Fp2=Q->x+T.x
	Fp2_sqr_lazy(&T->x,&tmp3_Fp2);			//next_T.x=tmp3_Fp2^2-tmp4_Fp2
	Fp2_sub(&T->x,&T->x,&tmp4_Fp2);
	Fp2_mul_lazy(&tmp5_Fp2,&tmp3_Fp2,&Tmp_T.x);		//tmp5_Fp2=tmp3_Fp2*T.x-T.y
	Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&Tmp_T.y);
	Fp2_mul_lazy(&T->y,&tmp3_Fp2,&T->x);		//next_T.y=tmp5_Fp2-tmp3_Fp2*next_T.x
	Fp2_sub(&T->y,&tmp5_Fp2,&T->y);
	
	//set ltq
	Fp_set_ui(&tmp1_Fp12.x0.x0.x0,1);
	Fp2_set_neg(&tmp1_Fp12.x1.x0,&tmp3_Fp2);
	Fp2_mul_mpn(&tmp1_Fp12.x1.x1,&tmp5_Fp2,L->x0);
	
	Pseudo_8_sparse_mul_lazy(f,f,&tmp1_Fp12);

}
void f_ltq_lazy_montgomery(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    static Fp12 tmp1_Fp12;
    EFp2 Tmp_T;
	EFp2_init(&Tmp_T);
	EFp2_set(&Tmp_T,T);
	
	//ltq
	Fp2_sub(&tmp1_Fp2,&Q->x,&Tmp_T.x);		//tmp1_Fp2=(Q->x-T.x)^-1
	Fp2_inv_lazy_montgomery(&tmp1_Fp2,&tmp1_Fp2);
	Fp2_sub_lazy(&tmp2_Fp2,&Q->y,&Tmp_T.y);		//tmp2_Fp2=(Q->y-T.y)
	Fp2_mul_lazy_montgomery(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);			//tmp3_Fp2=tmp1_Fp2*tmp2_Fp2
	Fp2_add(&tmp4_Fp2,&Tmp_T.x,&Q->x);		//tmp4_Fp2=Q->x+T.x
	Fp2_sqr_lazy_montgomery(&T->x,&tmp3_Fp2);			//next_T.x=tmp3_Fp2^2-tmp4_Fp2
	Fp2_sub(&T->x,&T->x,&tmp4_Fp2);
	Fp2_mul_lazy_montgomery(&tmp5_Fp2,&tmp3_Fp2,&Tmp_T.x);		//tmp5_Fp2=tmp3_Fp2*T.x-T.y
	Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&Tmp_T.y);
	Fp2_mul_lazy_montgomery(&T->y,&tmp3_Fp2,&T->x);		//next_T.y=tmp5_Fp2-tmp3_Fp2*next_T.x
	Fp2_sub(&T->y,&tmp5_Fp2,&T->y);
	
	//set ltq
	Fp_set_ui(&tmp1_Fp12.x0.x0.x0,1);
	Fp2_set_neg(&tmp1_Fp12.x1.x0,&tmp3_Fp2);
	Fp2_mul_mpn_montgomery(&tmp1_Fp12.x1.x1,&tmp5_Fp2,L->x0);
	
	Pseudo_8_sparse_mul_lazy_montgomery(f,f,&tmp1_Fp12);

}
