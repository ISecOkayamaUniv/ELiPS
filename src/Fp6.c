#include <ELiPS/Fp6.h>
//Fp6
void Fp6_init(Fp6 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
    Fp2_init(&A->x2);
}

void Fp6_printf(Fp6 *A,char *str){
    gmp_printf("%s(",str);
    Fp2_printf(&A->x0,"");
    gmp_printf(",");
    Fp2_printf(&A->x1,"");
    gmp_printf(",");
    Fp2_printf(&A->x2,"");
    gmp_printf(")");
}

void Fp6_set(Fp6 *ANS,Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
    Fp2_set(&ANS->x2,&A->x2);
}

void Fp6_set_ui(Fp6 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x0,UI);
    Fp2_set_ui(&ANS->x1,0);
    Fp2_set_ui(&ANS->x2,0);
}

void Fp6_set_ui_ui(Fp6 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x0,UI);
    Fp2_set_ui(&ANS->x1,UI);
    Fp2_set_ui(&ANS->x2,UI);
}
void Fp6_set_mpn(Fp6 *ANS,mp_limb_t *A){
    Fp2_set_mpn(&ANS->x0,A);
    Fp2_set_mpn(&ANS->x1,A);
    Fp2_set_mpn(&ANS->x2,A);
}

void Fp6_set_neg(Fp6 *ANS,Fp6 *A){
    Fp2_set_neg(&ANS->x0,&A->x0);
    Fp2_set_neg(&ANS->x1,&A->x1);
    Fp2_set_neg(&ANS->x2,&A->x2);
}

void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state){
    Fp2_set_random(&ANS->x0,state);
    Fp2_set_random(&ANS->x1,state);
    Fp2_set_random(&ANS->x2,state);
}

void Fp6_mul(Fp6 *ANS,Fp6 *A,Fp6 *B){
	static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2,tmp6_Fp2,tmp7_Fp2;
    //set
    Fp2_mul(&tmp1_Fp2,&A->x0,&B->x0);//x0*y0
    Fp2_mul(&tmp2_Fp2,&A->x1,&B->x1);//x1*y1
    Fp2_mul(&tmp3_Fp2,&A->x2,&B->x2);//x2*y2
    
    Fp2_add(&tmp5_Fp2,&A->x0,&A->x1);//x0+x1
    Fp2_add(&tmp4_Fp2,&B->x0,&B->x1);//y0+y1
    Fp2_mul(&tmp5_Fp2,&tmp5_Fp2,&tmp4_Fp2);//(x0+x1)(y0+y1)
    
    Fp2_add(&tmp6_Fp2,&A->x1,&A->x2);//x1+x2
    Fp2_add(&tmp4_Fp2,&B->x1,&B->x2);//y1+y2
    Fp2_mul(&tmp6_Fp2,&tmp6_Fp2,&tmp4_Fp2);//(x1+x2)(y1+y2)
    
    Fp2_add(&tmp7_Fp2,&B->x0,&B->x2);//y2+y0
    Fp2_add(&tmp4_Fp2,&A->x0,&A->x2);//x2+x0
    Fp2_mul(&tmp7_Fp2,&tmp7_Fp2,&tmp4_Fp2);//(x2+x0)(y2+y0)
    //x0
    Fp2_sub(&tmp6_Fp2,&tmp6_Fp2,&tmp2_Fp2);
    Fp2_sub(&tmp6_Fp2,&tmp6_Fp2,&tmp3_Fp2);//(x1+x2)(y1+y2)-x1y1-x2y2
    Fp2_mul_basis(&tmp4_Fp2,&tmp6_Fp2);
    Fp2_add(&ANS->x0,&tmp1_Fp2,&tmp4_Fp2);
    //x1
    Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&tmp1_Fp2);
    Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&tmp2_Fp2);
    Fp2_mul_basis(&tmp4_Fp2,&tmp3_Fp2);
    Fp2_add(&ANS->x1,&tmp4_Fp2,&tmp5_Fp2);
    //x2
    Fp2_sub(&tmp7_Fp2,&tmp7_Fp2,&tmp1_Fp2);
    Fp2_sub(&tmp7_Fp2,&tmp7_Fp2,&tmp3_Fp2);
    Fp2_add(&ANS->x2,&tmp2_Fp2,&tmp7_Fp2);
}
void Fp6_mul_lazy(Fp6 *ANS,Fp6 *A,Fp6 *B){
	static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2,tmp6_Fp2,tmp7_Fp2;
    //set
    Fp2_mul_lazy(&tmp1_Fp2,&A->x0,&B->x0);//tmp1
    Fp2_mul_lazy(&tmp2_Fp2,&A->x1,&B->x1);//tmp2
    Fp2_mul_lazy(&tmp3_Fp2,&A->x2,&B->x2);//tmp3
    
    Fp2_add_lazy(&tmp5_Fp2,&A->x0,&A->x1);
    Fp2_add_lazy(&tmp4_Fp2,&B->x0,&B->x1);
    Fp2_mul_lazy(&tmp5_Fp2,&tmp5_Fp2,&tmp4_Fp2);//tmp5

    Fp2_add_lazy(&tmp6_Fp2,&A->x1,&A->x2);
    Fp2_add_lazy(&tmp4_Fp2,&B->x1,&B->x2);
    Fp2_mul_lazy(&tmp6_Fp2,&tmp6_Fp2,&tmp4_Fp2);//tmp6
    
    Fp2_add_lazy(&tmp7_Fp2,&B->x0,&B->x2);
    Fp2_add_lazy(&tmp4_Fp2,&A->x0,&A->x2);
    Fp2_mul_lazy(&tmp7_Fp2,&tmp7_Fp2,&tmp4_Fp2);//tmp7
    //x0
    Fp2_sub(&tmp6_Fp2,&tmp6_Fp2,&tmp2_Fp2);
    Fp2_sub(&tmp6_Fp2,&tmp6_Fp2,&tmp3_Fp2);
    Fp2_mul_basis(&tmp4_Fp2,&tmp6_Fp2);
    Fp2_add(&ANS->x0,&tmp1_Fp2,&tmp4_Fp2);
    //x1
    Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&tmp1_Fp2);
    Fp2_sub(&tmp5_Fp2,&tmp5_Fp2,&tmp2_Fp2);
    Fp2_mul_basis(&tmp4_Fp2,&tmp3_Fp2);
    Fp2_add(&ANS->x1,&tmp4_Fp2,&tmp5_Fp2);
    //x2
    Fp2_sub(&tmp7_Fp2,&tmp7_Fp2,&tmp1_Fp2);
    Fp2_sub(&tmp7_Fp2,&tmp7_Fp2,&tmp3_Fp2);
    Fp2_add(&ANS->x2,&tmp2_Fp2,&tmp7_Fp2);
}
void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_mul_ui(&ANS->x0,&A->x0,UI);
    Fp2_mul_ui(&ANS->x1,&A->x1,UI);
    Fp2_mul_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_mul_mpn(Fp6 *ANS,Fp6 *A,mp_limb_t *B){
    Fp2_mul_mpn(&ANS->x0,&A->x0,B);
    Fp2_mul_mpn(&ANS->x1,&A->x1,B);
    Fp2_mul_mpn(&ANS->x2,&A->x2,B);
}

void Fp6_mul_basis(Fp6 *ANS,Fp6 *A){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    Fp2_set(&tmp1_Fp2,&A->x0);
    Fp2_set(&tmp2_Fp2,&A->x1);
    Fp2_set(&tmp3_Fp2,&A->x2);
	
    Fp_sub(&ANS->x0.x0,&tmp3_Fp2.x0,&tmp3_Fp2.x1);
    Fp_add(&ANS->x0.x1,&tmp3_Fp2.x0,&tmp3_Fp2.x1);
    Fp_set(&ANS->x1.x0,&tmp1_Fp2.x0);
    Fp_set(&ANS->x1.x1,&tmp1_Fp2.x1);
    Fp_set(&ANS->x2.x0,&tmp2_Fp2.x0);
    Fp_set(&ANS->x2.x1,&tmp2_Fp2.x1);
}

void Fp6_mul_basis_lazy(Fp6 *ANS,Fp6 *A){
    static mp_limb_t tmp1[FPLIMB],tmp2[FPLIMB];
    static Fp2 tmp1_Fp2,tmp2_Fp2;
    Fp2_set(&tmp1_Fp2,&A->x0);
    Fp2_set(&tmp2_Fp2,&A->x1);
    mpn_copyd(tmp1,A->x2.x0.x0,FPLIMB);
    mpn_copyd(tmp2,A->x2.x1.x0,FPLIMB);
	
    Lazy_sub(ANS->x0.x0.x0,FPLIMB,tmp1,FPLIMB,tmp2,FPLIMB);
    Lazy_add(ANS->x0.x1.x0,FPLIMB,tmp1,FPLIMB,tmp2,FPLIMB);

    Fp_set(&ANS->x1.x0,&tmp1_Fp2.x0);
    Fp_set(&ANS->x1.x1,&tmp1_Fp2.x1);
    Fp_set(&ANS->x2.x0,&tmp2_Fp2.x0);
    Fp_set(&ANS->x2.x1,&tmp2_Fp2.x1);
}
void Fp6_sqr(Fp6 *ANS,Fp6 *A){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    Fp2_sqr(&tmp1_Fp2,&A->x0);        //x0^2
    Fp2_sqr(&tmp4_Fp2,&A->x2);        //x2^2
    Fp2_add(&tmp5_Fp2,&A->x1,&A->x1);        //2x1
    Fp2_mul(&tmp2_Fp2,&tmp5_Fp2,&A->x2);  //2x1x2
    Fp2_mul(&tmp3_Fp2,&A->x0,&tmp5_Fp2);  //2x0x1
    Fp2_add(&tmp5_Fp2,&A->x0,&A->x1);        //x0+x1+x2
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&A->x2);
    
    //x0
    Fp2_mul_basis(&ANS->x0,&tmp2_Fp2);
    Fp2_add(&ANS->x0,&ANS->x0,&tmp1_Fp2);
    //x1
    Fp2_mul_basis(&ANS->x1,&tmp4_Fp2);
    Fp2_add(&ANS->x1,&ANS->x1,&tmp3_Fp2);
    //x2
    Fp2_sqr(&ANS->x2,&tmp5_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp1_Fp2,&tmp4_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp2_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp3_Fp2);
    Fp2_sub(&ANS->x2,&ANS->x2,&tmp5_Fp2);
}
void Fp6_sqr_lazy(Fp6 *ANS,Fp6 *A){
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2;
    Fp2_sqr_lazy(&tmp1_Fp2,&A->x0);        //x0^2
    Fp2_sqr_lazy(&tmp4_Fp2,&A->x2);        //x2^2
    Fp2_add_lazy(&tmp5_Fp2,&A->x1,&A->x1);        //2x1

    Fp2_mul_lazy(&tmp2_Fp2,&tmp5_Fp2,&A->x2);  //2x1x2
    Fp2_mul_lazy(&tmp3_Fp2,&A->x0,&tmp5_Fp2);  //2x0x1
    Fp2_add(&tmp5_Fp2,&A->x0,&A->x1);        //x0+x1+x2
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&A->x2);
    
    //x0
    Fp2_mul_basis(&ANS->x0,&tmp2_Fp2);
    Fp2_add(&ANS->x0,&ANS->x0,&tmp1_Fp2);

    //x1
    Fp2_mul_basis(&ANS->x1,&tmp4_Fp2);
    Fp2_add(&ANS->x1,&ANS->x1,&tmp3_Fp2);

    //x2
    Fp2_sqr_lazy(&ANS->x2,&tmp5_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp1_Fp2,&tmp4_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp2_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp3_Fp2);
    Fp2_sub(&ANS->x2,&ANS->x2,&tmp5_Fp2);
}
void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_add(&ANS->x0,&A->x0,&B->x0);
    Fp2_add(&ANS->x1,&A->x1,&B->x1);
    Fp2_add(&ANS->x2,&A->x2,&B->x2);
}
void Fp6_add_final(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_add_final(&ANS->x0,&A->x0,&B->x0);
    Fp2_add_final(&ANS->x1,&A->x1,&B->x1);
    Fp2_add_final(&ANS->x2,&A->x2,&B->x2);
}
void Fp6_add_lazy(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_add_lazy(&ANS->x0,&A->x0,&B->x0);
    Fp2_add_lazy(&ANS->x1,&A->x1,&B->x1);
    Fp2_add_lazy(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_add_ui(&ANS->x0,&A->x0,UI);
    Fp2_add_ui(&ANS->x1,&A->x1,0);
    Fp2_add_ui(&ANS->x2,&A->x2,0);
}

void Fp6_add_ui_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_add_ui_ui(&ANS->x0,&A->x0,UI);
    Fp2_add_ui_ui(&ANS->x1,&A->x1,UI);
    Fp2_add_ui_ui(&ANS->x2,&A->x2,UI);
}
void Fp6_add_mpn(Fp6 *ANS,Fp6 *A,mp_limb_t *B){
    Fp2_add_mpn(&ANS->x0,&A->x0,B);
    Fp2_add_mpn(&ANS->x1,&A->x1,B);
    Fp2_add_mpn(&ANS->x2,&A->x2,B);
}

void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_sub(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub(&ANS->x1,&A->x1,&B->x1);
    Fp2_sub(&ANS->x2,&A->x2,&B->x2);
}
void Fp6_sub_final(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_sub_final(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub_final(&ANS->x1,&A->x1,&B->x1);
    Fp2_sub_final(&ANS->x2,&A->x2,&B->x2);
}
void Fp6_sub_lazy(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_sub_lazy(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub_lazy(&ANS->x1,&A->x1,&B->x1);
    Fp2_sub_lazy(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_sub_ui(&ANS->x0,&A->x0,UI);
    Fp2_sub_ui(&ANS->x1,&A->x1,0);
    Fp2_sub_ui(&ANS->x2,&A->x2,0);
}

void Fp6_sub_ui_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_sub_ui_ui(&ANS->x0,&A->x0,UI);
    Fp2_sub_ui_ui(&ANS->x1,&A->x1,UI);
    Fp2_sub_ui_ui(&ANS->x2,&A->x2,UI);
}
void Fp6_sub_mpn(Fp6 *ANS,Fp6 *A,mp_limb_t *B){
    Fp2_sub_mpn(&ANS->x0,&A->x0,B);
    Fp2_sub_mpn(&ANS->x1,&A->x1,B);
    Fp2_sub_mpn(&ANS->x2,&A->x2,B);
}

void Fp6_inv(Fp6 *ANS,Fp6 *A){
	static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2,tmp6_Fp2,tmp7_Fp2,tmp8_Fp2;

    Fp2_sqr(&tmp1_Fp2,&A->x0);
    Fp2_sqr(&tmp2_Fp2,&A->x1);
    Fp2_sqr(&tmp3_Fp2,&A->x2);
    
    Fp2_mul(&tmp4_Fp2,&A->x1,&A->x2);
    Fp2_mul_basis(&tmp4_Fp2,&tmp4_Fp2);
    Fp2_sub(&tmp6_Fp2,&tmp1_Fp2,&tmp4_Fp2);
    
    Fp2_mul(&tmp4_Fp2,&A->x0,&A->x1);
    Fp2_mul_basis(&tmp7_Fp2,&tmp3_Fp2);
    Fp2_sub(&tmp7_Fp2,&tmp7_Fp2,&tmp4_Fp2);
    
    Fp2_mul(&tmp4_Fp2,&A->x0,&A->x2);
    Fp2_sub(&tmp8_Fp2,&tmp2_Fp2,&tmp4_Fp2);
    
    Fp2_mul(&tmp1_Fp2,&tmp1_Fp2,&A->x0);
    Fp2_mul(&tmp3_Fp2,&tmp3_Fp2,&A->x2);
    Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);
    
    Fp2_add(&tmp5_Fp2,&tmp4_Fp2,&tmp4_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp4_Fp2);
    Fp2_sub(&tmp5_Fp2,&tmp2_Fp2,&tmp5_Fp2);
    Fp2_mul(&tmp5_Fp2,&tmp5_Fp2,&A->x1);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp3_Fp2);
    Fp2_mul_basis(&tmp5_Fp2,&tmp5_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp1_Fp2);
    
    Fp2_inv(&tmp5_Fp2,&tmp5_Fp2);
    
    Fp2_mul(&ANS->x0,&tmp6_Fp2,&tmp5_Fp2);
    Fp2_mul(&ANS->x1,&tmp7_Fp2,&tmp5_Fp2);
    Fp2_mul(&ANS->x2,&tmp8_Fp2,&tmp5_Fp2);
}
void Fp6_inv_lazy(Fp6 *ANS,Fp6 *A){
	static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2,tmp4_Fp2,tmp5_Fp2,tmp6_Fp2,tmp7_Fp2,tmp8_Fp2;

    Fp2_sqr_lazy(&tmp1_Fp2,&A->x0);
    Fp2_sqr_lazy(&tmp2_Fp2,&A->x1);
    Fp2_sqr_lazy(&tmp3_Fp2,&A->x2);
    
    Fp2_mul_lazy(&tmp4_Fp2,&A->x1,&A->x2);
    Fp2_mul_basis(&tmp4_Fp2,&tmp4_Fp2);
    Fp2_sub(&tmp6_Fp2,&tmp1_Fp2,&tmp4_Fp2);//tmp6
    
    Fp2_mul_lazy(&tmp4_Fp2,&A->x0,&A->x1);
    Fp2_mul_basis(&tmp7_Fp2,&tmp3_Fp2);
    Fp2_sub(&tmp7_Fp2,&tmp7_Fp2,&tmp4_Fp2);//tmp7
    
    Fp2_mul_lazy(&tmp4_Fp2,&A->x0,&A->x2);//tmp4
    Fp2_sub(&tmp8_Fp2,&tmp2_Fp2,&tmp4_Fp2);//tmp8
    
    Fp2_mul_lazy(&tmp1_Fp2,&tmp1_Fp2,&A->x0);//tmp1
    Fp2_mul_lazy(&tmp3_Fp2,&tmp3_Fp2,&A->x2);
    Fp2_mul_basis(&tmp3_Fp2,&tmp3_Fp2);//tmp3
    
    Fp2_add_lazy(&tmp5_Fp2,&tmp4_Fp2,&tmp4_Fp2);
    Fp2_add_lazy(&tmp5_Fp2,&tmp5_Fp2,&tmp4_Fp2);
    Fp2_sub_lazy(&tmp5_Fp2,&tmp2_Fp2,&tmp5_Fp2);
    Fp2_mul_lazy(&tmp5_Fp2,&tmp5_Fp2,&A->x1);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp3_Fp2);
    Fp2_mul_basis(&tmp5_Fp2,&tmp5_Fp2);
    Fp2_add(&tmp5_Fp2,&tmp5_Fp2,&tmp1_Fp2);//mod
    
    Fp2_inv_lazy(&tmp5_Fp2,&tmp5_Fp2);
    
    Fp2_mul_lazy(&ANS->x0,&tmp6_Fp2,&tmp5_Fp2);
    Fp2_mul_lazy(&ANS->x1,&tmp7_Fp2,&tmp5_Fp2);
    Fp2_mul_lazy(&ANS->x2,&tmp8_Fp2,&tmp5_Fp2);
}
int  Fp6_legendre(Fp6 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp6 tmp;
    Fp6_init(&tmp);
    
    mpz_pow_ui(exp,prime_z,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&tmp,A,exp);
    
    if(Fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

int  Fp6_isCNR(Fp6 *A){
    Fp6 tmp;
    Fp6_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime_z,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp6_pow(&tmp,A,exp);
    
    if(Fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

void Fp6_sqrt(Fp6 *ANS,Fp6 *A){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp2_set(&tmp1.x0,&A->x0);
    Fp2_mul_mpn(&tmp1.x1,&A->x1,frobenius_constant[f_p4][1].x0.x0);
    Fp2_mul_mpn(&tmp1.x2,&A->x2,frobenius_constant[f_p4][2].x0.x0);
    
    Fp2_set(&tmp2.x0,&A->x0);
    Fp2_mul_mpn(&tmp2.x1,&A->x1,frobenius_constant[f_p2][1].x0.x0);
    Fp2_mul_mpn(&tmp2.x2,&A->x2,frobenius_constant[f_p2][2].x0.x0);
    
    Fp6_mul(&tmp1,&tmp1,&tmp2);
    Fp6_mul(&tmp1,&tmp1,A);
    Fp6_set_ui(&tmp2,0);
    Fp2_sqrt(&tmp2.x0,&tmp1.x0);
    Fp2_inv(&tmp2.x0,&tmp2.x0);
    Fp2_set(&tmp2.x0,&tmp2.x0);
    mpz_pow_ui(exp,prime_z,8);
    mpz_pow_ui(buf,prime_z,4);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&tmp1,A,exp);
    Fp6_mul(&tmp1,&tmp1,&tmp2);
    Fp6_set(ANS,&tmp1);
    
    mpz_clear(exp);
    mpz_clear(buf);
}

void Fp6_pow(Fp6 *ANS,Fp6 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    
    for(i=1;i<length; i++){
        Fp6_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp6_mul(&tmp,A,&tmp);
        }
    }
    
    Fp6_set(ANS,&tmp);
}

int  Fp6_cmp(Fp6 *A,Fp6 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0 && Fp2_cmp(&A->x2,&B->x2)==0){
        return 0;   
    }
    return 1;
}

int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI){
    if(Fp2_cmp_ui(&A->x0,UI)==0 && Fp2_cmp_ui(&A->x1,UI)==0 && Fp2_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_mpn(Fp6 *A,mp_limb_t *B){
    if(Fp2_cmp_mpn(&A->x0,B)==0 && Fp2_cmp_mpn(&A->x1,B)==0 && Fp2_cmp_mpn(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_zero(Fp6 *A){
    if(Fp2_cmp_zero(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_one(Fp6 *A){
    if(Fp2_cmp_one(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}
