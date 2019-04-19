#include <ELiPS/Fp2.h>
//Fp2
void Fp2_init(Fp2 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
}

void Fp2_printf(char *str,Fp2 *A){
    gmp_printf("%s(",str);
    Fp_printf("",&A->x0);
    gmp_printf(",");
    Fp_printf("",&A->x1);
    gmp_printf(")");
}

void Fp2_set(Fp2 *ANS,Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
}

void Fp2_set_ui(Fp2 *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x0,UI);
    Fp_set_ui(&ANS->x1,0);
}
void Fp2_set_ui_ui(Fp2 *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x0,UI);
    Fp_set_ui(&ANS->x1,UI);
}  

void Fp2_set_mpn(Fp2 *ANS,mp_limb_t *A){
    Fp_set_mpn(&ANS->x0,A);
    Fp_set_mpn(&ANS->x1,A);
}

void Fp2_set_neg(Fp2 *ANS,Fp2 *A){
    Fp_set_neg(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
}

void Fp2_lshift(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_lshift(&ANS->x0,&A->x0,UI);
    Fp_lshift(&ANS->x1,&A->x1,UI);
}

void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state){
    Fp_set_random(&ANS->x0,state);
    Fp_set_random(&ANS->x1,state);
}
void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B){
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp,tmp4_Fp;
	
    //set
    Fp_mul(&tmp1_Fp,&A->x0,&B->x0);//a*c
    Fp_mul(&tmp2_Fp,&A->x1,&B->x1);//b*d
    Fp_add(&tmp3_Fp,&A->x0,&A->x1);//a+b
    Fp_add(&tmp4_Fp,&B->x0,&B->x1);//c+d
    //x0
    Fp_sub(&ANS->x0,&tmp1_Fp,&tmp2_Fp);//a*c+b*d*v
    //x1
    Fp_mul(&ANS->x1,&tmp3_Fp,&tmp4_Fp);//(a+b)(c+d)
    Fp_sub(&ANS->x1,&ANS->x1,&tmp1_Fp);
    Fp_sub(&ANS->x1,&ANS->x1,&tmp2_Fp);
}
void Fp2_mul_lazy(Fp2 *ANS,Fp2 *A,Fp2 *B){
    static mp_limb_t buf1L[FPLIMB2],buf2L[FPLIMB2],tmpL1[FPLIMB2],tmpL2[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp3[FPLIMB],tmp4[FPLIMB];

    //set
    Lazy_mul(tmpL1,A->x0.x0,B->x0.x0);//a*c
    Lazy_mul(tmpL2,A->x1.x0,B->x1.x0);//b*d
	
    Lazy_add(tmp3,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);
    Lazy_add(tmp4,FPLIMB,B->x0.x0,FPLIMB,B->x1.x0,FPLIMB);

    //x0
    Lazy_sub_mod(&ANS->x0,tmpL1,tmpL2);

    //x1
    Lazy_mul(buf1L,tmp3,tmp4);
    Lazy_sub(buf2L,FPLIMB2,buf1L,FPLIMB2,tmpL1,FPLIMB2);
    Lazy_sub_mod(&ANS->x1,buf2L,tmpL2);
}
void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_mul_ui(&ANS->x0,&A->x0,UI);
    Fp_mul_ui(&ANS->x1,&A->x1,UI);
}
void Fp2_mul_mpn(Fp2 *ANS,Fp2 *A,mp_limb_t *B){
    Fp_mul_mpn(&ANS->x0,&A->x0,B);
    Fp_mul_mpn(&ANS->x1,&A->x1,B);
}
void Fp2_mul_basis(Fp2 *ANS,Fp2 *A){
    static Fp tmp1_Fp;
    Fp_set(&tmp1_Fp,&A->x0);
    
    Fp_sub(&ANS->x0,&tmp1_Fp,&A->x1);
    Fp_add(&ANS->x1,&tmp1_Fp,&A->x1);
}
void Fp2_mul_basis_lazy(Fp2 *ANS,Fp2 *A){
    static mp_limb_t tmp1[FPLIMB];
    mpn_copyd(tmp1,A->x0.x0,FPLIMB);
    
    Lazy_sub(ANS->x0.x0,FPLIMB,tmp1,FPLIMB,A->x1.x0,FPLIMB);
    Lazy_add(ANS->x1.x0,FPLIMB,tmp1,FPLIMB,A->x1.x0,FPLIMB);
}
void Fp2_add_basis(Fp2 *ANS,Fp2 *A,Fp2 *B){
    static Fp2 tmp1_Fp2;
    Fp2_mul_basis(&tmp1_Fp2,B);
    
    Fp2_add(ANS,A,&tmp1_Fp2);
}
void Fp2_sub_basis(Fp2 *ANS,Fp2 *A,Fp2 *B){
    static Fp2 tmp1_Fp2;
    Fp2_mul_basis(&tmp1_Fp2,B);
    
    Fp2_sub(ANS,A,&tmp1_Fp2);
}
void Fp2_inv_basis(Fp2 *ANS,Fp2 *A){
    static Fp tmp1_Fp,tmp2_Fp;
    Fp_set(&tmp1_Fp,&A->x0);
    Fp_set(&tmp2_Fp,&A->x1);
    
    Fp_add(&ANS->x0,&tmp1_Fp,&tmp2_Fp);
    Fp_mul_mpn(&ANS->x0,&ANS->x0,Alpha_1_inv.x0.x0);
    Fp_sub(&ANS->x1,&tmp2_Fp,&tmp1_Fp);
    Fp_mul_mpn(&ANS->x1,&ANS->x1,Alpha_1_inv.x0.x0);
}


//Karat
void Fp2_sqr(Fp2 *ANS,Fp2 *A){
    static Fp tmp1_Fp,tmp2_Fp;
    Fp_add(&tmp1_Fp,&A->x0,&A->x1);
    Fp_sub(&tmp2_Fp,&A->x0,&A->x1);
    //x1
    Fp_mul(&ANS->x1,&A->x0,&A->x1);
    Fp_add(&ANS->x1,&ANS->x1,&ANS->x1);
    //x0
    Fp_mul(&ANS->x0,&tmp1_Fp,&tmp2_Fp);
}

/*
//complex
void Fp2_sqr(Fp2 *ANS,Fp2 *A){
    static Fp tmp1,tmp2,tmp3;
    
    //Cost = 2M + 4Ap + 3m in Fp
    Fp_mul(&tmp1,&A->x0,&A->x1);//t=a0a1
    Fp_add(&tmp2,&A->x0,&A->x1);//(a0+a1)
    Fp_set_neg(&tmp3,&A->x1);//a1*i
    Fp_add(&tmp3,&tmp3,&A->x0);//(a0+a1*i)
    Fp_mul(&ANS->x0,&tmp2,&tmp3);// (a0+a1)(a0+a1*i)
    Fp_sub(&ANS->x0, &ANS->x0,&tmp1);
    Fp_set_neg(&tmp2,&tmp1);//
    
    Fp_sub(&ANS->x0,&ANS->x0,&tmp2);//
    //Fp_mul_ui(&ANS->x1,&tmp1,c1);//
    Fp_add(&ANS->x1,&tmp1,&tmp1);
}
*/
void Fp2_sqr_lazy(Fp2 *ANS,Fp2 *A){
    static mp_limb_t bufL[FPLIMB2],tmpL3[FPLIMB2];
    static mp_limb_t tmp1[FPLIMB],tmp2[FPLIMB];

    Lazy_add(tmp1,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);
    Lazy_sub(tmp2,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);

    //x1
    Lazy_mul(tmpL3,A->x0.x0,A->x1.x0);
    Lazy_add(bufL,FPLIMB2,tmpL3,FPLIMB2,tmpL3,FPLIMB2);
    mpn_mod(&ANS->x1,bufL,FPLIMB2);

    //x0
    Lazy_mul(bufL,tmp1,tmp2);
    mpn_mod(&ANS->x0,bufL,FPLIMB2);
}
void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_add_lazy(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Lazy_add(ANS->x0.x0,FPLIMB,A->x0.x0,FPLIMB,B->x0.x0,FPLIMB);
    Lazy_add(ANS->x1.x0,FPLIMB,A->x1.x0,FPLIMB,B->x1.x0,FPLIMB);
}
void Fp2_add_final(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_add_final(&ANS->x0,&A->x0,&B->x0);
    Fp_add_final(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_add_ui(&ANS->x0,&A->x0,UI);
    Fp_add_ui(&ANS->x1,&A->x1,0);
}
void Fp2_add_ui_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_add_ui(&ANS->x0,&A->x0,UI);
    Fp_add_ui(&ANS->x1,&A->x1,UI);
}
void Fp2_add_mpn(Fp2 *ANS,Fp2 *A,mp_limb_t *B){
    Fp_add_mpn(&ANS->x0,&A->x0,B);
    Fp_add_mpn(&ANS->x1,&A->x1,B);
}
void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
}
void Fp2_sub_final(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_sub_final(&ANS->x0,&A->x0,&B->x0);
    Fp_sub_final(&ANS->x1,&A->x1,&B->x1);
}
void Fp2_sub_lazy(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Lazy_sub(ANS->x0.x0,FPLIMB,A->x0.x0,FPLIMB,B->x0.x0,FPLIMB);
    Lazy_sub(ANS->x1.x0,FPLIMB,A->x1.x0,FPLIMB,B->x1.x0,FPLIMB);
}  
void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_sub_ui(&ANS->x0,&A->x0,UI);
    Fp_sub_ui(&ANS->x1,&A->x1,0);
}
void Fp2_sub_ui_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_sub_ui(&ANS->x0,&A->x0,UI);
    Fp_sub_ui(&ANS->x1,&A->x1,UI);
}
void Fp2_sub_mpn(Fp2 *ANS,Fp2 *A,mp_limb_t *B){
    Fp_sub_mpn(&ANS->x0,&A->x0,B);
    Fp_sub_mpn(&ANS->x1,&A->x1,B);
}

void Fp2_inv(Fp2 *ANS,Fp2 *A){
	static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp,tmp4_Fp;
    Fp_set(&tmp1_Fp,&A->x0);
    Fp_set_neg(&tmp2_Fp,&A->x1);
    
    Fp_mul(&tmp3_Fp,&tmp1_Fp,&A->x0);
    Fp_mul(&tmp4_Fp,&tmp2_Fp,&A->x1);
    Fp_sub(&tmp3_Fp,&tmp3_Fp,&tmp4_Fp);
    Fp_inv(&tmp3_Fp,&tmp3_Fp);
    Fp_mul(&ANS->x0,&tmp1_Fp,&tmp3_Fp);
    Fp_mul(&ANS->x1,&tmp2_Fp,&tmp3_Fp);
}
void Fp2_inv_lazy(Fp2 *ANS,Fp2 *A){
	static Fp tmp1_Fp,tmp2_Fp;
	static Fp tmp3;
	static mp_limb_t tmp1[FPLIMB2],tmp2[FPLIMB2];
    Fp_set(&tmp1_Fp,&A->x0);
    Fp_set_neg(&tmp2_Fp,&A->x1);
    
    Lazy_mul(tmp1,tmp1_Fp.x0,A->x0.x0);
    Lazy_mul(tmp2,tmp2_Fp.x0,A->x1.x0);
    Lazy_sub_mod(&tmp3,tmp1,tmp2);
    Fp_inv(&tmp3,&tmp3);
    Fp_mul(&ANS->x0,&tmp1_Fp,&tmp3);
    Fp_mul(&ANS->x1,&tmp2_Fp,&tmp3);
}
int  Fp2_legendre(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime_z,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

int  Fp2_isCNR(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime_z,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }

}

void Fp2_sqrt(Fp2 *ANS,Fp2 *A){
    Fp2 x,y,t,k,n,tmp;
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&t);
    Fp2_init(&k);
    Fp2_init(&n);
    Fp2_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp2_set_random(&n,state);
    while(Fp2_legendre(&n)!=-1){
        Fp2_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,2);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&x,A,exp);
    Fp2_mul(&tmp,&x,&x);
    Fp2_mul(&k,&tmp,A);
    Fp2_mul(&x,&x,A);
    while(Fp2_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp2_pow(&tmp,&k,exp);
        while(Fp2_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp2_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp2_pow(&t,&y,result);
        Fp2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&k,&k,&y);
    }
    Fp2_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}

void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp,A);
    
    for(i=1;i<length; i++){
        Fp2_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp2_mul(&tmp,A,&tmp);
        }
    }
    
    Fp2_set(ANS,&tmp);
}

int  Fp2_cmp(Fp2 *A,Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI){
    if(Fp_cmp_ui(&A->x0,UI)==0 && Fp_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_mpn(Fp2 *A,mp_limb_t *B){
    if(Fp_cmp_mpn(&A->x0,B)==0 && Fp_cmp_mpn(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_zero(Fp2 *A){
    if(Fp_cmp_zero(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_one(Fp2 *A){
    if(Fp_cmp_one(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
