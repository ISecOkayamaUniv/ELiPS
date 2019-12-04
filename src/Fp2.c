#include <ELiPS/fp2.h>
//fp2
void fp2_init(fp2_t *A){
    fp_init(&A->x0);
    fp_init(&A->x1);
}

void fp2_printf(char *str,fp2_t *A){
    gmp_printf("%s(",str);
    fp_printf("",&A->x0);
    gmp_printf(",");
    fp_printf("",&A->x1);
    gmp_printf(")");
}

void fp2_println(char *str,fp2_t *A){
    gmp_printf("%s(",str);
    fp_printf("",&A->x0);
    gmp_printf(",");
    fp_printf("",&A->x1);
    gmp_printf(")\n");
}
void fp2_printf_montgomery(char *str,fp2_t *A){
    gmp_printf("%s(",str);
    fp_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp_printf_montgomery("",&A->x1);
    gmp_printf(")");
}
void fp2_set(fp2_t *ANS,fp2_t *A){
    fp_set(&ANS->x0,&A->x0);
    fp_set(&ANS->x1,&A->x1);
}

void fp2_set_ui(fp2_t *ANS,unsigned long int UI){
    fp_set_ui(&ANS->x0,UI);
    fp_set_ui(&ANS->x1,0);
}
void fp2_set_ui_ui(fp2_t *ANS,unsigned long int UI){
    fp_set_ui(&ANS->x0,UI);
    fp_set_ui(&ANS->x1,UI);
}  
//TODO:set_0
void fp2_set_mpn(fp2_t *ANS,mp_limb_t *A){
    fp_set_mpn(&ANS->x0,A);
    fp_set_ui(&ANS->x1,0);
}

void fp2_set_neg(fp2_t *ANS,fp2_t *A){
    fp_set_neg(&ANS->x0,&A->x0);
    fp_set_neg(&ANS->x1,&A->x1);
}
void fp2_to_montgomery(fp2_t *ANS,fp2_t *A){
    fp_to_montgomery(&ANS->x0,&A->x0);
    fp_to_montgomery(&ANS->x1,&A->x1);
}
void fp2_mod_montgomery(fp2_t *ANS,fp2_t *A){
    fp_mod_montgomery(&ANS->x0,&A->x0);
    fp_mod_montgomery(&ANS->x1,&A->x1);
}
void fp2_lshift(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_lshift(&ANS->x0,&A->x0,UI);
    fp_lshift(&ANS->x1,&A->x1,UI);
}
void fp2_lshift2(fp2_t *ANS,fp2_t *A){
    fp_lshift2(&ANS->x0,&A->x0);
    fp_lshift2(&ANS->x1,&A->x1);
}

void fp2_set_random(fp2_t *ANS,gmp_randstate_t state){
    fp_set_random(&ANS->x0,state);
    fp_set_random(&ANS->x1,state);
}
void fp2_mul(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static fp_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;
	
    //set
    fp_mul(&tmp1_fp,&A->x0,&B->x0);//a*c
    fp_mul(&tmp2_fp,&A->x1,&B->x1);//b*d
    fp_add(&tmp3_fp,&A->x0,&A->x1);//a+b
    fp_add(&tmp4_fp,&B->x0,&B->x1);//c+d
    //x0
    fp_sub(&ANS->x0,&tmp1_fp,&tmp2_fp);//a*c+b*d*v
    //x1
    fp_mul(&ANS->x1,&tmp3_fp,&tmp4_fp);//(a+b)(c+d)
    fp_sub(&ANS->x1,&ANS->x1,&tmp1_fp);
    fp_sub(&ANS->x1,&ANS->x1,&tmp2_fp);
}
void fp2_mul_lazy(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static mp_limb_t buf1L[FPLIMB2],buf2L[FPLIMB2],tmpL1[FPLIMB2],tmpL2[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp3[FPLIMB],tmp4[FPLIMB];

    //set
    fp_mul_lazy(tmpL1,A->x0.x0,B->x0.x0);//a*c
    fp_mul_lazy(tmpL2,A->x1.x0,B->x1.x0);//b*d
	
    fp_add_lazy(tmp3,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);
    fp_add_lazy(tmp4,FPLIMB,B->x0.x0,FPLIMB,B->x1.x0,FPLIMB);

    //x0
    //fp_t out;
    //fp_mod(&out,tmpL1,FPLIMB2);
    //fp_printf("tmp1=",&out);printf("\n");
    //fp_mod(&out,tmpL2,FPLIMB2);
    //fp_printf("tmp2=",&out);printf("\n");
    fp_sub_lazy_mod(&ANS->x0,tmpL1,tmpL2);

    //x1
    fp_mul_lazy(buf1L,tmp3,tmp4);
    fp_sub_lazy(buf2L,FPLIMB2,buf1L,FPLIMB2,tmpL1,FPLIMB2);
    fp_sub_lazy_mod(&ANS->x1,buf2L,tmpL2);
}
void fp2_mul_lazy_montgomery(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static fp_t buf1,buf2,tmp1,tmp2;
    static fp_t buf,tmp3,tmp4;
    
    //set
    fp_mulmod_montgomery(&tmp1,&A->x0,&B->x0);//a*c
    fp_mulmod_montgomery(&tmp2,&A->x1,&B->x1);//b*d
	
    fp_add_lazy(tmp3.x0,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);
    fp_add_lazy(tmp4.x0,FPLIMB,B->x0.x0,FPLIMB,B->x1.x0,FPLIMB);

    //x0
    //fp_printf_montgomery("tmp1=",&tmp1);printf("\n");
    //fp_printf_montgomery("tmp2=",&tmp2);printf("\n");
    fp_sub(&ANS->x0,&tmp1,&tmp2);

    //x1
    fp_mulmod_montgomery(&buf1,&tmp3,&tmp4);
    fp_sub(&buf2,&buf1,&tmp1);
    fp_sub(&ANS->x1,&buf2,&tmp2);
}
void fp2_mul_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_mul_ui(&ANS->x0,&A->x0,UI);
    fp_mul_ui(&ANS->x1,&A->x1,UI);
}
void fp2_mul_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    fp_mul_mpn(&ANS->x0,&A->x0,B);
    fp_mul_mpn(&ANS->x1,&A->x1,B);
}
void fp2_mul_mpn_montgomery(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    mpn_mulmod_montgomery(ANS->x0.x0,FPLIMB,A->x0.x0,FPLIMB,B,FPLIMB);
    mpn_mulmod_montgomery(ANS->x1.x0,FPLIMB,A->x1.x0,FPLIMB,B,FPLIMB);
}
void fp2_mul_basis(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp;
    fp_set(&tmp1_fp,&A->x0);
    
    fp_sub(&ANS->x0,&tmp1_fp,&A->x1);
    fp_add(&ANS->x1,&tmp1_fp,&A->x1);
}
void fp2_mul_basis_lazy(fp2_t *ANS,fp2_t *A){
    static mp_limb_t tmp1[FPLIMB];
    mpn_copyd(tmp1,A->x0.x0,FPLIMB);
    
    fp_sub_lazy(ANS->x0.x0,FPLIMB,tmp1,FPLIMB,A->x1.x0,FPLIMB);
    fp_add_lazy(ANS->x1.x0,FPLIMB,tmp1,FPLIMB,A->x1.x0,FPLIMB);
}
void fp2_add_basis(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static fp2_t tmp1_fp2;
    fp2_mul_basis(&tmp1_fp2,B);
    
    fp2_add(ANS,A,&tmp1_fp2);
}
void fp2_sub_basis(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static fp2_t tmp1_fp2;
    fp2_mul_basis(&tmp1_fp2,B);
    
    fp2_sub(ANS,A,&tmp1_fp2);
}
void fp2_inv_basis(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp,tmp2_fp;
    fp_set(&tmp1_fp,&A->x0);
    fp_set(&tmp2_fp,&A->x1);
    
    fp_add(&ANS->x0,&tmp1_fp,&tmp2_fp);
    fp_mul_mpn(&ANS->x0,&ANS->x0,Alpha_1_inv.x0.x0);
    fp_sub(&ANS->x1,&tmp2_fp,&tmp1_fp);
    fp_mul_mpn(&ANS->x1,&ANS->x1,Alpha_1_inv.x0.x0);
}


//Karat
void fp2_sqr(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp,tmp2_fp;
    fp_add(&tmp1_fp,&A->x0,&A->x1);
    fp_sub(&tmp2_fp,&A->x0,&A->x1);
    //x1
    fp_mul(&ANS->x1,&A->x0,&A->x1);
    fp_add(&ANS->x1,&ANS->x1,&ANS->x1);
    //x0
    fp_mul(&ANS->x0,&tmp1_fp,&tmp2_fp);
}

/*
//complex
void fp2_sqr(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1,tmp2,tmp3;
    
    //Cost = 2M + 4Ap + 3m in fp_t
    fp_mul(&tmp1,&A->x0,&A->x1);//t=a0a1
    fp_add(&tmp2,&A->x0,&A->x1);//(a0+a1)
    fp_set_neg(&tmp3,&A->x1);//a1*i
    fp_add(&tmp3,&tmp3,&A->x0);//(a0+a1*i)
    fp_mul(&ANS->x0,&tmp2,&tmp3);// (a0+a1)(a0+a1*i)
    fp_sub(&ANS->x0, &ANS->x0,&tmp1);
    fp_set_neg(&tmp2,&tmp1);//
    
    fp_sub(&ANS->x0,&ANS->x0,&tmp2);//
    //fp_mul_ui(&ANS->x1,&tmp1,c1);//
    fp_add(&ANS->x1,&tmp1,&tmp1);
}
*/
void fp2_sqr_lazy(fp2_t *ANS,fp2_t *A){
    static mp_limb_t bufL[FPLIMB2],tmpL3[FPLIMB2];
    static mp_limb_t tmp1[FPLIMB],tmp2[FPLIMB];

    fp_add_lazy(tmp1,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);
    fp_sub_lazy(tmp2,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);

    //x1
    fp_mul_lazy(tmpL3,A->x0.x0,A->x1.x0);
    fp_add_lazy(bufL,FPLIMB2,tmpL3,FPLIMB2,tmpL3,FPLIMB2);
    fp_mod(&ANS->x1,bufL,FPLIMB2);

    //x0
    fp_mul_lazy(bufL,tmp1,tmp2);
    fp_mod(&ANS->x0,bufL,FPLIMB2);
}
void fp2_sqr_lazy_montgomery(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1,tmp2,tmp3,buf;

    fp_add_lazy(tmp1.x0,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);
    fp_sub_lazy(tmp2.x0,FPLIMB,A->x0.x0,FPLIMB,A->x1.x0,FPLIMB);

    //x1
    fp_mulmod_montgomery(&tmp3,&A->x0,&A->x1);
    fp_add(&ANS->x1,&tmp3,&tmp3);

    //x0
    fp_mulmod_montgomery(&ANS->x0,&tmp1,&tmp2);
}
void fp2_add(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_add(&ANS->x0,&A->x0,&B->x0);
    fp_add(&ANS->x1,&A->x1,&B->x1);
}

void fp2_add_lazy(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_add_lazy(ANS->x0.x0,FPLIMB,A->x0.x0,FPLIMB,B->x0.x0,FPLIMB);
    fp_add_lazy(ANS->x1.x0,FPLIMB,A->x1.x0,FPLIMB,B->x1.x0,FPLIMB);
}
void fp2_add_final(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_add_final(&ANS->x0,&A->x0,&B->x0);
    fp_add_final(&ANS->x1,&A->x1,&B->x1);
}

void fp2_add_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_add_ui(&ANS->x0,&A->x0,UI);
    fp_add_ui(&ANS->x1,&A->x1,0);
}
void fp2_add_ui_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_add_ui(&ANS->x0,&A->x0,UI);
    fp_add_ui(&ANS->x1,&A->x1,UI);
}
void fp2_add_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    fp_add_mpn(&ANS->x0,&A->x0,B);
    fp_add_mpn(&ANS->x1,&A->x1,B);
}
void fp2_sub(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_sub(&ANS->x0,&A->x0,&B->x0);
    fp_sub(&ANS->x1,&A->x1,&B->x1);
}
void fp2_sub_final(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_sub_final(&ANS->x0,&A->x0,&B->x0);
    fp_sub_final(&ANS->x1,&A->x1,&B->x1);
}
void fp2_sub_lazy(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_sub_lazy(ANS->x0.x0,FPLIMB,A->x0.x0,FPLIMB,B->x0.x0,FPLIMB);
    fp_sub_lazy(ANS->x1.x0,FPLIMB,A->x1.x0,FPLIMB,B->x1.x0,FPLIMB);
}  
void fp2_sub_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_sub_ui(&ANS->x0,&A->x0,UI);
    fp_sub_ui(&ANS->x1,&A->x1,0);
}
void fp2_sub_ui_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_sub_ui(&ANS->x0,&A->x0,UI);
    fp_sub_ui(&ANS->x1,&A->x1,UI);
}
void fp2_sub_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    fp_sub_mpn(&ANS->x0,&A->x0,B);
    fp_sub_mpn(&ANS->x1,&A->x1,B);
}

void fp2_inv(fp2_t *ANS,fp2_t *A){
	static fp_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;
    fp_set(&tmp1_fp,&A->x0);
    fp_set_neg(&tmp2_fp,&A->x1);
    
    fp_sqr(&tmp3_fp,&tmp1_fp);
    fp_mul(&tmp4_fp,&tmp2_fp,&A->x1);
    fp_sub(&tmp3_fp,&tmp3_fp,&tmp4_fp);
    fp_inv(&tmp3_fp,&tmp3_fp);
    fp_mul(&ANS->x0,&tmp1_fp,&tmp3_fp);
    fp_mul(&ANS->x1,&tmp2_fp,&tmp3_fp);
}
void fp2_inv_lazy(fp2_t *ANS,fp2_t *A){
	static fp_t tmp1_fp,tmp2_fp;
	static fp_t tmp3;
	static mp_limb_t tmp1[FPLIMB2],tmp2[FPLIMB2];
    fp_set(&tmp1_fp,&A->x0);
    fp_set_neg(&tmp2_fp,&A->x1);
    
    fp_sqr_lazy(tmp1,tmp1_fp.x0);
    fp_mul_lazy(tmp2,tmp2_fp.x0,A->x1.x0);
    fp_sub_lazy_mod(&tmp3,tmp1,tmp2);
    fp_inv(&tmp3,&tmp3);
    fp_mul(&ANS->x0,&tmp1_fp,&tmp3);
    fp_mul(&ANS->x1,&tmp2_fp,&tmp3);
}
void fp2_inv_lazy_montgomery(fp2_t *ANS,fp2_t *A){
	static fp_t tmp1_fp,tmp2_fp;
	static fp_t tmp3;
	static fp_t tmp1,tmp2;
    fp_set(&tmp1_fp,&A->x0);
    fp_set_neg(&tmp2_fp,&A->x1);
    
    fp_mulmod_montgomery(&tmp1,&tmp1_fp,&tmp1_fp);
    fp_mulmod_montgomery(&tmp2,&tmp2_fp,&A->x1);
    fp_sub(&tmp3,&tmp1,&tmp2);
    fp_inv_montgomery(&tmp3,&tmp3);
    fp_mulmod_montgomery(&ANS->x0,&tmp1_fp,&tmp3);
    fp_mulmod_montgomery(&ANS->x1,&tmp2_fp,&tmp3);
}
int  fp2_legendre(fp2_t *A){
    fp2_t tmp;
    fp2_init(&tmp);
    
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime_z,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp2_pow(&tmp,A,exp);
    
    if(fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

int  fp2_isCNR(fp2_t *A){
    fp2_t tmp;
    fp2_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime_z,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    fp2_pow(&tmp,A,exp);
    
    if(fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }

}

void fp2_sqrt(fp2_t *ANS,fp2_t *A){
    fp2_t x,y,t,k,n,tmp;
    fp2_init(&x);
    fp2_init(&y);
    fp2_init(&t);
    fp2_init(&k);
    fp2_init(&n);
    fp2_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    //gmp_randstate_t state;
	//gmp_randinit_default(state);
	//gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    fp2_set_random(&n,state);
    while(fp2_legendre(&n)!=-1){
        fp2_set_random(&n,state);
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
    fp2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp2_pow(&x,A,exp);
    fp2_mul(&tmp,&x,&x);
    fp2_mul(&k,&tmp,A);
    fp2_mul(&x,&x,A);
    while(fp2_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fp2_pow(&tmp,&k,exp);
        while(fp2_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fp2_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fp2_pow(&t,&y,result);
        fp2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fp2_mul(&x,&x,&t);
        fp2_mul(&k,&k,&y);
    }
    fp2_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}

void fp2_pow(fp2_t *ANS,fp2_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp2_t tmp;
    fp2_init(&tmp);
    fp2_set(&tmp,A);
    
    for(i=1;i<length; i++){
        fp2_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            fp2_mul(&tmp,A,&tmp);
        }
    }
    
    fp2_set(ANS,&tmp);
}

int  fp2_cmp(fp2_t *A,fp2_t *B){
    if(fp_cmp(&A->x0,&B->x0)==0 && fp_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  fp2_cmp_ui(fp2_t *A,unsigned long int UI){
    if(fp_cmp_ui(&A->x0,UI)==0 && fp_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  fp2_cmp_mpn(fp2_t *A,mp_limb_t *B){
    if(fp_cmp_mpn(&A->x0,B)==0 && fp_cmp_mpn(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  fp2_cmp_zero(fp2_t *A){
    if(fp_cmp_zero(&A->x0)==0 && fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  fp2_cmp_one(fp2_t *A){
    if(fp_cmp_one(&A->x0)==0 && fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
int fp2_montgomery_trick(fp2_t *A_inv,fp2_t *A,int n){
    int i;
    fp2_t ANS[n],ALL_inv;
	fp2_set(&ANS[0],&A[0]);
	fp2_t check;
	
	for(i=1;i<n;i++){
	fp2_mul_lazy(&ANS[i],&ANS[i-1],&A[i]);
	}
	fp2_inv_lazy(&ALL_inv,&ANS[n-1]);	
	for(i=n-1;i>0;i--){
    fp2_mul_lazy(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    fp2_mul_lazy(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    fp2_set(&A_inv[0],&ALL_inv);
    /*
    for(i=0;i<n;i++){
    fp2_mul(&check,&A[i],&A_inv[i]);
    printf("check:%d",i);	
	fp2_println("=",&check);
    }
    */
    return 0;
}
int fp2_montgomery_trick_montgomery(fp2_t *A_inv,fp2_t *A,int n){
    int i;
    fp2_t ANS[n],ALL_inv;
	fp2_set(&ANS[0],&A[0]);
	fp2_t check;
	
	for(i=1;i<n;i++){
	fp2_mul_lazy_montgomery(&ANS[i],&ANS[i-1],&A[i]);
	}
	fp2_inv_lazy_montgomery(&ALL_inv,&ANS[n-1]);	
	for(i=n-1;i>0;i--){
    fp2_mul_lazy_montgomery(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    fp2_mul_lazy_montgomery(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    fp2_set(&A_inv[0],&ALL_inv);
    /*
    for(i=0;i<n;i++){
    fp2_mul(&check,&A[i],&A_inv[i]);
    printf("check:%d",i);	
	fp2_println("=",&check);
    }
    */
    return 0;
}