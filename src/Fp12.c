#include <ELiPS/Fp12.h>
//Fp12
void Fp12_init(Fp12 *A){
    Fp6_init(&A->x0);
    Fp6_init(&A->x1);
}

void Fp12_printf(Fp12 *A,char *str){
    gmp_printf("%s(",str);
    Fp6_printf(&A->x0,"");
    gmp_printf(",");
    Fp6_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp12_set(Fp12 *ANS,Fp12 *A){
    Fp6_set(&ANS->x0,&A->x0);
    Fp6_set(&ANS->x1,&A->x1);
}
void Fp12_set_ui(Fp12 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x0,UI);
    Fp6_set_ui(&ANS->x1,0);
}
void Fp12_set_ui_ui(Fp12 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x0,UI);
    Fp6_set_ui(&ANS->x1,UI);
}

void Fp12_set_mpn(Fp12 *ANS,mp_limb_t *A){
    Fp6_set_mpn(&ANS->x0,A);
    Fp6_set_mpn(&ANS->x1,A);    
}

void Fp12_set_neg(Fp12 *ANS,Fp12 *A){
    Fp6_set_neg(&ANS->x0,&A->x0);
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state){
    Fp6_set_random(&ANS->x0,state);
    Fp6_set_random(&ANS->x1,state);
}

void Fp12_mul(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    
    //set
    Fp6_mul(&tmp2_Fp6,&A->x1,&B->x1);//b*d
    Fp6_add(&tmp1_Fp6,&A->x0,&A->x1);//a+b
    Fp6_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp6_mul(&ANS->x1,&tmp1_Fp6,&ANS->x1);//(a+b)(c+d)
    Fp6_mul(&tmp1_Fp6,&A->x0,&B->x0);//a*c
    
    //x0
    Fp6_mul_basis(&ANS->x0,&tmp2_Fp6);//b*d*v
    Fp6_add(&ANS->x0,&ANS->x0,&tmp1_Fp6);//a*c+b*d*v
    
    //x1
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp1_Fp6);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
}
void Fp12_mul_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    //set
    Fp6_mul_lazy(&tmp2_Fp6,&A->x1,&B->x1);//b*d
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);//a+b
    Fp6_add_lazy(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp6_mul_lazy(&ANS->x1,&tmp1_Fp6,&ANS->x1);//(a+b)(c+d)
    Fp6_mul_lazy(&tmp1_Fp6,&A->x0,&B->x0);//a*c
    
    //x0
    Fp6_mul_basis(&ANS->x0,&tmp2_Fp6);//b*d*v
    Fp6_add(&ANS->x0,&ANS->x0,&tmp1_Fp6);//a*c+b*d*v
    
    //x1
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp1_Fp6);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
}
void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_mul_ui(&ANS->x0,&A->x0,UI);
    Fp6_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_mul_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B){
    Fp6_mul_mpn(&ANS->x0,&A->x0,B);
    Fp6_mul_mpn(&ANS->x1,&A->x1,B);
}

void Fp12_sqr(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    Fp6_add(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_mul_basis(&tmp2_Fp6,&A->x1);
    Fp6_add(&tmp2_Fp6,&tmp2_Fp6,&A->x0);
    Fp6_mul(&tmp3_Fp6,&A->x0,&A->x1);
	
    //x0
    Fp6_mul(&ANS->x0,&tmp1_Fp6,&tmp2_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3_Fp6);
    Fp6_mul_basis(&tmp1_Fp6,&tmp3_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp1_Fp6);

    //x1
    Fp6_add(&ANS->x1,&tmp3_Fp6,&tmp3_Fp6);
}

void Fp12_sqr_karat(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    Fp6_mul_basis(&tmp2_Fp6,&A->x1);
    Fp6_add(&tmp1_Fp6,&A->x0,&tmp2_Fp6);
    Fp6_mul(&tmp2_Fp6,&A->x0,&A->x1);
    Fp6_add(&tmp2_Fp6,&tmp2_Fp6,&tmp2_Fp6);
    Fp6_mul_basis(&tmp3_Fp6,&tmp2_Fp6);
	
    //x0
    Fp6_sqr(&ANS->x0,&tmp1_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3_Fp6);

    //x1
    Fp6_set(&ANS->x1,&tmp2_Fp6);
}
void Fp12_sqr_lazy(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_mul_basis_lazy(&tmp2_Fp6,&A->x1);
    Fp6_add_lazy(&tmp2_Fp6,&tmp2_Fp6,&A->x0);
    Fp6_mul_lazy(&tmp3_Fp6,&A->x0,&A->x1);
	
    //x0
    Fp6_mul_lazy(&ANS->x0,&tmp1_Fp6,&tmp2_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3_Fp6);
    Fp6_mul_basis(&tmp1_Fp6,&tmp3_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp1_Fp6);
    
    //x1
    Fp6_add(&ANS->x1,&tmp3_Fp6,&tmp3_Fp6);
}

void Fp12_sqr_cyclotomic(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    //A=a+b*gamma in G3
    //A^2=(1+2b^2*beta)+((a+b)^2-1-b^2*beta-b^2)
    Fp6_add(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_sqr(&tmp1_Fp6,&tmp1_Fp6);           //(a+b)^2
    Fp6_sqr(&tmp2_Fp6,&A->x1);          //b^2
    Fp6_mul_basis(&ANS->x1,&tmp2_Fp6);  //b^2*beta
    Fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    Fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);   //1+2b^2*beta
    
    Fp6_sub(&ANS->x1,&tmp1_Fp6,&ANS->x1);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
    Fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);   //(a+b)^2-1-b^2*beta-b^2
}
void Fp12_sqr_cyclotomic_lazy(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_sqr_lazy(&tmp1_Fp6,&tmp1_Fp6);           //(a+b)^2
    Fp6_sqr_lazy(&tmp2_Fp6,&A->x1);          //b^2
    Fp6_mul_basis(&ANS->x1,&tmp2_Fp6);  //b^2*beta
    Fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    Fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);   //1+2b^2*beta
    
    Fp6_sub(&ANS->x1,&tmp1_Fp6,&ANS->x1);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
    Fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);   //(a+b)^2-1-b^2*beta-b^2
}
void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_add(&ANS->x0,&A->x0,&B->x0);
    Fp6_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp12_add_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_add_lazy(&ANS->x0,&A->x0,&B->x0);
    Fp6_add_lazy(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_add_ui(&ANS->x0,&A->x0,UI);
    Fp6_add_ui(&ANS->x1,&A->x1,0);
}

void Fp12_add_ui_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_add_ui_ui(&ANS->x0,&A->x0,UI);
    Fp6_add_ui_ui(&ANS->x1,&A->x1,UI);
}
void Fp12_add_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B){
    Fp6_add_mpn(&ANS->x0,&ANS->x0,B);
    Fp6_add_mpn(&ANS->x1,&ANS->x1,B);
}

void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_sub(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub(&ANS->x1,&A->x1,&B->x1);
}
void Fp12_sub_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_sub_lazy(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub_lazy(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_sub_ui(&ANS->x0,&ANS->x0,UI);
    Fp6_sub_ui(&ANS->x1,&ANS->x1,0);
}

void Fp12_sub_ui_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_sub_ui_ui(&ANS->x0,&ANS->x0,UI);
    Fp6_sub_ui_ui(&ANS->x1,&ANS->x1,UI);
}
void Fp12_sub_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B){
    Fp6_sub_mpn(&ANS->x0,&ANS->x0,B);
    Fp6_sub_mpn(&ANS->x1,&ANS->x1,B);
}

void Fp12_inv(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6,tmp4_Fp6;
    Fp6_set(&tmp1_Fp6,&A->x0);
    Fp6_set_neg(&tmp2_Fp6,&A->x1);
    
    Fp6_mul(&tmp3_Fp6,&tmp1_Fp6,&A->x0);
    Fp6_mul(&tmp4_Fp6,&tmp2_Fp6,&A->x1);
    Fp6_mul_basis(&tmp4_Fp6,&tmp4_Fp6);
    Fp6_add(&tmp3_Fp6,&tmp3_Fp6,&tmp4_Fp6);
    Fp6_inv(&tmp3_Fp6,&tmp3_Fp6);
    Fp6_mul(&ANS->x0,&tmp1_Fp6,&tmp3_Fp6);
    Fp6_mul(&ANS->x1,&tmp2_Fp6,&tmp3_Fp6);
}
void Fp12_inv_lazy(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6,tmp4_Fp6;
    Fp6_set(&tmp1_Fp6,&A->x0);
    Fp6_set_neg(&tmp2_Fp6,&A->x1);
    
    Fp6_mul_lazy(&tmp3_Fp6,&tmp1_Fp6,&A->x0);
    Fp6_mul_lazy(&tmp4_Fp6,&tmp2_Fp6,&A->x1);
    Fp6_mul_basis(&tmp4_Fp6,&tmp4_Fp6);
    Fp6_add(&tmp3_Fp6,&tmp3_Fp6,&tmp4_Fp6);
    Fp6_inv_lazy(&tmp3_Fp6,&tmp3_Fp6);
    Fp6_mul_lazy(&ANS->x0,&tmp1_Fp6,&tmp3_Fp6);
    Fp6_mul_lazy(&ANS->x1,&tmp2_Fp6,&tmp3_Fp6);
}

int  Fp12_legendre(Fp12 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp12 tmp;
    Fp12_init(&tmp);
    
    mpz_pow_ui(exp,prime_z,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&tmp,A,exp);
    
    if(Fp12_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

void Fp12_sqrt(Fp12 *ANS,Fp12 *A){
    Fp12 x,y,t,k,n,tmp;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp12_set_random(&n,state);
    while(Fp12_legendre(&n)!=-1){
        Fp12_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp12_pow(&y,&n,q);
    mpz_set_ui(z,e);    
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&x,A,exp);
    Fp12_mul(&tmp,&x,&x);
    Fp12_mul(&k,&tmp,A);
    Fp12_mul(&x,&x,A);
    while(Fp12_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp12_pow(&tmp,&k,exp);
        while(Fp12_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp12_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp12_pow(&t,&y,result);
        Fp12_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp12_mul(&x,&x,&t); 
        Fp12_mul(&k,&k,&y);
    }
    Fp12_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
}

void Fp12_pow(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    Fp12 tmp;
    Fp12_init(&tmp);
    Fp12_set(&tmp,A);
    
    for(i=1;i<length; i++){
        Fp12_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp12_mul(&tmp,A,&tmp);
        }
    }
    
    Fp12_set(ANS,&tmp);
}

int  Fp12_cmp(Fp12 *A,Fp12 *B){
    if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI){
    if(Fp6_cmp_ui(&A->x0,UI)==0 && Fp6_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_mpn(Fp12 *A,mp_limb_t *B){
    if(Fp6_cmp_mpn(&A->x0,B)==0 && Fp6_cmp_mpn(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_zero(Fp12 *A){
    if(Fp6_cmp_zero(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_one(Fp12 *A){
    if(Fp6_cmp_one(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

void Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp2_mul_mpn(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant[f_p1][1].x1.x0);
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    Fp2_mul_mpn(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant[f_p1][2].x0.x0);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p1][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p1][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p1][5]);
}
void Fp12_frobenius_map_p1_lazy(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp2_mul_mpn(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant[f_p1][1].x1.x0);
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    Fp2_mul_mpn(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant[f_p1][2].x0.x0);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul_lazy(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p1][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p1][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul_lazy(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p1][5]);
}
void Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p2][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p2][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p2][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p2][5].x0.x0);
}

void Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p3][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p3][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p3][5]);
}
void Fp12_frobenius_map_p3_lazy(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul_lazy(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p3][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p3][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul_lazy(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p3][5]);
}

void Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p4][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p4][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p4][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p4][5].x0.x0);
}

void Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A){
    //x0
    Fp6_set(&ANS->x0,&A->x0);
    //x1
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p8][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p8][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p8][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p8][5].x0.x0);
}

void Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p10][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p10][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p10][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p10][5].x0.x0);
}

