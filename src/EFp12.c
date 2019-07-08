#include <ELiPS/EFp12.h>
//EFp12
void EFp12_init(EFp12 *P){
    Fp12_init(&P->x);
    Fp12_init(&P->y);
    P->infinity=0;
}
void EFp12_printf(char *str,EFp12 *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp12_printf("",&P->x);
        printf(",");
        Fp12_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}
void EFp12_println(char *str,EFp12 *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp12_printf("",&P->x);
        printf(",");
        Fp12_printf("",&P->y);
        printf(")\n");
    }else{
        printf("0\n");
    }
}
void EFp12_set(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_set_ui(EFp12 *ANS,unsigned long int UI1,unsigned long int UI2){
    Fp12_set_ui(&ANS->x,UI1);
    Fp12_set_ui(&ANS->y,UI2);
    ANS->infinity=0;
}

void EFp12_set_mpn(EFp12 *ANS,mp_limb_t *A){
    Fp12_set_mpn(&ANS->x,A);
    Fp12_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}

void EFp12_set_neg(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

int  EFp12_cmp(EFp12 *A,EFp12 *B){
    if(Fp12_cmp(&A->x,&B->x)==0 && Fp12_cmp(&A->y,&B->y)==0){
        return 0;   
    }else if(A->infinity==1&&B->infinity==1){
	return 0;
    }else{
    return 1;
    }
}

void EFp12_rational_point(EFp12 *P){
    Fp12 tmp1,tmp2;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    gmp_randseed_ui(state,1);
    
    while(1){
        Fp12_set_random(&P->x,state);
        Fp12_sqr(&tmp1,&P->x);
        Fp12_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpn(&tmp2.x0.x0.x0,&tmp2.x0.x0.x0,curve_b);
        if(Fp12_legendre(&tmp2)==1){
            Fp12_sqrt(&P->y,&tmp2);
            break;
        }
    }
}

void BN12_EFp12_generate_G1(EFp12 *P){
    EFp tmp_P;
    EFp_init(&tmp_P);
    
    EFp_rational_point(&tmp_P);
    EFp12_set_ui(P,0,0);
    Fp_set(&P->x.x0.x0.x0,&tmp_P.x);
    Fp_set(&P->y.x0.x0.x0,&tmp_P.y);
    P->infinity=tmp_P.infinity;
}
void BLS12_EFp12_generate_G1(EFp12 *P){
    EFp tmp_P;
    EFp_init(&tmp_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp_rational_point(&tmp_P);
    EFp12_set_ui(P,0,0);
    mpz_tdiv_q(exp,EFp_total,order_z);
    EFp_SCM(&tmp_P,&tmp_P,exp);
    Fp_set(&P->x.x0.x0.x0,&tmp_P.x);
    Fp_set(&P->y.x0.x0.x0,&tmp_P.y);
    P->infinity=tmp_P.infinity;
    
    mpz_clear(exp);
}

void EFp12_generate_G2(EFp12 *Q){
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point(&random_P);
    mpz_pow_ui(exp,order_z,2);
    mpz_tdiv_q(exp,EFp12_total,exp);
    EFp12_SCM(&P,&random_P,exp);
    Fp12_frobenius_map_p1(&frobenius_P.x,&P.x);
    Fp12_frobenius_map_p1(&frobenius_P.y,&P.y);
    EFp12_set_neg(&P,&P);
    EFp12_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
}

void EFp12_ECD(EFp12 *ANS,EFp12 *P){
    static EFp12 tmp1_EFp12;
    static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
    if(Fp12_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp12_set(&tmp1_EFp12,P);
    
    Fp12_add(&tmp1_Fp12,&tmp1_EFp12.y,&tmp1_EFp12.y);
    
    Fp12_inv(&tmp1_Fp12,&tmp1_Fp12);
    Fp12_sqr(&tmp2_Fp12,&tmp1_EFp12.x);
    Fp12_add(&tmp3_Fp12,&tmp2_Fp12,&tmp2_Fp12);
    Fp12_add(&tmp2_Fp12,&tmp2_Fp12,&tmp3_Fp12);
    Fp12_mul(&tmp3_Fp12,&tmp1_Fp12,&tmp2_Fp12);
    
    Fp12_sqr(&tmp1_Fp12,&tmp3_Fp12);
    Fp12_add(&tmp2_Fp12,&tmp1_EFp12.x,&tmp1_EFp12.x);
    Fp12_sub(&ANS->x,&tmp1_Fp12,&tmp2_Fp12);
    
    Fp12_sub(&tmp1_Fp12,&tmp1_EFp12.x,&ANS->x);
    Fp12_mul(&tmp2_Fp12,&tmp3_Fp12,&tmp1_Fp12);
    Fp12_sub(&ANS->y,&tmp2_Fp12,&tmp1_EFp12.y);
}

void EFp12_ECD_lazy(EFp12 *ANS,EFp12 *P){
    static EFp12 tmp1_EFp12;
    static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
    if(Fp12_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp12_set(&tmp1_EFp12,P);
    
    Fp12_add(&tmp1_Fp12,&tmp1_EFp12.y,&tmp1_EFp12.y);
    
    Fp12_inv(&tmp1_Fp12,&tmp1_Fp12);
    Fp12_sqr_lazy(&tmp2_Fp12,&tmp1_EFp12.x);
    Fp12_add_lazy(&tmp3_Fp12,&tmp2_Fp12,&tmp2_Fp12);
    Fp12_add_lazy(&tmp2_Fp12,&tmp2_Fp12,&tmp3_Fp12);
    Fp12_mul_lazy(&tmp3_Fp12,&tmp1_Fp12,&tmp2_Fp12);
    
    Fp12_sqr_lazy(&tmp1_Fp12,&tmp3_Fp12);
    Fp12_add(&tmp2_Fp12,&tmp1_EFp12.x,&tmp1_EFp12.x);
    Fp12_sub(&ANS->x,&tmp1_Fp12,&tmp2_Fp12);
    
    Fp12_sub_lazy(&tmp1_Fp12,&tmp1_EFp12.x,&ANS->x);
    Fp12_mul_lazy(&tmp2_Fp12,&tmp3_Fp12,&tmp1_Fp12);
    Fp12_sub(&ANS->y,&tmp2_Fp12,&tmp1_EFp12.y);
}

void EFp12_ECA(EFp12 *ANS,EFp12 *P1,EFp12 *P2){
    static EFp12 tmp1_EFp12,tmp2_EFp12;
    static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
    if(P1->infinity==1){
        EFp12_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp12_set(ANS,P1);
        return;
    }else if(Fp12_cmp(&P1->x,&P2->x)==0){
        if(Fp12_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp12_ECD(ANS,P1);
            return;
        }
    }
    
    EFp12_set(&tmp1_EFp12,P1);
    EFp12_set(&tmp2_EFp12,P2);
    
    Fp12_sub(&tmp1_Fp12,&tmp2_EFp12.x,&tmp1_EFp12.x);
    Fp12_inv(&tmp1_Fp12,&tmp1_Fp12);
    Fp12_sub(&tmp2_Fp12,&tmp2_EFp12.y,&tmp1_EFp12.y);
    Fp12_mul(&tmp3_Fp12,&tmp1_Fp12,&tmp2_Fp12);
    Fp12_sqr(&tmp1_Fp12,&tmp3_Fp12);
    Fp12_sub(&tmp2_Fp12,&tmp1_Fp12,&tmp1_EFp12.x);
    Fp12_sub(&ANS->x,&tmp2_Fp12,&tmp2_EFp12.x);
    Fp12_sub(&tmp1_Fp12,&tmp1_EFp12.x,&ANS->x);
    Fp12_mul(&tmp2_Fp12,&tmp3_Fp12,&tmp1_Fp12);
    Fp12_sub(&ANS->y,&tmp2_Fp12,&tmp1_EFp12.y);
}
void EFp12_ECA_lazy(EFp12 *ANS,EFp12 *P1,EFp12 *P2){
    static EFp12 tmp1_EFp12,tmp2_EFp12;
    static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
    if(P1->infinity==1){
        EFp12_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp12_set(ANS,P1);
        return;
    }else if(Fp12_cmp(&P1->x,&P2->x)==0){
        if(Fp12_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp12_ECD_lazy(ANS,P1);
            return;
        }
    }
    
    EFp12_set(&tmp1_EFp12,P1);
    EFp12_set(&tmp2_EFp12,P2);
    
    Fp12_sub(&tmp1_Fp12,&tmp2_EFp12.x,&tmp1_EFp12.x);
    Fp12_inv(&tmp1_Fp12,&tmp1_Fp12);
    Fp12_sub_lazy(&tmp2_Fp12,&tmp2_EFp12.y,&tmp1_EFp12.y);
    Fp12_mul_lazy(&tmp3_Fp12,&tmp1_Fp12,&tmp2_Fp12);
    Fp12_sqr_lazy(&tmp1_Fp12,&tmp3_Fp12);
    Fp12_sub_lazy(&tmp2_Fp12,&tmp1_Fp12,&tmp1_EFp12.x);
    Fp12_sub(&ANS->x,&tmp2_Fp12,&tmp2_EFp12.x);
    Fp12_sub_lazy(&tmp1_Fp12,&tmp1_EFp12.x,&ANS->x);
    Fp12_mul_lazy(&tmp2_Fp12,&tmp3_Fp12,&tmp1_Fp12);
    Fp12_sub(&ANS->y,&tmp2_Fp12,&tmp1_EFp12.y);
}

void EFp12_SCM(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp12_set(ANS,P);
        return;
    }
    
    EFp12 Tmp_P,Next_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    EFp12_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFp12_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp12_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp12_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp12_set(ANS,&Next_P);
}
void EFp12_SCM_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp12_set(ANS,P);
        return;
    }
    
    EFp12 Tmp_P,Next_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    EFp12_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFp12_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp12_ECD_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp12_ECA_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp12_set(ANS,&Next_P);
}
