#include <ELiPS/EFp6.h>
//EFp6
void EFp6_init(EFp6 *P){
    Fp6_init(&P->x);
    Fp6_init(&P->y);
    P->infinity=0;
}

void EFp6_printf(char *str,EFp6 *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp6_printf("",&P->x);
        printf(",");
        Fp6_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}
void EFp6_println(char *str,EFp6 *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp6_printf("",&P->x);
        printf(",");
        Fp6_printf("",&P->y);
        printf(")\n");
    }else{
        printf("0\n");
    }
}
void EFp6_set(EFp6 *ANS,EFp6 *A){
    Fp6_set(&ANS->x,&A->x);
    Fp6_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp6_set_ui(EFp6 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x,UI);
    Fp6_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp6_set_mpn(EFp6 *ANS,mp_limb_t *A){
    Fp6_set_mpn(&ANS->x,A);
    Fp6_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}

void EFp6_set_neg(EFp6 *ANS,EFp6 *A){
    Fp6_set(&ANS->x,&A->x);
    Fp6_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}


void EFp6_rational_point(EFp6 *P){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp6_set_random(&P->x,state);
        Fp6_sqr(&tmp1,&P->x);
        Fp6_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpn(&tmp2.x0.x0,&tmp2.x0.x0,curve_b);
        if(Fp6_legendre(&tmp2)==1){
            Fp6_sqrt(&P->y,&tmp2);
            break;
        }
    }
}

void EFp6_ECD(EFp6 *ANS,EFp6 *P){
	static EFp6 tmp1_EFp6;
	static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    if(Fp6_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp6_set(&tmp1_EFp6,P);
    
    Fp6_add(&tmp1_Fp6,&tmp1_EFp6.y,&tmp1_EFp6.y);
    
    Fp6_inv(&tmp1_Fp6,&tmp1_Fp6);
    Fp6_sqr(&tmp2_Fp6,&tmp1_EFp6.x);
    Fp6_add(&tmp3_Fp6,&tmp2_Fp6,&tmp2_Fp6);
    Fp6_add(&tmp2_Fp6,&tmp2_Fp6,&tmp3_Fp6);
    Fp6_mul(&tmp3_Fp6,&tmp1_Fp6,&tmp2_Fp6);
    
    Fp6_sqr(&tmp1_Fp6,&tmp3_Fp6);
    Fp6_add(&tmp2_Fp6,&tmp1_EFp6.x,&tmp1_EFp6.x);
    Fp6_sub(&ANS->x,&tmp1_Fp6,&tmp2_Fp6);
    
    Fp6_sub(&tmp1_Fp6,&tmp1_EFp6.x,&ANS->x);
    Fp6_mul(&tmp2_Fp6,&tmp3_Fp6,&tmp1_Fp6);
    Fp6_sub(&ANS->y,&tmp2_Fp6,&tmp1_EFp6.y);
}

void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2){
	static EFp6 tmp1_EFp6,tmp2_EFp6;
	static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    if(P1->infinity==1){
        EFp6_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp6_set(ANS,P1);
        return;
    }else if(Fp6_cmp(&P1->x,&P2->x)==0){
        if(Fp6_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp6_ECD(ANS,P1);
            return;
        }
    }
    
    EFp6_set(&tmp1_EFp6,P1);
    EFp6_set(&tmp2_EFp6,P2);
    
    Fp6_sub(&tmp1_Fp6,&tmp2_EFp6.x,&tmp1_EFp6.x);
    Fp6_inv(&tmp1_Fp6,&tmp1_Fp6);
    Fp6_sub(&tmp2_Fp6,&tmp2_EFp6.y,&tmp1_EFp6.y);
    Fp6_mul(&tmp3_Fp6,&tmp1_Fp6,&tmp2_Fp6);
    Fp6_sqr(&tmp1_Fp6,&tmp3_Fp6);
    Fp6_sub(&tmp2_Fp6,&tmp1_Fp6,&tmp1_EFp6.x);
    Fp6_sub(&ANS->x,&tmp2_Fp6,&tmp2_EFp6.x);
    Fp6_sub(&tmp1_Fp6,&tmp1_EFp6.x,&ANS->x);
    Fp6_mul(&tmp2_Fp6,&tmp3_Fp6,&tmp1_Fp6);
    Fp6_sub(&ANS->y,&tmp2_Fp6,&tmp1_EFp6.y);
}

void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp6_set(ANS,P);
        return;
    }
    
    EFp6 Tmp_P,Next_P;
    EFp6_init(&Tmp_P);
    EFp6_set(&Tmp_P,P);
    EFp6_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFp6_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp6_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp6_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp6_set(ANS,&Next_P);
}
