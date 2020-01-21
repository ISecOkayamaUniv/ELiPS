#include <ELiPS/efp6.h>
//efp6
void efp6_init(efp6_t *P){
    fp6_init(&P->x);
    fp6_init(&P->y);
    P->infinity=0;
}

void efp6_printf(char *str,efp6_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp6_printf("",&P->x);
        printf(",");
        fp6_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}
void efp6_println(char *str,efp6_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp6_printf("",&P->x);
        printf(",");
        fp6_printf("",&P->y);
        printf(")\n");
    }else{
        printf("0\n");
    }
}
void efp6_set(efp6_t *ANS,efp6_t *A){
    fp6_set(&ANS->x,&A->x);
    fp6_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void efp6_set_ui(efp6_t *ANS,unsigned long int UI){
    fp6_set_ui(&ANS->x,UI);
    fp6_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void efp6_set_mpn(efp6_t *ANS,mp_limb_t *A){
    fp6_set_mpn(&ANS->x,A);
    fp6_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}

void efp6_set_neg(efp6_t *ANS,efp6_t *A){
    fp6_set(&ANS->x,&A->x);
    fp6_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}


void efp6_rational_point(efp6_t *P){
    fp6_t tmp1,tmp2;
    fp6_init(&tmp1);
    fp6_init(&tmp2);
	gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        fp6_set_random(&P->x,state);
        fp6_sqr(&tmp1,&P->x);
        fp6_mul(&tmp2,&tmp1,&P->x);
        fp_add_mpn(&tmp2.x0.x0,&tmp2.x0.x0,curve_b);
        if(fp6_legendre(&tmp2)==1){
            fp6_sqrt(&P->y,&tmp2);
            break;
        }
    }
}

void efp6_ecd(efp6_t *ANS,efp6_t *P){
	static efp6_t tmp1_efp6;
	static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
    if(fp6_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    efp6_set(&tmp1_efp6,P);
    
    fp6_add(&tmp1_fp6,&tmp1_efp6.y,&tmp1_efp6.y);
    
    fp6_inv(&tmp1_fp6,&tmp1_fp6);
    fp6_sqr(&tmp2_fp6,&tmp1_efp6.x);
    fp6_add(&tmp3_fp6,&tmp2_fp6,&tmp2_fp6);
    fp6_add(&tmp2_fp6,&tmp2_fp6,&tmp3_fp6);
    fp6_mul(&tmp3_fp6,&tmp1_fp6,&tmp2_fp6);
    
    fp6_sqr(&tmp1_fp6,&tmp3_fp6);
    fp6_add(&tmp2_fp6,&tmp1_efp6.x,&tmp1_efp6.x);
    fp6_sub(&ANS->x,&tmp1_fp6,&tmp2_fp6);
    
    fp6_sub(&tmp1_fp6,&tmp1_efp6.x,&ANS->x);
    fp6_mul(&tmp2_fp6,&tmp3_fp6,&tmp1_fp6);
    fp6_sub(&ANS->y,&tmp2_fp6,&tmp1_efp6.y);
}

void efp6_eca(efp6_t *ANS,efp6_t *P1,efp6_t *P2){
	static efp6_t tmp1_efp6,tmp2_efp6;
	static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
    if(P1->infinity==1){
        efp6_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp6_set(ANS,P1);
        return;
    }else if(fp6_cmp(&P1->x,&P2->x)==0){
        if(fp6_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp6_ecd(ANS,P1);
            return;
        }
    }
    
    efp6_set(&tmp1_efp6,P1);
    efp6_set(&tmp2_efp6,P2);
    
    fp6_sub(&tmp1_fp6,&tmp2_efp6.x,&tmp1_efp6.x);
    fp6_inv(&tmp1_fp6,&tmp1_fp6);
    fp6_sub(&tmp2_fp6,&tmp2_efp6.y,&tmp1_efp6.y);
    fp6_mul(&tmp3_fp6,&tmp1_fp6,&tmp2_fp6);
    fp6_sqr(&tmp1_fp6,&tmp3_fp6);
    fp6_sub(&tmp2_fp6,&tmp1_fp6,&tmp1_efp6.x);
    fp6_sub(&ANS->x,&tmp2_fp6,&tmp2_efp6.x);
    fp6_sub(&tmp1_fp6,&tmp1_efp6.x,&ANS->x);
    fp6_mul(&tmp2_fp6,&tmp3_fp6,&tmp1_fp6);
    fp6_sub(&ANS->y,&tmp2_fp6,&tmp1_efp6.y);
}

void efp6_scm(efp6_t *ANS,efp6_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp6_set(ANS,P);
        return;
    }
    
    efp6_t Tmp_P,Next_P;
    efp6_init(&Tmp_P);
    efp6_set(&Tmp_P,P);
    efp6_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    efp6_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp6_ecd(&Next_P,&Next_P);
        if(binary[i]=='1'){
            efp6_eca(&Next_P,&Next_P,&Tmp_P);
        }
    }
    efp6_set(ANS,&Next_P);
}
