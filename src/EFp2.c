#include <ELiPS/EFp2.h>
//EFp2
void EFp2_init(EFp2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    P->infinity=0;
}

void EFpZ2_init(EFpZ2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    Fp2_init(&P->z);
    P->infinity=0;
}
void EFpZT2_init(EFpZT2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    Fp2_init(&P->z);
    Fp2_init(&P->zz);
    Fp2_init(&P->zzz);
    P->infinity=0;
}
void EFp2_printf(char *str,EFp2 *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp2_printf("",&P->x);
        printf(",");
        Fp2_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}
void EFpZ2_printf(char *str,EFpZ2 *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp2_printf("",&P->x);
        printf(",");
        Fp2_printf("",&P->y);
        printf(",");
        Fp2_printf("",&P->z);
        printf(")");
    }else{
        printf("0");
    }
}
void EFp2_set(EFp2 *ANS,EFp2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void EFpZ2_set(EFpZ2 *ANS,EFpZ2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    Fp2_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpZT2_set(EFpZT2 *ANS,EFpZT2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    Fp2_set(&ANS->z,&A->z);
    Fp2_set(&ANS->zz,&A->zz);
    Fp2_set(&ANS->zzz,&A->zzz);
    ANS->infinity=A->infinity;
}
void EFpZT2_to_EFpZ2(EFpZ2 *ANS,EFpZT2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    Fp2_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpZ2_to_EFpZT2(EFpZT2 *ANS,EFpZ2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    Fp2_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFp2_to_EFpZ2(EFpZ2 *ANS,EFp2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set(&ANS->y,&A->y);
    Fp2_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}
void EFp2_Jacobian(EFp2 *ANS,EFpZ2 *A){
    static Fp2 Zi,Zt;
    Fp2_inv(&Zi,&A->z);
    Fp2_mul(&Zt,&Zi,&Zi);
    Fp2_mul(&ANS->x,&A->x,&Zt);
    Fp2_mul(&Zt,&Zt,&Zi);
    Fp2_mul(&ANS->y,&A->y,&Zt);
    ANS->infinity=A->infinity;
}

void EFp2_set_ui(EFp2 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x,UI);
    Fp2_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp2_set_mpn(EFp2 *ANS,mp_limb_t *A){
    Fp2_set_mpn(&ANS->x,A);
    Fp2_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}

void EFp2_set_neg(EFp2 *ANS,EFp2 *A){
    Fp2_set(&ANS->x,&A->x);
    Fp2_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}


int  EFp2_cmp(EFp2 *A,EFp2 *B){
    if(Fp2_cmp(&A->x,&B->x)==0 && Fp2_cmp(&A->y,&B->y)==0){
        return 0;   
    }else if(A->infinity==1&&B->infinity==1){
	return 0;
    }else{
    return 1;
    }
}

void EFp2_rational_point(EFp2 *P){
    Fp2 tmp1,tmp2;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp2_set_random(&P->x,state);
        Fp2_sqr(&tmp1,&P->x);
        Fp2_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpn(&tmp2.x0,&tmp2.x0,curve_b);
        if(Fp2_legendre(&tmp2)==1){
            Fp2_sqrt(&P->y,&tmp2);
            break;
        }
    }
}

void EFp2_ECD(EFp2 *ANS,EFp2 *P){
    static EFp2 tmp1_EFp2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp2_set(&tmp1_EFp2,P);
    
    Fp2_add(&tmp1_Fp2,&tmp1_EFp2.y,&tmp1_EFp2.y);
    
    Fp2_inv(&tmp1_Fp2,&tmp1_Fp2);
    Fp2_sqr(&tmp2_Fp2,&tmp1_EFp2.x);
    Fp2_add(&tmp3_Fp2,&tmp2_Fp2,&tmp2_Fp2);
    Fp2_add(&tmp2_Fp2,&tmp2_Fp2,&tmp3_Fp2);
    Fp2_mul(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);
    
    Fp2_sqr(&tmp1_Fp2,&tmp3_Fp2);
    Fp2_add(&tmp2_Fp2,&tmp1_EFp2.x,&tmp1_EFp2.x);
    Fp2_sub(&ANS->x,&tmp1_Fp2,&tmp2_Fp2);
    
    Fp2_sub(&tmp1_Fp2,&tmp1_EFp2.x,&ANS->x);
    Fp2_mul(&tmp2_Fp2,&tmp3_Fp2,&tmp1_Fp2);
    Fp2_sub(&ANS->y,&tmp2_Fp2,&tmp1_EFp2.y);
}

void EFp2_ECD_Jacobian(EFpZ2 *ANS,EFpZ2 *P){
    static EFpZ2 Pt;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    static Fp2 tmpy_Fp2;
    static Fp2 s,m,T;
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFpZ2_set(&Pt,P);
    
    //s
    Fp2_sqr(&tmpy_Fp2,&Pt.y);
    Fp2_mul(&tmp1_Fp2,&tmpy_Fp2,&Pt.x);
    Fp2_lshift(&s,&tmp1_Fp2,2);

    //m
    Fp2_sqr(&tmp1_Fp2,&Pt.x);
    Fp2_lshift(&tmp2_Fp2,&tmp1_Fp2,1);
    Fp2_add(&m,&tmp2_Fp2,&tmp1_Fp2);

    //T
    Fp2_sqr(&T,&m);
    Fp2_add(&tmp1_Fp2,&s,&s);
    Fp2_sub(&T,&T,&tmp1_Fp2);


    //ANS->x
    Fp2_set(&ANS->x,&T);

    //ANS->y
    Fp2_sub(&tmp1_Fp2,&s,&T);
    Fp2_mul(&ANS->y,&tmp1_Fp2,&m);

    Fp2_sqr(&tmp1_Fp2,&tmpy_Fp2);
    Fp2_lshift(&tmp1_Fp2,&tmp1_Fp2,3);
    Fp2_sub(&ANS->y,&ANS->y,&tmp1_Fp2);

    
    //ANS->z
    Fp2_add(&tmp1_Fp2,&Pt.y,&Pt.y);
    Fp2_mul(&ANS->z,&tmp1_Fp2,&Pt.z);

}
void EFp2_ECD_lazy(EFp2 *ANS,EFp2 *P){
    static EFp2 tmp1_EFp2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp2_set(&tmp1_EFp2,P);
    
    Fp2_add(&tmp1_Fp2,&tmp1_EFp2.y,&tmp1_EFp2.y);
    
    Fp2_inv(&tmp1_Fp2,&tmp1_Fp2);
    Fp2_sqr(&tmp2_Fp2,&tmp1_EFp2.x);
    Fp2_add_lazy(&tmp3_Fp2,&tmp2_Fp2,&tmp2_Fp2);
    Fp2_add_lazy(&tmp2_Fp2,&tmp2_Fp2,&tmp3_Fp2);
    Fp2_mul_lazy(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);
    
    Fp2_sqr(&tmp1_Fp2,&tmp3_Fp2);
    Fp2_add(&tmp2_Fp2,&tmp1_EFp2.x,&tmp1_EFp2.x);
    Fp2_sub(&ANS->x,&tmp1_Fp2,&tmp2_Fp2);
    
    Fp2_sub_lazy(&tmp1_Fp2,&tmp1_EFp2.x,&ANS->x);
    Fp2_mul_lazy(&tmp2_Fp2,&tmp3_Fp2,&tmp1_Fp2);
    Fp2_sub(&ANS->y,&tmp2_Fp2,&tmp1_EFp2.y);
}

void EFp2_ECD_Jacobian_lazy(EFpZ2 *ANS,EFpZ2 *P){
    static EFpZ2 Pt;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    static Fp2 tmpy_Fp2;
    static Fp2 s,m,T;
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFpZ2_set(&Pt,P);
    
    //s
    Fp2_sqr_lazy(&tmpy_Fp2,&Pt.y);
    Fp2_mul_lazy(&tmp1_Fp2,&tmpy_Fp2,&Pt.x);
    Fp2_lshift(&s,&tmp1_Fp2,2);

    //m
    Fp2_sqr_lazy(&tmp1_Fp2,&Pt.x);
    Fp2_lshift(&tmp2_Fp2,&tmp1_Fp2,1);
    Fp2_add(&m,&tmp2_Fp2,&tmp1_Fp2);

    //T
    Fp2_sqr_lazy(&T,&m);
    Fp2_add(&tmp1_Fp2,&s,&s);
    Fp2_sub(&T,&T,&tmp1_Fp2);

    //ANS->x
    Fp2_set(&ANS->x,&T);

    //ANS->y
    Fp2_sub_lazy(&tmp1_Fp2,&s,&T);
    Fp2_mul_lazy(&ANS->y,&tmp1_Fp2,&m);

    Fp2_sqr_lazy(&tmp1_Fp2,&tmpy_Fp2);
    Fp2_lshift(&tmp1_Fp2,&tmp1_Fp2,3);
    Fp2_sub(&ANS->y,&ANS->y,&tmp1_Fp2);

    //ANS->z
    Fp2_add_lazy(&tmp1_Fp2,&Pt.y,&Pt.y);
    Fp2_mul(&ANS->z,&tmp1_Fp2,&Pt.z);

}
void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2){
    static EFp2 tmp1_EFp2,tmp2_EFp2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    if(P1->infinity==1){
        EFp2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD(ANS,P1);
            return;
        }
    }
    
    EFp2_set(&tmp1_EFp2,P1);
    EFp2_set(&tmp2_EFp2,P2);
    
    Fp2_sub(&tmp1_Fp2,&tmp2_EFp2.x,&tmp1_EFp2.x);
    Fp2_inv(&tmp1_Fp2,&tmp1_Fp2);
    Fp2_sub(&tmp2_Fp2,&tmp2_EFp2.y,&tmp1_EFp2.y);
    Fp2_mul(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);
    Fp2_sqr(&tmp1_Fp2,&tmp3_Fp2);
    Fp2_sub(&tmp2_Fp2,&tmp1_Fp2,&tmp1_EFp2.x);
    Fp2_sub(&ANS->x,&tmp2_Fp2,&tmp2_EFp2.x);
    Fp2_sub(&tmp1_Fp2,&tmp1_EFp2.x,&ANS->x);
    Fp2_mul(&tmp2_Fp2,&tmp3_Fp2,&tmp1_Fp2);
    Fp2_sub(&ANS->y,&tmp2_Fp2,&tmp1_EFp2.y);
}
void EFp2_ECA_Jacobian(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2){
    static EFpZ2 Pt1,Pt2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    static Fp2 tmpZ1_Fp2,tmpZ2_Fp2,tmpH2_Fp2,tmpH3_Fp2,tmpU1H2_Fp2;
    static Fp2 U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        EFpZ2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->z,&P2->z)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0&&Fp2_cmp(&P1->z,&P2->z)==0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD_Jacobian(ANS,P1);
            return;
        }
    }
    
    EFpZ2_set(&Pt1,P1);
    EFpZ2_set(&Pt2,P2);
    
    //U1
    Fp2_sqr(&tmpZ2_Fp2,&Pt2.z);
    Fp2_mul(&U1,&tmpZ2_Fp2,&Pt1.x);

    //U2
    Fp2_sqr(&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul(&U2,&tmpZ1_Fp2,&Pt2.x);

    //S1
    Fp2_mul(&tmp1_Fp2,&tmpZ2_Fp2,&Pt2.z);
    Fp2_mul(&S1,&tmp1_Fp2,&Pt1.y);

    //S2
    Fp2_mul(&tmp1_Fp2,&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul(&S2,&tmp1_Fp2,&Pt2.y);

    //H
    Fp2_sub(&H,&U2,&U1);
    //r
    Fp2_sub(&r,&S2,&S1);

    //ANS->x
    Fp2_sqr(&ANS->x,&r);

    Fp2_sqr(&tmpH2_Fp2,&H);
    Fp2_mul(&tmpH3_Fp2,&tmpH2_Fp2,&H);
    Fp2_sub(&ANS->x,&ANS->x,&tmpH3_Fp2);

    Fp2_mul(&tmpU1H2_Fp2,&tmpH2_Fp2,&U1);
    Fp2_add(&tmp1_Fp2,&tmpU1H2_Fp2,&tmpU1H2_Fp2);
    Fp2_sub(&ANS->x,&ANS->x,&tmp1_Fp2);


    //ANS->y
    Fp2_sub(&tmp1_Fp2,&tmpU1H2_Fp2,&ANS->x);
    Fp2_mul(&ANS->y,&tmp1_Fp2,&r);

    Fp2_mul(&tmp1_Fp2,&tmpH3_Fp2,&S1);
    Fp2_sub(&ANS->y,&ANS->y,&tmp1_Fp2);

    
    //ANS->z
    Fp2_mul(&tmp1_Fp2,&Pt1.z,&Pt2.z);
    Fp2_mul(&ANS->z,&tmp1_Fp2,&H);
}

void EFp2_ECA_Jacobian_scm(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2){
    static EFpZ2 Pt1,Pt2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    static Fp2 tmpZ1_Fp2,tmpZ2_Fp2,tmpH2_Fp2,tmpH3_Fp2,tmpU1H2_Fp2;
    static Fp2 U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        EFpZ2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->z,&P2->z)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0&&Fp2_cmp(&P1->z,&P2->z)==0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD_Jacobian(ANS,P1);
            return;
        }
    }
    
    EFpZ2_set(&Pt1,P1);
    EFpZ2_set(&Pt2,P2);
    
    //U1
    Fp2_set(&U1,&Pt1.x);

    //U2
    Fp2_sqr(&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul(&U2,&tmpZ1_Fp2,&Pt2.x);

    //S1
    Fp2_set(&S1,&Pt1.y);

    //S2
    Fp2_mul(&tmp1_Fp2,&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul(&S2,&tmp1_Fp2,&Pt2.y);

    //H
    Fp2_sub(&H,&U2,&U1);
    //r
    Fp2_sub(&r,&S2,&S1);

    //ANS->x
    Fp2_sqr(&ANS->x,&r);

    Fp2_sqr(&tmpH2_Fp2,&H);
    Fp2_mul(&tmpH3_Fp2,&tmpH2_Fp2,&H);
    Fp2_sub(&ANS->x,&ANS->x,&tmpH3_Fp2);

    Fp2_mul(&tmpU1H2_Fp2,&tmpH2_Fp2,&U1);
    Fp2_add(&tmp1_Fp2,&tmpU1H2_Fp2,&tmpU1H2_Fp2);
    Fp2_sub(&ANS->x,&ANS->x,&tmp1_Fp2);


    //ANS->y
    Fp2_sub(&tmp1_Fp2,&tmpU1H2_Fp2,&ANS->x);
    Fp2_mul(&ANS->y,&tmp1_Fp2,&r);

    Fp2_mul(&tmp1_Fp2,&tmpH3_Fp2,&S1);
    Fp2_sub(&ANS->y,&ANS->y,&tmp1_Fp2);

    
    //ANS->z
    Fp2_mul(&ANS->z,&Pt1.z,&H);
}
void EFp2_ECA_lazy(EFp2 *ANS,EFp2 *P1,EFp2 *P2){
    static EFp2 tmp1_EFp2,tmp2_EFp2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    if(P1->infinity==1){
        EFp2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD_lazy(ANS,P1);
            return;
        }
    }
    
    EFp2_set(&tmp1_EFp2,P1);
    EFp2_set(&tmp2_EFp2,P2);
    
    Fp2_sub_lazy(&tmp1_Fp2,&tmp2_EFp2.x,&tmp1_EFp2.x);
    Fp2_inv(&tmp1_Fp2,&tmp1_Fp2);
    Fp2_sub_lazy(&tmp2_Fp2,&tmp2_EFp2.y,&tmp1_EFp2.y);
    Fp2_mul_lazy(&tmp3_Fp2,&tmp1_Fp2,&tmp2_Fp2);
    Fp2_sqr(&tmp1_Fp2,&tmp3_Fp2);
    Fp2_sub_lazy(&tmp2_Fp2,&tmp1_Fp2,&tmp1_EFp2.x);
    Fp2_sub(&ANS->x,&tmp2_Fp2,&tmp2_EFp2.x);
    Fp2_sub_lazy(&tmp1_Fp2,&tmp1_EFp2.x,&ANS->x);
    Fp2_mul_lazy(&tmp2_Fp2,&tmp3_Fp2,&tmp1_Fp2);
    Fp2_sub(&ANS->y,&tmp2_Fp2,&tmp1_EFp2.y);
}
void EFp2_ECA_Jacobian_lazy(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2){
    static EFpZ2 Pt1,Pt2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    static Fp2 tmpZ1_Fp2,tmpZ2_Fp2,tmpH2_Fp2,tmpH3_Fp2,tmpU1H2_Fp2;
    static Fp2 U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        EFpZ2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->z,&P2->z)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0&&Fp2_cmp(&P1->z,&P2->z)==0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD_Jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    EFpZ2_set(&Pt1,P1);
    EFpZ2_set(&Pt2,P2);
    
    //U1
    Fp2_sqr_lazy(&tmpZ2_Fp2,&Pt2.z);
    Fp2_mul_lazy(&U1,&tmpZ2_Fp2,&Pt1.x);

    //U2
    Fp2_sqr_lazy(&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul_lazy(&U2,&tmpZ1_Fp2,&Pt2.x);

    //S1
    Fp2_mul_lazy(&tmp1_Fp2,&tmpZ2_Fp2,&Pt2.z);
    Fp2_mul_lazy(&S1,&tmp1_Fp2,&Pt1.y);

    //S2
    Fp2_mul_lazy(&tmp1_Fp2,&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul_lazy(&S2,&tmp1_Fp2,&Pt2.y);

    //H
    Fp2_sub(&H,&U2,&U1);
    //r
    Fp2_sub(&r,&S2,&S1);

    //ANS->x
    Fp2_sqr_lazy(&ANS->x,&r);

    Fp2_sqr_lazy(&tmpH2_Fp2,&H);
    Fp2_mul_lazy(&tmpH3_Fp2,&tmpH2_Fp2,&H);
    Fp2_sub(&ANS->x,&ANS->x,&tmpH3_Fp2);

    Fp2_mul_lazy(&tmpU1H2_Fp2,&tmpH2_Fp2,&U1);
    Fp2_add(&tmp1_Fp2,&tmpU1H2_Fp2,&tmpU1H2_Fp2);
    Fp2_sub(&ANS->x,&ANS->x,&tmp1_Fp2);


    //ANS->y
    Fp2_sub_lazy(&tmp1_Fp2,&tmpU1H2_Fp2,&ANS->x);
    Fp2_mul_lazy(&ANS->y,&tmp1_Fp2,&r);

    Fp2_mul_lazy(&tmp1_Fp2,&tmpH3_Fp2,&S1);
    Fp2_sub(&ANS->y,&ANS->y,&tmp1_Fp2);

    
    //ANS->z
    Fp2_mul_lazy(&tmp1_Fp2,&Pt1.z,&Pt2.z);
    Fp2_mul_lazy(&ANS->z,&tmp1_Fp2,&H);
}
void EFp2_ECA_Jacobian_table(EFpZ2 *ANS,EFpZ2 *P1,EFpZT2 *P2){
    static EFpZ2 Pt1;
    static EFpZT2 Pt2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    static Fp2 tmpZ1_Fp2,tmpZ2_Fp2,tmpH2_Fp2,tmpH3_Fp2,tmpU1H2_Fp2;
    static Fp2 U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        EFpZT2_to_EFpZ2(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->z,&P2->z)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0&&Fp2_cmp(&P1->z,&P2->z)==0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD_Jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    EFpZ2_set(&Pt1,P1);
    EFpZT2_set(&Pt2,P2);
    
    //U1
    Fp2_mul_lazy(&U1,&Pt2.zz,&Pt1.x);

    //U2
    Fp2_sqr_lazy(&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul_lazy(&U2,&tmpZ1_Fp2,&Pt2.x);

    //S1
    Fp2_mul_lazy(&S1,&Pt2.zzz,&Pt1.y);

    //S2
    Fp2_mul_lazy(&tmp1_Fp2,&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul_lazy(&S2,&tmp1_Fp2,&Pt2.y);

    //H
    Fp2_sub(&H,&U2,&U1);
    //r
    Fp2_sub(&r,&S2,&S1);

    //ANS->x
    Fp2_sqr_lazy(&ANS->x,&r);

    Fp2_sqr_lazy(&tmpH2_Fp2,&H);
    Fp2_mul_lazy(&tmpH3_Fp2,&tmpH2_Fp2,&H);
    Fp2_sub(&ANS->x,&ANS->x,&tmpH3_Fp2);

    Fp2_mul_lazy(&tmpU1H2_Fp2,&tmpH2_Fp2,&U1);
    Fp2_add(&tmp1_Fp2,&tmpU1H2_Fp2,&tmpU1H2_Fp2);
    Fp2_sub(&ANS->x,&ANS->x,&tmp1_Fp2);


    //ANS->y
    Fp2_sub_lazy(&tmp1_Fp2,&tmpU1H2_Fp2,&ANS->x);
    Fp2_mul_lazy(&ANS->y,&tmp1_Fp2,&r);

    Fp2_mul_lazy(&tmp1_Fp2,&tmpH3_Fp2,&S1);
    Fp2_sub(&ANS->y,&ANS->y,&tmp1_Fp2);

    
    //ANS->z
    Fp2_mul_lazy(&tmp1_Fp2,&Pt1.z,&Pt2.z);
    Fp2_mul_lazy(&ANS->z,&tmp1_Fp2,&H);
}
void EFp2_ECA_Jacobian_scm_lazy(EFpZ2 *ANS,EFpZ2 *P1,EFpZ2 *P2){
    static EFpZ2 Pt1,Pt2;
    static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
    static Fp2 tmpZ1_Fp2,tmpZ2_Fp2,tmpH2_Fp2,tmpH3_Fp2,tmpU1H2_Fp2;
    static Fp2 U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        EFpZ2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->z,&P2->z)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0&&Fp2_cmp(&P1->z,&P2->z)==0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD_Jacobian(ANS,P1);
            return;
        }
    }
    
    EFpZ2_set(&Pt1,P1);
    EFpZ2_set(&Pt2,P2);
    
    //U1
    Fp2_set(&U1,&Pt1.x);

    //U2
    Fp2_sqr_lazy(&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul_lazy(&U2,&tmpZ1_Fp2,&Pt2.x);

    //S1
    Fp2_set(&S1,&Pt1.y);

    //S2
    Fp2_mul_lazy(&tmp1_Fp2,&tmpZ1_Fp2,&Pt1.z);
    Fp2_mul_lazy(&S2,&tmp1_Fp2,&Pt2.y);

    //H
    Fp2_sub_lazy(&H,&U2,&U1);
    //r
    Fp2_sub_lazy(&r,&S2,&S1);

    //ANS->x
    Fp2_sqr_lazy(&ANS->x,&r);

    Fp2_sqr_lazy(&tmpH2_Fp2,&H);
    Fp2_mul_lazy(&tmpH3_Fp2,&tmpH2_Fp2,&H);
    Fp2_sub_lazy(&ANS->x,&ANS->x,&tmpH3_Fp2);

    Fp2_mul_lazy(&tmpU1H2_Fp2,&tmpH2_Fp2,&U1);
    Fp2_add_lazy(&tmp1_Fp2,&tmpU1H2_Fp2,&tmpU1H2_Fp2);
    Fp2_sub(&ANS->x,&ANS->x,&tmp1_Fp2);


    //ANS->y
    Fp2_sub_lazy(&tmp1_Fp2,&tmpU1H2_Fp2,&ANS->x);
    Fp2_mul_lazy(&ANS->y,&tmp1_Fp2,&r);

    Fp2_mul_lazy(&tmp1_Fp2,&tmpH3_Fp2,&S1);
    Fp2_sub(&ANS->y,&ANS->y,&tmp1_Fp2);

    
    //ANS->z
    Fp2_mul(&ANS->z,&Pt1.z,&H);
}
void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    EFp2 Tmp_P,Next_P;
    EFp2_init(&Tmp_P);
    EFp2_set(&Tmp_P,P);
    EFp2_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFp2_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp2_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp2_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp2_set(ANS,&Next_P);
}
void EFp2_SCM_Jacobian(EFp2 *ANS,EFp2 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    EFpZ2 Tmp_P,Next_P;
    EFp2_to_EFpZ2(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFpZ2_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp2_ECD_Jacobian(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp2_ECA_Jacobian_scm(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp2_Jacobian(ANS,&Next_P);
}
void EFp2_SCM_lazy(EFp2 *ANS,EFp2 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    EFp2 Tmp_P,Next_P;
    EFp2_init(&Tmp_P);
    EFp2_set(&Tmp_P,P);
    EFp2_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFp2_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp2_ECD_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp2_ECA_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp2_set(ANS,&Next_P);
}
void EFp2_SCM_Jacobian_lazy(EFp2 *ANS,EFp2 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    EFpZ2 Tmp_P,Next_P;
    EFp2_to_EFpZ2(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFpZ2_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp2_ECD_Jacobian_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp2_ECA_Jacobian_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp2_Jacobian(ANS,&Next_P);
}


//skew_frobenius_map
void EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
    Fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p1][1]);
    //y
    Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
    Fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p1][4]);
}

void EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A){
    //x
    Fp2_mul(&ANS->x,&A->x,&frobenius_constant[f_p2][1]);
    //y
    Fp2_mul(&ANS->y,&A->y,&frobenius_constant[f_p2][4]);
}

void EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
    Fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p3][1]);
    //y
    Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
    Fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p3][4]);
}

void EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A){
    //x
    Fp2_mul(&ANS->x,&A->x,&frobenius_constant[f_p10][1]);
    //y
    Fp2_mul(&ANS->y,&A->y,&frobenius_constant[f_p10][4]);
}
