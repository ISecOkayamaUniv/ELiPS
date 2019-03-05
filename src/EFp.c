#include <ELiPS/EFp.h>
//EFp
void EFp_init(EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->infinity=0;
}
void EFpZ_init(EFpZ *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    Fp_init(&P->z);
    P->infinity=0;
}
void EFpZT_init(EFpZT *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    Fp_init(&P->z);
    Fp_init(&P->zz);
    Fp_init(&P->zzz);
    P->infinity=0;
}
void EFp_printf(EFp *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}
void EFpZ_printf(EFpZ *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(",");
        Fp_printf(&P->z,"");
        printf(")");
    }else{
        printf("Infinity");
    }
}
void EFpZT_printf(EFpZT *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(",");
        Fp_printf(&P->z,"");
        printf(",");
        Fp_printf(&P->zz,"");
        printf(",");
        Fp_printf(&P->zzz,"");
        printf(")");
    }else{
        printf("Infinity");
    }
}
void EFp_set(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void EFpZ_set(EFpZ *ANS,EFpZ *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpZT_set(EFpZT *ANS,EFpZT *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    Fp_set(&ANS->zz,&A->zz);
    Fp_set(&ANS->zzz,&A->zzz);
    ANS->infinity=A->infinity;
}
void EFp_to_EFpZ(EFpZ *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}
void EFpZT_to_EFpZ(EFpZ *ANS,EFpZT *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpZ_to_EFpZT(EFpZT *ANS,EFpZ *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFp_Jacobian(EFp *ANS,EFpZ *A){
    static Fp Zi,Zt;
    Fp_inv(&Zi,&A->z);
    Fp_mul(&Zt,&Zi,&Zi);
    Fp_mul(&ANS->x,&A->x,&Zt);
    Fp_mul(&Zt,&Zt,&Zi);
    Fp_mul(&ANS->y,&A->y,&Zt);
    ANS->infinity=A->infinity;
}
void EFp_set_ui(EFp *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x,UI);
    Fp_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}
void EFp_set_mpn(EFp *ANS,mp_limb_t *A){
    Fp_set_mpn(&ANS->x,A);
    Fp_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}
void EFp_set_neg(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}


int  EFp_cmp(EFp *A,EFp *B){
    if(Fp_cmp(&A->x,&B->x)==0 && Fp_cmp(&A->y,&B->y)==0){
        return 0;   
    }else if(A->infinity==1&&B->infinity==1){
	return 0;
    }else{
    return 1;
    }
}
void EFp_rational_point(EFp *P){
    Fp tmp1,tmp2,tmp_x;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp_x);
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
	
    while(1){
        Fp_set_random(&P->x,state);
        Fp_mul(&tmp1,&P->x,&P->x);
        Fp_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpn(&tmp_x,&tmp2,curve_b);
        if(Fp_legendre(&tmp_x)==1){
            Fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
}

void EFp_ECD(EFp *ANS,EFp *P){
    static EFp tmp1_EFp;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    EFp_set(&tmp1_EFp,P);
    
    Fp_add(&tmp1_Fp,&tmp1_EFp.y,&tmp1_EFp.y);
    Fp_inv(&tmp1_Fp,&tmp1_Fp);


    Fp_mul(&tmp2_Fp,&tmp1_EFp.x,&tmp1_EFp.x);
    Fp_add(&tmp3_Fp,&tmp2_Fp,&tmp2_Fp);
    Fp_add(&tmp2_Fp,&tmp2_Fp,&tmp3_Fp);

    Fp_mul(&tmp3_Fp,&tmp1_Fp,&tmp2_Fp);
    Fp_mul(&tmp1_Fp,&tmp3_Fp,&tmp3_Fp);

    Fp_add(&tmp2_Fp,&tmp1_EFp.x,&tmp1_EFp.x);
    Fp_sub(&ANS->x,&tmp1_Fp,&tmp2_Fp);

    Fp_sub(&tmp1_Fp,&tmp1_EFp.x,&ANS->x);
    Fp_mul(&tmp2_Fp,&tmp3_Fp,&tmp1_Fp);
    Fp_sub(&ANS->y,&tmp2_Fp,&tmp1_EFp.y);
}

void EFp_ECD_Jacobian(EFpZ *ANS,EFpZ *P){
    static EFpZ Pt;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    static Fp s,m,T;
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFpZ_set(&Pt,P);
    
    //s
    Fp_mul_ui(&tmp1_Fp,&Pt.x,4);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&Pt.y);
    Fp_mul(&s,&tmp1_Fp,&Pt.y);

    //m
    Fp_mul_ui(&tmp2_Fp,&Pt.x,3);
    Fp_mul(&m,&tmp2_Fp,&Pt.x);

    //T
    Fp_add(&tmp1_Fp,&s,&s);
    Fp_set_neg(&T,&tmp1_Fp);
    Fp_mul(&tmp1_Fp,&m,&m);
    Fp_add(&T,&T,&tmp1_Fp);


    //ANS->x
    Fp_set(&ANS->x,&T);

    //ANS->y
    Fp_mul_ui(&tmp1_Fp,&Pt.y,8);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&Pt.y);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&Pt.y);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&Pt.y);
    Fp_set_neg(&ANS->y,&tmp1_Fp);

    Fp_sub(&tmp1_Fp,&s,&T);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&m);
    Fp_add(&ANS->y,&ANS->y,&tmp1_Fp);
    
    //ANS->z
    Fp_add(&tmp1_Fp,&Pt.y,&Pt.y);
    Fp_mul(&ANS->z,&tmp1_Fp,&Pt.z);
}

void EFp_ECD_lazy(EFp *ANS,EFp *P){
    static Fp tmpF2,tmpF3;
    static mp_limb_t bufL[FPLIMB2],tmpL2[FPLIMB2],tmpL3[FPLIMB2],tmpL5[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB],tmp4[FPLIMB];
    static EFp P1t;
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    EFp_set(&P1t,P);
    
    Lazy_add(tmp1,FPLIMB,P1t.y.x0,FPLIMB,P1t.y.x0,FPLIMB);
    mpn_invert(tmp1,tmp1,prime);

    Lazy_mul(tmpL2,P1t.x.x0,P1t.x.x0);
    Lazy_add(tmpL3,FPLIMB2,tmpL2,FPLIMB2,tmpL2,FPLIMB2);
    Lazy_add(tmpL2,FPLIMB2,tmpL3,FPLIMB2,tmpL2,FPLIMB2);//3x^2
    mpn_mod(&tmpF2,tmpL2,FPLIMB2);

    Lazy_mul(tmpL3,tmp1,tmpF2.x0);
    mpn_mod(&tmpF3,tmpL3,FPLIMB2);
    Lazy_mul(tmpL2,tmpF3.x0,tmpF3.x0);

    Lazy_add(tmp4,FPLIMB,P1t.x.x0,FPLIMB,P1t.x.x0,FPLIMB);
    Lazy_sub(bufL,FPLIMB2,tmpL2,FPLIMB2,tmp4,FPLIMB);
    mpn_mod(&ANS->x,bufL,FPLIMB2);

    Lazy_sub(tmp1,FPLIMB,P1t.x.x0,FPLIMB,ANS->x.x0,FPLIMB);
    Lazy_mul(tmpL2,tmpF3.x0,tmp1);
    Lazy_sub(bufL,FPLIMB2,tmpL2,FPLIMB2,P1t.y.x0,FPLIMB);
    mpn_mod(&ANS->y,bufL,FPLIMB2);

}

void EFp_ECD_Jacobian_lazy(EFpZ *ANS,EFpZ *P){
    static mp_limb_t s[FPLIMB],m[FPLIMB],T[FPLIMB];

    static Fp bufF;
    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpy[FPLIMB];
    static EFpZ Pt;
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFpZ_set(&Pt,P);
    
    //s
    Lazy_sqr(bufL,Pt.y.x0);
    Lazy_mod(tmpy,bufL,FPLIMB2);

    Lazy_mul(bufL,tmpy,Pt.x.x0);
    Lazy_mod(tmp1,bufL,FPLIMB2);
    
    Lazy_add(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Lazy_add(s,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);


    //m
    Lazy_add(tmp1,FPLIMB,Pt.x.x0,FPLIMB,Pt.x.x0,FPLIMB);
    Lazy_add(tmp1,FPLIMB,tmp1,FPLIMB,Pt.x.x0,FPLIMB);
    Lazy_mul(bufL,tmp1,Pt.x.x0);
    Lazy_mod(m,bufL,FPLIMB2);

    //T
    Lazy_sqr(bufL,m);
    Lazy_mod(T,bufL,FPLIMB2);//heraseru
    Lazy_add(tmp1,FPLIMB,s,FPLIMB,s,FPLIMB);
    Lazy_sub(buf,FPLIMB,T,FPLIMB,tmp1,FPLIMB);
    Lazy_mod(T,buf,FPLIMB);

    //ANS->x
    mpn_copyd(ANS->x.x0,T,FPLIMB);

    //ANS->y
    Lazy_sub(tmp1,FPLIMB,s,FPLIMB,T,FPLIMB);
    Lazy_mul(bufL,tmp1,m);

    Lazy_sqr(tmpL1,tmpy);
    Lazy_mod(tmp1,tmpL1,FPLIMB2);
    Lazy_add(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Lazy_add(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Lazy_add(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Lazy_sub(bufL,FPLIMB,bufL,FPLIMB2,tmp1,FPLIMB);

    mpn_mod(&ANS->y,bufL,FPLIMB2);

    
    //ANS->z
    Lazy_add(tmp1,FPLIMB,Pt.y.x0,FPLIMB,Pt.y.x0,FPLIMB);
    Lazy_mul(bufL,tmp1,Pt.z.x0);
    mpn_mod(&ANS->z,bufL,FPLIMB2);
}
void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2){
    static EFp tmp1_EFp,tmp2_EFp;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    if(P1->infinity==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD(ANS,P1);
            return;
        }
    }
    
    EFp_set(&tmp1_EFp,P1);
    EFp_set(&tmp2_EFp,P2);
    
    Fp_sub(&tmp1_Fp,&tmp2_EFp.x,&tmp1_EFp.x);
    Fp_inv(&tmp1_Fp,&tmp1_Fp);
    Fp_sub(&tmp2_Fp,&tmp2_EFp.y,&tmp1_EFp.y);
    Fp_mul(&tmp3_Fp,&tmp1_Fp,&tmp2_Fp);
    Fp_mul(&tmp1_Fp,&tmp3_Fp,&tmp3_Fp);


    Fp_sub(&tmp2_Fp,&tmp1_Fp,&tmp1_EFp.x);
    Fp_sub(&ANS->x,&tmp2_Fp,&tmp2_EFp.x);

    Fp_sub(&tmp1_Fp,&tmp1_EFp.x,&ANS->x);
    Fp_mul(&tmp2_Fp,&tmp3_Fp,&tmp1_Fp);
    Fp_sub(&ANS->y,&tmp2_Fp,&tmp1_EFp.y);
}
void EFp_ECA_Jacobian(EFpZ *ANS,EFpZ *P1,EFpZ *P2){
    static EFpZ Pt1,Pt2;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    static Fp U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        EFpZ_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0&&Fp_cmp(&P1->z,&P2->z)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0&&Fp_cmp(&P1->z,&P2->z)==0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD_Jacobian(ANS,P1);
            return;
        }
    }
    
    EFpZ_set(&Pt1,P1);
    EFpZ_set(&Pt2,P2);
    
    //U1
    Fp_mul(&tmp1_Fp,&Pt1.x,&Pt2.z);
    Fp_mul(&U1,&tmp1_Fp,&Pt2.z);

    //U2
    Fp_mul(&tmp1_Fp,&Pt2.x,&Pt1.z);
    Fp_mul(&U2,&tmp1_Fp,&Pt1.z);

    //S1
    Fp_mul(&tmp1_Fp,&Pt1.y,&Pt2.z);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&Pt2.z);
    Fp_mul(&S1,&tmp1_Fp,&Pt2.z);

    //S2
    Fp_mul(&tmp1_Fp,&Pt2.y,&Pt1.z);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&Pt1.z);
    Fp_mul(&S2,&tmp1_Fp,&Pt1.z);

    //H
    Fp_sub(&H,&U2,&U1);
    //r
    Fp_sub(&r,&S2,&S1);

    //ANS->x
    Fp_mul(&ANS->x,&r,&r);

    Fp_mul(&tmp1_Fp,&H,&H);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&H);
    Fp_sub(&ANS->x,&ANS->x,&tmp1_Fp);

    Fp_add(&tmp1_Fp,&U1,&U1);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&H);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&H);
    Fp_sub(&ANS->x,&ANS->x,&tmp1_Fp);


    //ANS->y
    Fp_mul(&tmp1_Fp,&S1,&H);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&H);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&H);
    Fp_set_neg(&ANS->y,&tmp1_Fp);

    Fp_mul(&tmp1_Fp,&U1,&H);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&H);
    Fp_sub(&tmp1_Fp,&tmp1_Fp,&ANS->x);
    Fp_mul(&tmp1_Fp,&tmp1_Fp,&r);
    Fp_add(&ANS->y,&ANS->y,&tmp1_Fp);
    
    //ANS->z
    Fp_mul(&tmp1_Fp,&Pt1.z,&Pt2.z);
    Fp_mul(&ANS->z,&tmp1_Fp,&H);
}

void EFp_ECA_lazy(EFp *ANS,EFp *P1,EFp *P2){
    static Fp Fbuf,tmpF4,out;
    static mp_limb_t bufL[FPLIMB2],tmpL3[FPLIMB2],tmpL5[FPLIMB2],tmpL6[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB],tmp2[FPLIMB];
    static EFp P1t,P2t;
    if(P1->infinity==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD_lazy(ANS,P1);
            return;
        }
    }
    
    EFp_set(&P1t,P1);
    EFp_set(&P2t,P2);

    Lazy_sub(tmp1,FPLIMB,P2t.x.x0,FPLIMB,P1t.x.x0,FPLIMB);
    mpn_invert(tmp1,tmp1,prime);
    Lazy_sub(tmp2,FPLIMB,P2t.y.x0,FPLIMB,P1t.y.x0,FPLIMB);
    Lazy_mul(tmpL3,tmp1,tmp2);
    mpn_mod(&tmpF4,tmpL3,FPLIMB2);
    Lazy_mul(tmpL5,tmpF4.x0,tmpF4.x0);


    Lazy_sub(tmpL6,FPLIMB2,tmpL5,FPLIMB2,P1t.x.x0,FPLIMB);
    Lazy_sub(ANS->x.x0,FPLIMB2,tmpL6,FPLIMB2,P2t.x.x0,FPLIMB);
    mpn_mod(&ANS->x,ANS->x.x0,FPLIMB2);

    Lazy_sub(tmp1,FPLIMB,P1t.x.x0,FPLIMB,ANS->x.x0,FPLIMB);
    Lazy_mul(tmpL3,tmpF4.x0,tmp1);
    Lazy_sub(bufL,FPLIMB2,tmpL3,FPLIMB2,P1t.y.x0,FPLIMB);
    mpn_mod(&ANS->y,bufL,FPLIMB2);

}

void EFp_ECA_Jacobian_lazy(EFpZ *ANS,EFpZ *P1,EFpZ *P2){
    static EFpZ Pt1,Pt2;
    static mp_limb_t U1[FPLIMB],U2[FPLIMB],S1[FPLIMB],S2[FPLIMB],H[FPLIMB],r[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpZ1[FPLIMB],tmpZ2[FPLIMB],tmpH2[FPLIMB],tmpH3[FPLIMB],tmpU1H2[FPLIMB];
    Fp out;
    if(P1->infinity==1){
        EFpZ_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD_Jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    EFpZ_set(&Pt1,P1);
    EFpZ_set(&Pt2,P2);

    //U1
    //Lazy_sqr(bufL,Pt2.z.x0);
    Lazy_mul(bufL,Pt2.z.x0,Pt2.z.x0);
    Lazy_mod(tmpZ2,bufL,FPLIMB2);
    Lazy_mul(tmpL1,tmpZ2,Pt1.x.x0);
    Lazy_mod(U1,tmpL1,FPLIMB2);

    //U2
    //Lazy_sqr(bufL,Pt1.z.x0);
    Lazy_mul(bufL,Pt1.z.x0,Pt1.z.x0);
    Lazy_mod(tmpZ1,bufL,FPLIMB2);
    Lazy_mul(tmpL1,tmpZ1,Pt2.x.x0);
    Lazy_mod(U2,tmpL1,FPLIMB2);

    //S1
    Lazy_mul(bufL,tmpZ2,Pt2.z.x0);
    Lazy_mod(tmp1,bufL,FPLIMB2);
    Lazy_mul(tmpL1,tmp1,Pt1.y.x0);
    Lazy_mod(S1,tmpL1,FPLIMB2);

    //S2
    Lazy_mul(bufL,tmpZ1,Pt1.z.x0);
    Lazy_mod(tmp1,bufL,FPLIMB2);
    Lazy_mul(tmpL1,tmp1,Pt2.y.x0);
    Lazy_mod(S2,tmpL1,FPLIMB2);
    
    //H
    Lazy_sub(buf,FPLIMB,U2,FPLIMB,U1,FPLIMB);
    Lazy_mod(H,buf,FPLIMB);

    //r
    Lazy_sub(buf,FPLIMB,S2,FPLIMB,S1,FPLIMB);
    Lazy_mod(r,buf,FPLIMB);

    //ANS->x
    //Lazy_sqr(bufL,r);
    Lazy_mul(bufL,r,r);

    //Lazy_sqr(tmpL1,H);
    Lazy_mul(tmpL1,H,H);
    Lazy_mod(tmpH2,tmpL1,FPLIMB2);
    Lazy_mul(tmpL1,tmpH2,H);
    Lazy_mod(tmpH3,tmpL1,FPLIMB2);
    Lazy_sub(bufL,FPLIMB2,bufL,FPLIMB2,tmpL1,FPLIMB2);

    Lazy_mul(tmpL1,tmpH2,U1);
    Lazy_mod(tmpU1H2,tmpL1,FPLIMB2);
    Lazy_add(tmpL1,FPLIMB2,tmpL1,FPLIMB2,tmpL1,FPLIMB2);
    Lazy_sub_mod(&ANS->x,bufL,tmpL1);

    //ANS->y
    Lazy_sub(tmp1,FPLIMB,tmpU1H2,FPLIMB,ANS->x.x0,FPLIMB);
    Lazy_mul(bufL,tmp1,r);

    Lazy_mul(tmpL1,tmpH3,S1);
    Lazy_sub_mod(&ANS->y,bufL,tmpL1);
    
    //ANS->z
    Lazy_mul(tmpL1,Pt1.z.x0,Pt2.z.x0);
    Lazy_mod(tmp1,tmpL1,FPLIMB2);
    Lazy_mul(bufL,tmp1,H);
    mpn_mod(&ANS->z,bufL,FPLIMB2);
}
void EFp_ECA_Jacobian_table(EFpZ *ANS,EFpZ *P1,EFpZT *P2){
    static EFpZ Pt1;
    static EFpZT Pt2;
    static mp_limb_t U1[FPLIMB],U2[FPLIMB],S1[FPLIMB],S2[FPLIMB],H[FPLIMB],r[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpZ1[FPLIMB],tmpZ2[FPLIMB],tmpH2[FPLIMB],tmpH3[FPLIMB],tmpU1H2[FPLIMB];
    Fp out;
    if(P1->infinity==1){
        EFpZT_to_EFpZ(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpZ_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD_Jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    EFpZ_set(&Pt1,P1);
    EFpZT_set(&Pt2,P2);

    //U1
    Lazy_mul(tmpL1,Pt2.zz.x0,Pt1.x.x0);
    Lazy_mod(U1,tmpL1,FPLIMB2);

    //U2
    //Lazy_sqr(bufL,Pt1.z.x0);
    Lazy_mul(bufL,Pt1.z.x0,Pt1.z.x0);
    Lazy_mod(tmpZ1,bufL,FPLIMB2);
    Lazy_mul(tmpL1,tmpZ1,Pt2.x.x0);
    Lazy_mod(U2,tmpL1,FPLIMB2);

    //S1
    Lazy_mul(tmpL1,Pt2.zzz.x0,Pt1.y.x0);
    Lazy_mod(S1,tmpL1,FPLIMB2);

    //S2
    Lazy_mul(bufL,tmpZ1,Pt1.z.x0);
    Lazy_mod(tmp1,bufL,FPLIMB2);
    Lazy_mul(tmpL1,tmp1,Pt2.y.x0);
    Lazy_mod(S2,tmpL1,FPLIMB2);
    
    //H
    Lazy_sub(buf,FPLIMB,U2,FPLIMB,U1,FPLIMB);
    Lazy_mod(H,buf,FPLIMB);

    //r
    Lazy_sub(buf,FPLIMB,S2,FPLIMB,S1,FPLIMB);
    Lazy_mod(r,buf,FPLIMB);

    //ANS->x
    //Lazy_sqr(bufL,r);
    Lazy_mul(bufL,r,r);

    //Lazy_sqr(tmpL1,H);
    Lazy_mul(tmpL1,H,H);
    Lazy_mod(tmpH2,tmpL1,FPLIMB2);
    Lazy_mul(tmpL1,tmpH2,H);
    Lazy_mod(tmpH3,tmpL1,FPLIMB2);
    Lazy_sub(bufL,FPLIMB2,bufL,FPLIMB2,tmpL1,FPLIMB2);

    Lazy_mul(tmpL1,tmpH2,U1);
    Lazy_mod(tmpU1H2,tmpL1,FPLIMB2);
    Lazy_add(tmpL1,FPLIMB2,tmpL1,FPLIMB2,tmpL1,FPLIMB2);
    Lazy_sub_mod(&ANS->x,bufL,tmpL1);

    //ANS->y
    Lazy_sub(tmp1,FPLIMB,tmpU1H2,FPLIMB,ANS->x.x0,FPLIMB);
    Lazy_mul(bufL,tmp1,r);

    Lazy_mul(tmpL1,tmpH3,S1);
    Lazy_sub_mod(&ANS->y,bufL,tmpL1);
    
    //ANS->z
    Lazy_mul(tmpL1,Pt1.z.x0,Pt2.z.x0);
    Lazy_mod(tmp1,tmpL1,FPLIMB2);
    Lazy_mul(bufL,tmp1,H);
    mpn_mod(&ANS->z,bufL,FPLIMB2);

}
void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFp Tmp_P,Next_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    EFp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFp_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    EFp_set(ANS,&Next_P);

}
void EFp_SCM_Jacobian(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFpZ Tmp_P,Next_P;
    EFp_to_EFpZ(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFpZ_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp_ECD_Jacobian(&Next_P,&Next_P);
        if(binary[i]=='1'){
	    EFp_ECA_Jacobian(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp_Jacobian(ANS,&Next_P);
}
void EFp_SCM_lazy(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFp Tmp_P,Next_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    EFp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFp_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp_ECD_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp_ECA_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    EFp_set(ANS,&Next_P);

}
void EFp_SCM_Jacobian_lazy(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFpZ Tmp_P,Next_P;
    EFp_to_EFpZ(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFpZ_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        EFp_ECD_Jacobian_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
	    EFp_ECA_Jacobian_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp_Jacobian(ANS,&Next_P);
}
//skew frobenius map
void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A){
    Fp_mul_mpn(&ANS->x,&A->x,epsilon1);
    Fp_set_neg(&ANS->y,&A->y);
}

