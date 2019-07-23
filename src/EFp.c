#include <ELiPS/EFp.h>
//EFp
void EFp_init(EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->infinity=0;
}
void EFpP_init(EFpP *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    Fp_init(&P->z);
    P->infinity=0;
}
void EFpJ_init(EFpJ *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    Fp_init(&P->z);
    P->infinity=0;
}
void EFpJT_init(EFpJT *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    Fp_init(&P->z);
    Fp_init(&P->zz);
    Fp_init(&P->zzz);
    P->infinity=0;
}
void EFp_printf(char *str,EFp *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf("",&P->x);
        printf(",");
        Fp_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}
void EFp_println(char *str,EFp *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf("",&P->x);
        printf(",");
        Fp_printf("",&P->y);
        printf(")\n");
    }else{
        printf("0\n");
    }
}
void EFpP_printf(char *str,EFpP *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf("",&P->x);
        printf(",");
        Fp_printf("",&P->y);
        printf(",");
        Fp_printf("",&P->z);
        printf(")");
    }else{
        printf("Infinity");
    }
}
void EFpJ_printf(char *str,EFpJ *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf("",&P->x);
        printf(",");
        Fp_printf("",&P->y);
        printf(",");
        Fp_printf("",&P->z);
        printf(")");
    }else{
        printf("Infinity");
    }
}
void EFpJT_printf(char *str,EFpJT *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf("",&P->x);
        printf(",");
        Fp_printf("",&P->y);
        printf(",");
        Fp_printf("",&P->z);
        printf(",");
        Fp_printf("",&P->zz);
        printf(",");
        Fp_printf("",&P->zzz);
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
void EFpP_set(EFpP *ANS,EFpP *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpJ_set(EFpJ *ANS,EFpJ *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpJT_set(EFpJT *ANS,EFpJT *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    Fp_set(&ANS->zz,&A->zz);
    Fp_set(&ANS->zzz,&A->zzz);
    ANS->infinity=A->infinity;
}
void EFp_to_EFpP(EFpP *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}
void EFp_to_EFpJ(EFpJ *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}
void EFpJT_to_EFpJ(EFpJ *ANS,EFpJT *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpJ_to_EFpJT(EFpJT *ANS,EFpJ *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFp_Projective(EFp *ANS,EFpP *A){
    static Fp Zi,Zt;
    Fp_inv(&Zi,&A->z);
    Fp_mul(&ANS->x,&A->x,&Zi);
    Fp_mul(&ANS->y,&A->y,&Zi);
    ANS->infinity=A->infinity;
}
void EFp_Jacobian(EFp *ANS,EFpJ *A){
    static Fp Zi,Zt;
    Fp_inv(&Zi,&A->z);
    Fp_mul(&Zt,&Zi,&Zi);
    Fp_mul(&ANS->x,&A->x,&Zt);
    Fp_mul(&Zt,&Zt,&Zi);
    Fp_mul(&ANS->y,&A->y,&Zt);
    ANS->infinity=A->infinity;
}
void EFp_mix(EFpJ *ANS,EFpJ *A,Fp *Zi){
    static Fp Zt;
    Fp_mul(&Zt,Zi,Zi);
    Fp_mul(&ANS->x,&A->x,&Zt);
    Fp_mul(&Zt,&Zt,Zi);
    Fp_mul(&ANS->y,&A->y,&Zt);
    Fp_set_ui(&ANS->z,1);
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
void EFpP_set_neg(EFpP *ANS,EFpP *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set_neg(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void EFpJ_set_neg(EFpJ *ANS,EFpJ *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set_neg(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
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
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
	
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

void EFp_ECD_Jacobian(EFpJ *ANS,EFpJ *P){
    static EFpJ Pt;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    static Fp s,m,T;
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFpJ_set(&Pt,P);
    
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
    
    Fp_add_lazy(tmp1,FPLIMB,P1t.y.x0,FPLIMB,P1t.y.x0,FPLIMB);
    mpn_invert(tmp1,tmp1,prime);

    Fp_mul_lazy(tmpL2,P1t.x.x0,P1t.x.x0);
    Fp_add_lazy(tmpL3,FPLIMB2,tmpL2,FPLIMB2,tmpL2,FPLIMB2);
    Fp_add_lazy(tmpL2,FPLIMB2,tmpL3,FPLIMB2,tmpL2,FPLIMB2);//3x^2
    Fp_mod(&tmpF2,tmpL2,FPLIMB2);

    Fp_mul_lazy(tmpL3,tmp1,tmpF2.x0);
    Fp_mod(&tmpF3,tmpL3,FPLIMB2);
    Fp_mul_lazy(tmpL2,tmpF3.x0,tmpF3.x0);

    Fp_add_lazy(tmp4,FPLIMB,P1t.x.x0,FPLIMB,P1t.x.x0,FPLIMB);
    Fp_sub_lazy(bufL,FPLIMB2,tmpL2,FPLIMB2,tmp4,FPLIMB);
    Fp_mod(&ANS->x,bufL,FPLIMB2);

    Fp_sub_lazy(tmp1,FPLIMB,P1t.x.x0,FPLIMB,ANS->x.x0,FPLIMB);
    Fp_mul_lazy(tmpL2,tmpF3.x0,tmp1);
    Fp_sub_lazy(bufL,FPLIMB2,tmpL2,FPLIMB2,P1t.y.x0,FPLIMB);
    Fp_mod(&ANS->y,bufL,FPLIMB2);

}
void EFp_ECD_Projective_lazy(EFpP *ANS,EFpP *P){
    static EFpP Pt;
    static mp_limb_t XX[FPLIMB],ZZ[FPLIMB];
    static mp_limb_t w[FPLIMB],s[FPLIMB],ss[FPLIMB],sss[FPLIMB];
    static mp_limb_t R[FPLIMB],RR[FPLIMB],B[FPLIMB],h[FPLIMB];
    
    static mp_limb_t bufL[FPLIMB2],tmpL[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp[FPLIMB];
    
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFpP_set(&Pt,P);
    mpn_zero(bufL,FPLIMB2);
    mpn_zero(tmpL,FPLIMB2);
    

    //XX
    Fp_sqr_lazy(bufL,Pt.x.x0);
    mpn_mod(XX,bufL,FPLIMB2);

    //ZZ
    Fp_sqr_lazy(bufL,Pt.z.x0);
    mpn_mod(ZZ,bufL,FPLIMB2);
    
    //w
    Fp_add_lazy_mod(tmp,FPLIMB,XX,FPLIMB,XX,FPLIMB);
    Fp_add_lazy_mod(w,FPLIMB,tmp,FPLIMB,XX,FPLIMB);
    
    //s
    Fp_mul_lazy(tmpL,Pt.y.x0,Pt.z.x0);//Y2*Z1
    Fp_add_lazy_mod(bufL,FPLIMB2,tmpL,FPLIMB2,tmpL,FPLIMB2);
    mpn_mod(s,bufL,FPLIMB2);
    
    //ss
    Fp_sqr_lazy(bufL,s);
    mpn_mod(ss,bufL,FPLIMB2);
    
    //sss
    Fp_mul_lazy(bufL,s,ss);
    mpn_mod(sss,bufL,FPLIMB2);
    
    //R
    Fp_mul_lazy(bufL,Pt.y.x0,s);//Y2*Z1
    mpn_mod(R,bufL,FPLIMB2);
    
    //RR
    Fp_sqr_lazy(bufL,R);
    mpn_mod(RR,bufL,FPLIMB2);
    
    //B
    Fp_add_lazy(tmp,FPLIMB,Pt.x.x0,FPLIMB,R,FPLIMB);
    Fp_sqr_lazy(tmpL,tmp);
    Fp_sub_lazy(tmpL,FPLIMB2,tmpL,FPLIMB2,XX,FPLIMB);
    Fp_sub_lazy(bufL,FPLIMB2,tmpL,FPLIMB2,RR,FPLIMB);
    mpn_mod(B,bufL,FPLIMB2);

    //h
    Fp_sqr_lazy(tmpL,w);
    Fp_sub_lazy(tmpL,FPLIMB2,tmpL,FPLIMB2,B,FPLIMB);
    Fp_sub_lazy(bufL,FPLIMB2,tmpL,FPLIMB2,B,FPLIMB);
    mpn_mod(h,bufL,FPLIMB2);
    
    //ANS->x
    Fp_mul_lazy(bufL,h,s);
    Fp_mod(&ANS->x,bufL,FPLIMB2);

    //ANS->y
    Fp_sub_lazy(tmp,FPLIMB,B,FPLIMB,h,FPLIMB);
    Fp_mul_lazy(bufL,tmp,w);
    mpn_mod(tmp,bufL,FPLIMB2);
    Fp_sub_lazy(tmp,FPLIMB,tmp,FPLIMB,RR,FPLIMB);
    Fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,RR,FPLIMB);
    Fp_mod(&ANS->y,buf,FPLIMB);
    
    //ANS->z
    mpn_copyd(ANS->z.x0,sss,FPLIMB);
}
void EFp_ECD_Jacobian_lazy(EFpJ *ANS,EFpJ *P){
    static mp_limb_t s[FPLIMB],m[FPLIMB],T[FPLIMB];

    static Fp bufF;
    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpy[FPLIMB];
    static EFpJ Pt;
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFpJ_set(&Pt,P);
    
    //s
    Fp_sqr_lazy(bufL,Pt.y.x0);
    mpn_mod(tmpy,bufL,FPLIMB2);

    Fp_mul_lazy(bufL,tmpy,Pt.x.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    
    Fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Fp_add_lazy(s,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);


    //m
    Fp_add_lazy(tmp1,FPLIMB,Pt.x.x0,FPLIMB,Pt.x.x0,FPLIMB);
    Fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,Pt.x.x0,FPLIMB);
    Fp_mul_lazy(bufL,tmp1,Pt.x.x0);
    mpn_mod(m,bufL,FPLIMB2);

    //T
    Fp_sqr_lazy(bufL,m);
    mpn_mod(T,bufL,FPLIMB2);//heraseru
    Fp_add_lazy(tmp1,FPLIMB,s,FPLIMB,s,FPLIMB);
    Fp_sub_lazy(buf,FPLIMB,T,FPLIMB,tmp1,FPLIMB);
    mpn_mod(T,buf,FPLIMB);

    //ANS->x
    mpn_copyd(ANS->x.x0,T,FPLIMB);

    //ANS->y
    Fp_sub_lazy(tmp1,FPLIMB,s,FPLIMB,T,FPLIMB);
    Fp_mul_lazy(bufL,tmp1,m);

    Fp_sqr_lazy(tmpL1,tmpy);
    mpn_mod(tmp1,tmpL1,FPLIMB2);
    Fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    Fp_sub_lazy(bufL,FPLIMB,bufL,FPLIMB2,tmp1,FPLIMB);

    Fp_mod(&ANS->y,bufL,FPLIMB2);

    
    //ANS->z
    Fp_add_lazy(tmp1,FPLIMB,Pt.y.x0,FPLIMB,Pt.y.x0,FPLIMB);
    Fp_mul_lazy(bufL,tmp1,Pt.z.x0);
    Fp_mod(&ANS->z,bufL,FPLIMB2);
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
void EFp_ECA_Jacobian(EFpJ *ANS,EFpJ *P1,EFpJ *P2){
    static EFpJ Pt1,Pt2;
    static Fp tmp1_Fp,tmp2_Fp,tmp3_Fp;
    static Fp U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        EFpJ_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpJ_set(ANS,P1);
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
    
    EFpJ_set(&Pt1,P1);
    EFpJ_set(&Pt2,P2);
    
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
void EFp_ECA_Mixture(EFpJ *ANS,EFpJ *P1,EFpJ *P2){
    static EFpJ Pt1,Pt2;
    static Fp Z1Z1,HH,I,J,V;
    static Fp U1,U2,S1,S2,H,r;

    static Fp bufL,tmpL1,tmpL2;
    static Fp buf,tmp1,tmp2;
    Fp out;
    if(P1->infinity==1){
        EFpJ_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpJ_set(ANS,P1);
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
    
    EFpJ_set(&Pt1,P1);
    EFpJ_set(&Pt2,P2);
    
    //Z1Z1
    Fp_mul(&Z1Z1,&Pt1.z,&Pt1.z);
    
    //U2
    Fp_mul(&U2,&Pt2.x,&Z1Z1);
    
    //S2
    Fp_mul(&tmp1,&Z1Z1,&Pt1.z);
    Fp_mul(&S2,&tmp1,&Pt2.y);
    
    //H
    Fp_sub(&H,&U2,&Pt1.x);
    
    //HH
    Fp_mul(&HH,&H,&H);
    
    //I
    Fp_add(&I,&HH,&HH);
    Fp_add(&I,&I,&I);
    
    //J
    Fp_mul(&J,&HH,&H);
    
    //r
    Fp_sub(&r,&S2,&Pt1.y);
    
    //V
    Fp_mul(&V,&Pt1.x,&HH);
    
    //X3
    Fp_mul(&tmp1,&r,&r);
    Fp_add(&tmp2,&V,&V);
    Fp_sub(&tmp1,&tmp1,&J);
    Fp_sub(&ANS->x,&tmp1,&tmp2);
    
    //Y3
    Fp_sub(&tmp1,&V,&ANS->x);
    Fp_mul(&tmp1,&tmp1,&r);
    Fp_mul(&tmp2,&Pt1.y,&J);
    Fp_sub(&ANS->y,&tmp1,&tmp2);
    
    //Z3
    Fp_mul(&ANS->z,&Pt1.z,&H);
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

    Fp_sub_lazy(tmp1,FPLIMB,P2t.x.x0,FPLIMB,P1t.x.x0,FPLIMB);
    mpn_invert(tmp1,tmp1,prime);
    Fp_sub_lazy(tmp2,FPLIMB,P2t.y.x0,FPLIMB,P1t.y.x0,FPLIMB);
    Fp_mul_lazy(tmpL3,tmp1,tmp2);
    Fp_mod(&tmpF4,tmpL3,FPLIMB2);
    Fp_mul_lazy(tmpL5,tmpF4.x0,tmpF4.x0);


    Fp_sub_lazy(tmpL6,FPLIMB2,tmpL5,FPLIMB2,P1t.x.x0,FPLIMB);
    Fp_sub_lazy(ANS->x.x0,FPLIMB2,tmpL6,FPLIMB2,P2t.x.x0,FPLIMB);
    Fp_mod(&ANS->x,ANS->x.x0,FPLIMB2);

    Fp_sub_lazy(tmp1,FPLIMB,P1t.x.x0,FPLIMB,ANS->x.x0,FPLIMB);
    Fp_mul_lazy(tmpL3,tmpF4.x0,tmp1);
    Fp_sub_lazy(bufL,FPLIMB2,tmpL3,FPLIMB2,P1t.y.x0,FPLIMB);
    Fp_mod(&ANS->y,bufL,FPLIMB2);

}
void EFp_ECA_Projective_lazy(EFpP *ANS,EFpP *P1,EFpP *P2){
    static EFpP Pt1,Pt2;
    static mp_limb_t Y1Z2[FPLIMB],X1Z2[FPLIMB],Z1Z2[FPLIMB];
    static mp_limb_t u[FPLIMB],uu[FPLIMB];
    static mp_limb_t v[FPLIMB],vv[FPLIMB],vvv[FPLIMB],R[FPLIMB],A[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp[FPLIMB];
    
    Fp out;
    if(P1->infinity==1){
        EFpP_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpP_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD_Projective_lazy(ANS,P1);
            return;
        }
    }
    mpn_zero(bufL,FPLIMB2);
    mpn_zero(tmpL,FPLIMB2);
    
    EFpP_set(&Pt1,P1);
    EFpP_set(&Pt2,P2);

    //Y1Z2
    Fp_mul_lazy(bufL,Pt1.y.x0,Pt2.z.x0);
    mpn_mod(Y1Z2,bufL,FPLIMB2);

    //X1Z2
    Fp_mul_lazy(bufL,Pt1.x.x0,Pt2.z.x0);
    mpn_mod(X1Z2,bufL,FPLIMB2);
    
    //Z1Z2
    Fp_mul_lazy(bufL,Pt1.z.x0,Pt2.z.x0);
    mpn_mod(Z1Z2,bufL,FPLIMB2);
    
    //u
    Fp_mul_lazy(tmpL,Pt1.z.x0,Pt2.y.x0);//Y2*Z1
    mpn_mod(tmp,tmpL,FPLIMB2);
    Fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,Y1Z2,FPLIMB);
    mpn_mod(u,buf,FPLIMB);
    
    //uu
    Fp_sqr_lazy(bufL,u);
    mpn_mod(uu,bufL,FPLIMB2);
    
    //v
    Fp_mul_lazy(tmpL,Pt1.z.x0,Pt2.x.x0);//X2*Z1
    mpn_mod(tmp,tmpL,FPLIMB2);
    Fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,X1Z2,FPLIMB);
    mpn_mod(v,buf,FPLIMB);
    
    //vv
    Fp_sqr_lazy(bufL,v);
    mpn_mod(vv,bufL,FPLIMB2);
    
    //vv
    Fp_mul_lazy(bufL,v,vv);
    mpn_mod(vvv,bufL,FPLIMB2);
    
    //R
    Fp_mul_lazy(bufL,X1Z2,vv);
    mpn_mod(R,bufL,FPLIMB2);

    //A
    Fp_mul_lazy(tmpL,Z1Z2,uu);
    mpn_mod(tmp,tmpL,FPLIMB2);
    Fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,vvv,FPLIMB);
    Fp_sub_lazy(tmp,FPLIMB,buf,FPLIMB,R,FPLIMB);
    Fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,R,FPLIMB);
    mpn_mod(A,buf,FPLIMB);
    //ANS->x
    Fp_mul_lazy(bufL,A,v);
    Fp_mod(&ANS->x,bufL,FPLIMB2);

    //ANS->y
    Fp_sub_lazy(tmp,FPLIMB,R,FPLIMB,A,FPLIMB);
    Fp_mul_lazy(tmpL,tmp,u);
    Fp_mul_lazy(bufL,vvv,Y1Z2);

    Fp_sub_lazy_mod(&ANS->y,tmpL,bufL);
    
    //ANS->z
    Fp_mul_lazy(bufL,Z1Z2,vvv);
    Fp_mod(&ANS->z,bufL,FPLIMB2);
}
void EFp_ECA_Jacobian_lazy(EFpJ *ANS,EFpJ *P1,EFpJ *P2){
    static EFpJ Pt1,Pt2;
    static mp_limb_t U1[FPLIMB],U2[FPLIMB],S1[FPLIMB],S2[FPLIMB],H[FPLIMB],r[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2],tmpL2[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpZ1[FPLIMB],tmpZ2[FPLIMB],tmpH2[FPLIMB],tmpH3[FPLIMB],tmpU1H2[FPLIMB];
    Fp out;
    if(P1->infinity==1){
        EFpJ_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpJ_set(ANS,P1);
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
    
    EFpJ_set(&Pt1,P1);
    EFpJ_set(&Pt2,P2);

    //U1
    Fp_sqr_lazy(bufL,Pt2.z.x0);
    //Fp_mul_lazy(bufL,Pt2.z.x0,Pt2.z.x0);
    mpn_mod(tmpZ2,bufL,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmpZ2,Pt1.x.x0);
    mpn_mod(U1,tmpL1,FPLIMB2);

    //U2
    Fp_sqr_lazy(bufL,Pt1.z.x0);
    //Fp_mul_lazy(bufL,Pt1.z.x0,Pt1.z.x0);
    mpn_mod(tmpZ1,bufL,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmpZ1,Pt2.x.x0);
    mpn_mod(U2,tmpL1,FPLIMB2);

    //S1
    Fp_mul_lazy(bufL,tmpZ2,Pt2.z.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmp1,Pt1.y.x0);
    mpn_mod(S1,tmpL1,FPLIMB2);
    //gmp_printf("S1=%Nu\n",S1,FPLIMB);

    //S2
    Fp_mul_lazy(bufL,tmpZ1,Pt1.z.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmp1,Pt2.y.x0);
    mpn_mod(S2,tmpL1,FPLIMB2);
    //gmp_printf("S2=%Nu\n",S2,FPLIMB);
    
    //H
    Fp_sub_lazy(buf,FPLIMB,U2,FPLIMB,U1,FPLIMB);
    mpn_mod(H,buf,FPLIMB);

    //r
    Fp_sub_lazy(buf,FPLIMB,S2,FPLIMB,S1,FPLIMB);
    mpn_mod(r,buf,FPLIMB);
    //gmp_printf("r=%Nu\n",r,FPLIMB);

    //ANS->x
    Fp_sqr_lazy(tmpL2,r);
    //Fp_mul_lazy(tmpL2,r,r);

    Fp_sqr_lazy(tmpL1,H);
    //Fp_mul_lazy(tmpL1,H,H);
    mpn_mod(tmpH2,tmpL1,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmpH2,H);
    mpn_mod(tmpH3,tmpL1,FPLIMB2);
    Fp_sub_lazy(bufL,FPLIMB2,tmpL2,FPLIMB2,tmpL1,FPLIMB2);

    Fp_mul_lazy(tmpL1,tmpH2,U1);
    mpn_mod(tmpU1H2,tmpL1,FPLIMB2);
    Fp_add_lazy(tmpL1,FPLIMB2,tmpL1,FPLIMB2,tmpL1,FPLIMB2);
    Fp_sub_lazy_mod(&ANS->x,bufL,tmpL1);

    //ANS->y
    Fp_sub_lazy(tmp1,FPLIMB,tmpU1H2,FPLIMB,ANS->x.x0,FPLIMB);
    Fp_mul_lazy(bufL,tmp1,r);

    Fp_mul_lazy(tmpL1,tmpH3,S1);
    Fp_sub_lazy_mod(&ANS->y,bufL,tmpL1);
    
    //ANS->z
    Fp_mul_lazy(tmpL1,Pt1.z.x0,Pt2.z.x0);
    mpn_mod(tmp1,tmpL1,FPLIMB2);
    Fp_mul_lazy(bufL,tmp1,H);
    Fp_mod(&ANS->z,bufL,FPLIMB2);
}

void EFp_ECA_Mixture_lazy(EFpJ *ANS,EFpJ *P1,EFpJ *P2){
    static EFpJ Pt1,Pt2;
    static mp_limb_t Z1Z1[FPLIMB],HH[FPLIMB],I[FPLIMB],J[FPLIMB],V[FPLIMB];
    static mp_limb_t U1[FPLIMB],U2[FPLIMB],S1[FPLIMB],S2[FPLIMB],H[FPLIMB],r[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2],tmpL2[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB],tmp2[FPLIMB];
    Fp out;
    if(P1->infinity==1){
        EFpJ_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpJ_set(ANS,P1);
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
    
    EFpJ_set(&Pt1,P1);
    EFpJ_set(&Pt2,P2);
    /*
    Fp_printf("\nX1=",&Pt1.x);printf("\n");
    Fp_printf("Y1=",&Pt1.y);printf("\n");
    Fp_printf("Z1=",&Pt1.z);printf("\n");
    Fp_printf("X2=",&Pt2.x);printf("\n");
    Fp_printf("Y2=",&Pt2.y);printf("\n");
    Fp_printf("Z2=",&Pt2.z);printf("\n");
*/
    //Z1Z1
    Fp_mul_lazy(bufL,Pt1.z.x0,Pt1.z.x0);
    mpn_mod(Z1Z1,bufL,FPLIMB2);
    //gmp_printf("Z1Z1=%Nu\n",Z1Z1,FPLIMB);
    
    //U2
    Fp_mul_lazy(bufL,Pt2.x.x0,Z1Z1);
    mpn_mod(U2,bufL,FPLIMB2);
    //gmp_printf("U2=%Nu\n",U2,FPLIMB);
    
    //S2
    Fp_mul_lazy(bufL,Z1Z1,Pt1.z.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmp1,Pt2.y.x0);
    mpn_mod(S2,tmpL1,FPLIMB2);
    //gmp_printf("S2=%Nu\n",S2,FPLIMB);
    
    //H
    Fp_sub_lazy(buf,FPLIMB,U2,FPLIMB,Pt1.x.x0,FPLIMB);
    mpn_mod(H,buf,FPLIMB);
    //gmp_printf("H=%Nu\n",H,FPLIMB);
    
    //HH
    Fp_mul_lazy(tmpL1,H,H);
    mpn_mod(HH,tmpL1,FPLIMB2);
    //gmp_printf("HH=%Nu\n",HH,FPLIMB);
    
    //I
    Fp_add_lazy(I,FPLIMB,HH,FPLIMB,HH,FPLIMB);
    Fp_add_lazy(I,FPLIMB,I,FPLIMB,I,FPLIMB);
    mpn_mod(I,I,FPLIMB);
    //gmp_printf("I=%Nu\n",I,FPLIMB);
    
    //J
    //Fp_mul_lazy(bufL,I,H);
    Fp_mul_lazy(bufL,HH,H);
    mpn_mod(J,bufL,FPLIMB2);
    //gmp_printf("J=%Nu\n",J,FPLIMB);
    
    //r
    //gmp_printf("Y1=%Nu\n",Pt1.y.x0,FPLIMB);
    Fp_sub_lazy(buf,FPLIMB,S2,FPLIMB,Pt1.y.x0,FPLIMB);
    //Fp_add_lazy(buf,FPLIMB,buf,FPLIMB,buf,FPLIMB);
    mpn_mod(r,buf,FPLIMB);
    //gmp_printf("r=%Nu\n",r,FPLIMB);
    
    //V
    //Fp_mul_lazy(bufL,Pt1.x.x0,I);
    Fp_mul_lazy(bufL,Pt1.x.x0,HH);
    mpn_mod(V,bufL,FPLIMB2);
    //gmp_printf("V=%Nu\n",V,FPLIMB);
    
    //X3
    Fp_mul_lazy(tmpL1,r,r);
    Fp_add_lazy(tmp1,FPLIMB,V,FPLIMB,V,FPLIMB);
    Fp_sub_lazy(bufL,FPLIMB2,tmpL1,FPLIMB2,J,FPLIMB);
    Fp_sub_lazy(bufL,FPLIMB2,bufL,FPLIMB2,tmp1,FPLIMB);
    Fp_mod(&ANS->x,bufL,FPLIMB2);
    
    //Y3
    Fp_sub_lazy(tmp1,FPLIMB,V,FPLIMB,ANS->x.x0,FPLIMB);
    Fp_mul_lazy(tmpL1,tmp1,r);
    Fp_mul_lazy(tmpL2,Pt1.y.x0,J);
    //Fp_add_lazy(tmpL2,FPLIMB2,tmpL2,FPLIMB2,tmpL2,FPLIMB2);
    Fp_sub_lazy_mod(&ANS->y,tmpL1,tmpL2);
    
    //Z3
    /*
    Fp_add_lazy(tmp1,FPLIMB,Pt1.z.x0,FPLIMB,H,FPLIMB);
    Fp_mul_lazy(tmpL1,tmp1,tmp1);
    mpn_mod(tmp2,tmpL1,FPLIMB2);
    Fp_sub_lazy(buf,FPLIMB,tmp2,FPLIMB,Z1Z1,FPLIMB);
    Fp_sub_lazy(buf,FPLIMB,buf,FPLIMB,HH,FPLIMB);
    Fp_mod(&ANS->z,buf,FPLIMB);
    */
    
    //ANS->z
    Fp_mul_lazy(bufL,Pt1.z.x0,H);
    Fp_mod(&ANS->z,bufL,FPLIMB2);

}
void EFp_ECA_Jacobian_table(EFpJ *ANS,EFpJ *P1,EFpJT *P2){
    static EFpJ Pt1;
    static EFpJT Pt2;
    static mp_limb_t U1[FPLIMB],U2[FPLIMB],S1[FPLIMB],S2[FPLIMB],H[FPLIMB],r[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpZ1[FPLIMB],tmpZ2[FPLIMB],tmpH2[FPLIMB],tmpH3[FPLIMB],tmpU1H2[FPLIMB];
    Fp out;
    if(P1->infinity==1){
        EFpJT_to_EFpJ(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFpJ_set(ANS,P1);
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
    
    EFpJ_set(&Pt1,P1);
    EFpJT_set(&Pt2,P2);

    //U1
    Fp_mul_lazy(tmpL1,Pt2.zz.x0,Pt1.x.x0);
    mpn_mod(U1,tmpL1,FPLIMB2);

    //U2
    //Fp_sqr_lazy(bufL,Pt1.z.x0);
    Fp_mul_lazy(bufL,Pt1.z.x0,Pt1.z.x0);
    mpn_mod(tmpZ1,bufL,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmpZ1,Pt2.x.x0);
    mpn_mod(U2,tmpL1,FPLIMB2);

    //S1
    Fp_mul_lazy(tmpL1,Pt2.zzz.x0,Pt1.y.x0);
    mpn_mod(S1,tmpL1,FPLIMB2);

    //S2
    Fp_mul_lazy(bufL,tmpZ1,Pt1.z.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmp1,Pt2.y.x0);
    mpn_mod(S2,tmpL1,FPLIMB2);
    
    //H
    Fp_sub_lazy(buf,FPLIMB,U2,FPLIMB,U1,FPLIMB);
    mpn_mod(H,buf,FPLIMB);

    //r
    Fp_sub_lazy(buf,FPLIMB,S2,FPLIMB,S1,FPLIMB);
    mpn_mod(r,buf,FPLIMB);

    //ANS->x
    //Fp_sqr_lazy(bufL,r);
    Fp_mul_lazy(bufL,r,r);

    //Fp_sqr_lazy(tmpL1,H);
    Fp_mul_lazy(tmpL1,H,H);
    mpn_mod(tmpH2,tmpL1,FPLIMB2);
    Fp_mul_lazy(tmpL1,tmpH2,H);
    mpn_mod(tmpH3,tmpL1,FPLIMB2);
    Fp_sub_lazy(bufL,FPLIMB2,bufL,FPLIMB2,tmpL1,FPLIMB2);

    Fp_mul_lazy(tmpL1,tmpH2,U1);
    mpn_mod(tmpU1H2,tmpL1,FPLIMB2);
    Fp_add_lazy(tmpL1,FPLIMB2,tmpL1,FPLIMB2,tmpL1,FPLIMB2);
    Fp_sub_lazy_mod(&ANS->x,bufL,tmpL1);

    //ANS->y
    Fp_sub_lazy(tmp1,FPLIMB,tmpU1H2,FPLIMB,ANS->x.x0,FPLIMB);
    Fp_mul_lazy(bufL,tmp1,r);

    Fp_mul_lazy(tmpL1,tmpH3,S1);
    Fp_sub_lazy_mod(&ANS->y,bufL,tmpL1);
    
    //ANS->z
    Fp_mul_lazy(tmpL1,Pt1.z.x0,Pt2.z.x0);
    mpn_mod(tmp1,tmpL1,FPLIMB2);
    Fp_mul_lazy(bufL,tmp1,H);
    Fp_mod(&ANS->z,bufL,FPLIMB2);

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
    
    EFpJ Tmp_P,Next_P;
    EFp_to_EFpJ(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFpJ_set(&Next_P,&Tmp_P);
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
    
    EFpJ Tmp_P,Next_P;
    EFp_to_EFpJ(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    EFpJ_set(&Next_P,&Tmp_P);
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
void EFpJ_skew_frobenius_map_p2(EFpJ *ANS,EFpJ *A){
    Fp_mul_mpn(&ANS->x,&A->x,epsilon1);
    Fp_set_neg(&ANS->y,&A->y);
    Fp_set(&ANS->z,&A->z);
}

