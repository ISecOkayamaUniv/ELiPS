#include <ELiPS/efp.h>
//efp
void efp_init(efp_t *P){
    fp_init(&P->x);
    fp_init(&P->y);
    P->infinity=0;
}
void efp_projective_init(efp_projective_t *P){
    fp_init(&P->x);
    fp_init(&P->y);
    fp_init(&P->z);
    P->infinity=0;
}
void efp_jacobian_init(efp_jacobian_t *P){
    fp_init(&P->x);
    fp_init(&P->y);
    fp_init(&P->z);
    P->infinity=0;
}
void efp_printf(char *str,efp_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp_printf("",&P->x);
        printf(",");
        fp_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}
void efp_println(char *str,efp_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp_printf("",&P->x);
        printf(",");
        fp_printf("",&P->y);
        printf(")\n");
    }else{
        printf("0\n");
    }
}
void efp_projective_printf(char *str,efp_projective_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp_printf("",&P->x);
        printf(",");
        fp_printf("",&P->y);
        printf(",");
        fp_printf("",&P->z);
        printf(")");
    }else{
        printf("Infinity");
    }
}
void efp_jacobian_printf(char *str,efp_jacobian_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp_printf("",&P->x);
        printf(",");
        fp_printf("",&P->y);
        printf(",");
        fp_printf("",&P->z);
        printf(")");
    }else{
        printf("Infinity");
    }
}
void efp_set(efp_t *ANS,efp_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void efp_projective_set(efp_projective_t *ANS,efp_projective_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set(&ANS->y,&A->y);
    fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void efp_jacobian_set(efp_jacobian_t *ANS,efp_jacobian_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set(&ANS->y,&A->y);
    fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void efp_affine_to_projective(efp_projective_t *ANS,efp_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set(&ANS->y,&A->y);
    fp_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}
void efp_affine_to_jacobian(efp_jacobian_t *ANS,efp_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set(&ANS->y,&A->y);
    fp_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}
void efp_affine_to_jacobian_montgomery(efp_jacobian_t *ANS,efp_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set(&ANS->y,&A->y);
    mpn_copyd(ANS->z.x0,RmodP,FPLIMB);
    ANS->infinity=A->infinity;
}
void efp_projective_to_affine(efp_t *ANS,efp_projective_t *A){
    static fp_t Zi,Zt;
    //TODO:mul->mul_lazy
    fp_inv(&Zi,&A->z);
    fp_mul(&ANS->x,&A->x,&Zi);
    fp_mul(&ANS->y,&A->y,&Zi);
    ANS->infinity=A->infinity;
}
void efp_jacobian_to_affine(efp_t *ANS,efp_jacobian_t *A){
    static fp_t Zi,Zt;
    //TODO:mul->mul_lazy
    fp_inv(&Zi,&A->z);
    fp_mul(&Zt,&Zi,&Zi);
    fp_mul(&ANS->x,&A->x,&Zt);
    fp_mul(&Zt,&Zt,&Zi);
    fp_mul(&ANS->y,&A->y,&Zt);
    ANS->infinity=A->infinity;
}
void efp_jacobian_montgomery(efp_t *ANS,efp_jacobian_t *A){
    static fp_t Zi,Zt;
    fp_inv_montgomery(&Zi,&A->z);
    fp_mulmod_montgomery(&Zt,&Zi,&A->z);
    fp_mulmod_montgomery(&Zt,&Zi,&Zi);
    fp_mulmod_montgomery(&ANS->x,&A->x,&Zt);
    fp_mulmod_montgomery(&Zt,&Zt,&Zi);
    fp_mulmod_montgomery(&ANS->y,&A->y,&Zt);
    ANS->infinity=A->infinity;
}
void efp_mix(efp_jacobian_t *ANS,efp_jacobian_t *A,fp_t *Zi){
    static fp_t Zt;
    fp_mul(&Zt,Zi,Zi);
    fp_mul(&ANS->x,&A->x,&Zt);
    fp_mul(&Zt,&Zt,Zi);
    fp_mul(&ANS->y,&A->y,&Zt);
    fp_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}
void efp_mix_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *A,fp_t *Zi){
    static fp_t Zt;
    fp_mulmod_montgomery(&Zt,Zi,Zi);
    fp_mulmod_montgomery(&ANS->x,&A->x,&Zt);
    fp_mulmod_montgomery(&Zt,&Zt,Zi);
    fp_mulmod_montgomery(&ANS->y,&A->y,&Zt);
    //fp_set_ui(&ANS->z,1);
    mpn_copyd(ANS->z.x0,RmodP,FPLIMB);
    ANS->infinity=A->infinity;
}
void efp_to_montgomery(efp_t *ANS,efp_t *A){
    fp_to_montgomery(&ANS->x,&A->x);
    fp_to_montgomery(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void efp_mod_montgomery(efp_t *ANS,efp_t *A){
    fp_mod_montgomery(&ANS->x,&A->x);
    fp_mod_montgomery(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void efp_set_ui(efp_t *ANS,unsigned long int UI){
    fp_set_ui(&ANS->x,UI);
    fp_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}
void efp_set_mpn(efp_t *ANS,mp_limb_t *A){
    fp_set_mpn(&ANS->x,A);
    fp_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}
void efp_set_neg(efp_t *ANS,efp_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void efp_projective_set_neg(efp_projective_t *ANS,efp_projective_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set_neg(&ANS->y,&A->y);
    fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void efp_jacobian_set_neg(efp_jacobian_t *ANS,efp_jacobian_t *A){
    fp_set(&ANS->x,&A->x);
    fp_set_neg(&ANS->y,&A->y);
    fp_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}


int  efp_cmp(efp_t *A,efp_t *B){
    if(fp_cmp(&A->x,&B->x)==0 && fp_cmp(&A->y,&B->y)==0){
        return 0;   
    }else if(A->infinity==1&&B->infinity==1){
	return 0;
    }else{
    return 1;
    }
}
void efp_rational_point(efp_t *P){
    fp_t tmp1,tmp2,tmp_x;
    fp_init(&tmp1);
    fp_init(&tmp2);
    fp_init(&tmp_x);
	
    while(1){
        fp_set_random(&P->x,state);
        fp_mul(&tmp1,&P->x,&P->x);
        fp_mul(&tmp2,&tmp1,&P->x);
        fp_add_mpn(&tmp_x,&tmp2,curve_b);
        if(fp_legendre(&tmp_x)==1){
            fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
}

void efp_ecd(efp_t *ANS,efp_t *P){
    static efp_t tmp1_efp;
    static fp_t tmp1_fp,tmp2_fp,tmp3_fp;
    if(fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    efp_set(&tmp1_efp,P);
    
    fp_add(&tmp1_fp,&tmp1_efp.y,&tmp1_efp.y);
    fp_inv(&tmp1_fp,&tmp1_fp);


    fp_mul(&tmp2_fp,&tmp1_efp.x,&tmp1_efp.x);
    fp_add(&tmp3_fp,&tmp2_fp,&tmp2_fp);
    fp_add(&tmp2_fp,&tmp2_fp,&tmp3_fp);

    fp_mul(&tmp3_fp,&tmp1_fp,&tmp2_fp);
    fp_mul(&tmp1_fp,&tmp3_fp,&tmp3_fp);

    fp_add(&tmp2_fp,&tmp1_efp.x,&tmp1_efp.x);
    fp_sub(&ANS->x,&tmp1_fp,&tmp2_fp);

    fp_sub(&tmp1_fp,&tmp1_efp.x,&ANS->x);
    fp_mul(&tmp2_fp,&tmp3_fp,&tmp1_fp);
    fp_sub(&ANS->y,&tmp2_fp,&tmp1_efp.y);
}

void efp_ecd_jacobian(efp_jacobian_t *ANS,efp_jacobian_t *P){
    static efp_jacobian_t Pt;
    static fp_t tmp1_fp,tmp2_fp,tmp3_fp;
    static fp_t s,m,T;
    if(fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    efp_jacobian_set(&Pt,P);
    
    //s
    fp_mul_ui(&tmp1_fp,&Pt.x,4);
    fp_mul(&tmp1_fp,&tmp1_fp,&Pt.y);
    fp_mul(&s,&tmp1_fp,&Pt.y);

    //m
    fp_mul_ui(&tmp2_fp,&Pt.x,3);
    fp_mul(&m,&tmp2_fp,&Pt.x);

    //T
    fp_add(&tmp1_fp,&s,&s);
    fp_set_neg(&T,&tmp1_fp);
    fp_mul(&tmp1_fp,&m,&m);
    fp_add(&T,&T,&tmp1_fp);


    //ANS->x
    fp_set(&ANS->x,&T);

    //ANS->y
    fp_mul_ui(&tmp1_fp,&Pt.y,8);
    fp_mul(&tmp1_fp,&tmp1_fp,&Pt.y);
    fp_mul(&tmp1_fp,&tmp1_fp,&Pt.y);
    fp_mul(&tmp1_fp,&tmp1_fp,&Pt.y);
    fp_set_neg(&ANS->y,&tmp1_fp);

    fp_sub(&tmp1_fp,&s,&T);
    fp_mul(&tmp1_fp,&tmp1_fp,&m);
    fp_add(&ANS->y,&ANS->y,&tmp1_fp);
    
    //ANS->z
    fp_add(&tmp1_fp,&Pt.y,&Pt.y);
    fp_mul(&ANS->z,&tmp1_fp,&Pt.z);
}

void efp_ecd_lazy(efp_t *ANS,efp_t *P){
    static fp_t tmpF2,tmpF3;
    static mp_limb_t bufL[FPLIMB2],tmpL2[FPLIMB2],tmpL3[FPLIMB2],tmpL5[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB],tmp4[FPLIMB];
    static efp_t P1t;
    if(fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    efp_set(&P1t,P);
    
    fp_add_lazy(tmp1,FPLIMB,P1t.y.x0,FPLIMB,P1t.y.x0,FPLIMB);
    mpn_invert(tmp1,tmp1,prime);

    fp_mul_lazy(tmpL2,P1t.x.x0,P1t.x.x0);
    fp_add_lazy(tmpL3,FPLIMB2,tmpL2,FPLIMB2,tmpL2,FPLIMB2);
    fp_add_lazy(tmpL2,FPLIMB2,tmpL3,FPLIMB2,tmpL2,FPLIMB2);//3x^2
    fp_mod(&tmpF2,tmpL2,FPLIMB2);

    fp_mul_lazy(tmpL3,tmp1,tmpF2.x0);
    fp_mod(&tmpF3,tmpL3,FPLIMB2);
    fp_mul_lazy(tmpL2,tmpF3.x0,tmpF3.x0);

    fp_add_lazy(tmp4,FPLIMB,P1t.x.x0,FPLIMB,P1t.x.x0,FPLIMB);
    fp_sub_lazy(bufL,FPLIMB2,tmpL2,FPLIMB2,tmp4,FPLIMB);
    fp_mod(&ANS->x,bufL,FPLIMB2);

    fp_sub_lazy(tmp1,FPLIMB,P1t.x.x0,FPLIMB,ANS->x.x0,FPLIMB);
    fp_mul_lazy(tmpL2,tmpF3.x0,tmp1);
    fp_sub_lazy(bufL,FPLIMB2,tmpL2,FPLIMB2,P1t.y.x0,FPLIMB);
    fp_mod(&ANS->y,bufL,FPLIMB2);

}
void efp_ecd_projective_lazy(efp_projective_t *ANS,efp_projective_t *P){
    static efp_projective_t Pt;
    static mp_limb_t XX[FPLIMB],ZZ[FPLIMB];
    static mp_limb_t w[FPLIMB],s[FPLIMB],ss[FPLIMB],sss[FPLIMB];
    static mp_limb_t R[FPLIMB],RR[FPLIMB],B[FPLIMB],h[FPLIMB];
    
    static mp_limb_t bufL[FPLIMB2],tmpL[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp[FPLIMB];
    
    if(fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    efp_projective_set(&Pt,P);
    mpn_zero(bufL,FPLIMB2);
    mpn_zero(tmpL,FPLIMB2);
    

    //XX
    fp_sqr_lazy(bufL,Pt.x.x0);
    mpn_mod(XX,bufL,FPLIMB2);

    //ZZ
    fp_sqr_lazy(bufL,Pt.z.x0);
    mpn_mod(ZZ,bufL,FPLIMB2);
    
    //w
    fp_add_lazy_mod(tmp,FPLIMB,XX,FPLIMB,XX,FPLIMB);
    fp_add_lazy_mod(w,FPLIMB,tmp,FPLIMB,XX,FPLIMB);
    
    //s
    fp_mul_lazy(tmpL,Pt.y.x0,Pt.z.x0);//Y2*Z1
    fp_add_lazy_mod(bufL,FPLIMB2,tmpL,FPLIMB2,tmpL,FPLIMB2);
    mpn_mod(s,bufL,FPLIMB2);
    
    //ss
    fp_sqr_lazy(bufL,s);
    mpn_mod(ss,bufL,FPLIMB2);
    
    //sss
    fp_mul_lazy(bufL,s,ss);
    mpn_mod(sss,bufL,FPLIMB2);
    
    //R
    fp_mul_lazy(bufL,Pt.y.x0,s);//Y2*Z1
    mpn_mod(R,bufL,FPLIMB2);
    
    //RR
    fp_sqr_lazy(bufL,R);
    mpn_mod(RR,bufL,FPLIMB2);
    
    //B
    fp_add_lazy(tmp,FPLIMB,Pt.x.x0,FPLIMB,R,FPLIMB);
    fp_sqr_lazy(tmpL,tmp);
    fp_sub_lazy(tmpL,FPLIMB2,tmpL,FPLIMB2,XX,FPLIMB);
    fp_sub_lazy(bufL,FPLIMB2,tmpL,FPLIMB2,RR,FPLIMB);
    mpn_mod(B,bufL,FPLIMB2);

    //h
    fp_sqr_lazy(tmpL,w);
    fp_sub_lazy(tmpL,FPLIMB2,tmpL,FPLIMB2,B,FPLIMB);
    fp_sub_lazy(bufL,FPLIMB2,tmpL,FPLIMB2,B,FPLIMB);
    mpn_mod(h,bufL,FPLIMB2);
    
    //ANS->x
    fp_mul_lazy(bufL,h,s);
    fp_mod(&ANS->x,bufL,FPLIMB2);

    //ANS->y
    fp_sub_lazy(tmp,FPLIMB,B,FPLIMB,h,FPLIMB);
    fp_mul_lazy(bufL,tmp,w);
    mpn_mod(tmp,bufL,FPLIMB2);
    fp_sub_lazy(tmp,FPLIMB,tmp,FPLIMB,RR,FPLIMB);
    fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,RR,FPLIMB);
    fp_mod(&ANS->y,buf,FPLIMB);
    
    //ANS->z
    mpn_copyd(ANS->z.x0,sss,FPLIMB);
}
void efp_ecd_jacobian_lazy(efp_jacobian_t *ANS,efp_jacobian_t *P){
    static mp_limb_t s[FPLIMB],m[FPLIMB],T[FPLIMB];

    static fp_t bufF;
    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpy[FPLIMB];
    static efp_jacobian_t Pt;
    if(fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    efp_jacobian_set(&Pt,P);
    
    //s
    fp_sqr_lazy(bufL,Pt.y.x0);
    mpn_mod(tmpy,bufL,FPLIMB2);

    fp_mul_lazy(bufL,tmpy,Pt.x.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    
    fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    fp_add_lazy(s,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    

    //m
    fp_add_lazy(tmp1,FPLIMB,Pt.x.x0,FPLIMB,Pt.x.x0,FPLIMB);
    fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,Pt.x.x0,FPLIMB);
    fp_mul_lazy(bufL,tmp1,Pt.x.x0);
    mpn_mod(m,bufL,FPLIMB2);

    //T
    fp_sqr_lazy(bufL,m);
    mpn_mod(T,bufL,FPLIMB2);//heraseru
    fp_add_lazy(tmp1,FPLIMB,s,FPLIMB,s,FPLIMB);
    fp_sub_lazy(buf,FPLIMB,T,FPLIMB,tmp1,FPLIMB);
    mpn_mod(T,buf,FPLIMB);

    //ANS->x
    mpn_copyd(ANS->x.x0,T,FPLIMB);

    //ANS->y
    fp_sub_lazy(tmp1,FPLIMB,s,FPLIMB,T,FPLIMB);
    fp_mul_lazy(bufL,tmp1,m);

    fp_sqr_lazy(tmpL1,tmpy);
    mpn_mod(tmp1,tmpL1,FPLIMB2);
    fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    fp_add_lazy(tmp1,FPLIMB,tmp1,FPLIMB,tmp1,FPLIMB);
    fp_sub_lazy(bufL,FPLIMB,bufL,FPLIMB2,tmp1,FPLIMB);

    fp_mod(&ANS->y,bufL,FPLIMB2);

    
    //ANS->z
    fp_add_lazy(tmp1,FPLIMB,Pt.y.x0,FPLIMB,Pt.y.x0,FPLIMB);
    fp_mul_lazy(bufL,tmp1,Pt.z.x0);
    fp_mod(&ANS->z,bufL,FPLIMB2);
}
void efp_ecd_jacobian_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P){
    static fp_t s,m,T;

    static fp_t buf,tmp1;
    static fp_t tmpY2;
    static efp_jacobian_t Pt;
    if(fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    efp_jacobian_set(&Pt,P);
    
    //s
    fp_mulmod_montgomery(&tmpY2,&Pt.y,&Pt.y);

    fp_mulmod_montgomery(&tmp1,&tmpY2,&Pt.x);
    fp_add(&tmp1,&tmp1,&tmp1);
    fp_add(&s,&tmp1,&tmp1);
    //m
    fp_add_lazy(tmp1.x0,FPLIMB,Pt.x.x0,FPLIMB,Pt.x.x0,FPLIMB);
    fp_add_lazy(tmp1.x0,FPLIMB,tmp1.x0,FPLIMB,Pt.x.x0,FPLIMB);
    fp_mulmod_montgomery(&m,&tmp1,&Pt.x);

    //T
    fp_mulmod_montgomery(&T,&m,&m);
    fp_add(&tmp1,&s,&s);
    fp_sub(&T,&T,&tmp1);

    //ANS->x
    fp_set(&ANS->x,&T);

    //ANS->y
    fp_sub_lazy(tmp1.x0,FPLIMB,s.x0,FPLIMB,T.x0,FPLIMB);
    fp_mulmod_montgomery(&buf,&tmp1,&m);

    fp_mulmod_montgomery(&tmp1,&tmpY2,&tmpY2);
    fp_add(&tmp1,&tmp1,&tmp1);
    fp_add(&tmp1,&tmp1,&tmp1);
    fp_add(&tmp1,&tmp1,&tmp1);
    fp_sub(&ANS->y,&buf,&tmp1);
    
    //ANS->z
    fp_add_lazy(tmp1.x0,FPLIMB,Pt.y.x0,FPLIMB,Pt.y.x0,FPLIMB);
    fp_mulmod_montgomery(&ANS->z,&tmp1,&Pt.z);
}
void efp_eca(efp_t *ANS,efp_t *P1,efp_t *P2){
    static efp_t tmp1_efp,tmp2_efp;
    static fp_t tmp1_fp,tmp2_fp,tmp3_fp;
    if(P1->infinity==1){
        efp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd(ANS,P1);
            return;
        }
    }
    
    efp_set(&tmp1_efp,P1);
    efp_set(&tmp2_efp,P2);
    
    fp_sub(&tmp1_fp,&tmp2_efp.x,&tmp1_efp.x);
    fp_inv(&tmp1_fp,&tmp1_fp);
    fp_sub(&tmp2_fp,&tmp2_efp.y,&tmp1_efp.y);
    fp_mul(&tmp3_fp,&tmp1_fp,&tmp2_fp);
    fp_mul(&tmp1_fp,&tmp3_fp,&tmp3_fp);


    fp_sub(&tmp2_fp,&tmp1_fp,&tmp1_efp.x);
    fp_sub(&ANS->x,&tmp2_fp,&tmp2_efp.x);

    fp_sub(&tmp1_fp,&tmp1_efp.x,&ANS->x);
    fp_mul(&tmp2_fp,&tmp3_fp,&tmp1_fp);
    fp_sub(&ANS->y,&tmp2_fp,&tmp1_efp.y);
}
void efp_eca_jacobian(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2){
    static efp_jacobian_t Pt1,Pt2;
    static fp_t tmp1_fp,tmp2_fp,tmp3_fp;
    static fp_t U1,U2,S1,S2,H,r;
    if(P1->infinity==1){
        efp_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_jacobian_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0&&fp_cmp(&P1->z,&P2->z)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0&&fp_cmp(&P1->z,&P2->z)==0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_jacobian(ANS,P1);
            return;
        }
    }
    
    efp_jacobian_set(&Pt1,P1);
    efp_jacobian_set(&Pt2,P2);
    
    //U1
    fp_mul(&tmp1_fp,&Pt1.x,&Pt2.z);
    fp_mul(&U1,&tmp1_fp,&Pt2.z);

    //U2
    fp_mul(&tmp1_fp,&Pt2.x,&Pt1.z);
    fp_mul(&U2,&tmp1_fp,&Pt1.z);

    //S1
    fp_mul(&tmp1_fp,&Pt1.y,&Pt2.z);
    fp_mul(&tmp1_fp,&tmp1_fp,&Pt2.z);
    fp_mul(&S1,&tmp1_fp,&Pt2.z);

    //S2
    fp_mul(&tmp1_fp,&Pt2.y,&Pt1.z);
    fp_mul(&tmp1_fp,&tmp1_fp,&Pt1.z);
    fp_mul(&S2,&tmp1_fp,&Pt1.z);

    //H
    fp_sub(&H,&U2,&U1);
    //r
    fp_sub(&r,&S2,&S1);

    //ANS->x
    fp_mul(&ANS->x,&r,&r);

    fp_mul(&tmp1_fp,&H,&H);
    fp_mul(&tmp1_fp,&tmp1_fp,&H);
    fp_sub(&ANS->x,&ANS->x,&tmp1_fp);

    fp_add(&tmp1_fp,&U1,&U1);
    fp_mul(&tmp1_fp,&tmp1_fp,&H);
    fp_mul(&tmp1_fp,&tmp1_fp,&H);
    fp_sub(&ANS->x,&ANS->x,&tmp1_fp);


    //ANS->y
    fp_mul(&tmp1_fp,&S1,&H);
    fp_mul(&tmp1_fp,&tmp1_fp,&H);
    fp_mul(&tmp1_fp,&tmp1_fp,&H);
    fp_set_neg(&ANS->y,&tmp1_fp);

    fp_mul(&tmp1_fp,&U1,&H);
    fp_mul(&tmp1_fp,&tmp1_fp,&H);
    fp_sub(&tmp1_fp,&tmp1_fp,&ANS->x);
    fp_mul(&tmp1_fp,&tmp1_fp,&r);
    fp_add(&ANS->y,&ANS->y,&tmp1_fp);
    
    //ANS->z
    fp_mul(&tmp1_fp,&Pt1.z,&Pt2.z);
    fp_mul(&ANS->z,&tmp1_fp,&H);
}
void efp_eca_mixture(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2){
    static efp_jacobian_t Pt1,Pt2;
    static fp_t Z1Z1,HH,I,J,V;
    static fp_t U1,U2,S1,S2,H,r;

    static fp_t bufL,tmpL1,tmpL2;
    static fp_t buf,tmp1,tmp2;
    fp_t out;
    if(P1->infinity==1){
        efp_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_jacobian_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    efp_jacobian_set(&Pt1,P1);
    efp_jacobian_set(&Pt2,P2);
    
    //Z1Z1
    fp_mul(&Z1Z1,&Pt1.z,&Pt1.z);
    
    //U2
    fp_mul(&U2,&Pt2.x,&Z1Z1);
    
    //S2
    fp_mul(&tmp1,&Z1Z1,&Pt1.z);
    fp_mul(&S2,&tmp1,&Pt2.y);
    
    //H
    fp_sub(&H,&U2,&Pt1.x);
    
    //HH
    fp_mul(&HH,&H,&H);
    
    //I
    fp_add(&I,&HH,&HH);
    fp_add(&I,&I,&I);
    
    //J
    fp_mul(&J,&HH,&H);
    
    //r
    fp_sub(&r,&S2,&Pt1.y);
    
    //V
    fp_mul(&V,&Pt1.x,&HH);
    
    //X3
    fp_mul(&tmp1,&r,&r);
    fp_add(&tmp2,&V,&V);
    fp_sub(&tmp1,&tmp1,&J);
    fp_sub(&ANS->x,&tmp1,&tmp2);
    
    //Y3
    fp_sub(&tmp1,&V,&ANS->x);
    fp_mul(&tmp1,&tmp1,&r);
    fp_mul(&tmp2,&Pt1.y,&J);
    fp_sub(&ANS->y,&tmp1,&tmp2);
    
    //Z3
    fp_mul(&ANS->z,&Pt1.z,&H);
}
void efp_eca_lazy(efp_t *ANS,efp_t *P1,efp_t *P2){
    static fp_t Fbuf,tmpF4,out;
    static mp_limb_t bufL[FPLIMB2],tmpL3[FPLIMB2],tmpL5[FPLIMB2],tmpL6[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB],tmp2[FPLIMB];
    static efp_t P1t,P2t;
    if(P1->infinity==1){
        efp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_lazy(ANS,P1);
            return;
        }
    }
    
    efp_set(&P1t,P1);
    efp_set(&P2t,P2);

    fp_sub_lazy(tmp1,FPLIMB,P2t.x.x0,FPLIMB,P1t.x.x0,FPLIMB);
    mpn_invert(tmp1,tmp1,prime);
    fp_sub_lazy(tmp2,FPLIMB,P2t.y.x0,FPLIMB,P1t.y.x0,FPLIMB);
    fp_mul_lazy(tmpL3,tmp1,tmp2);
    fp_mod(&tmpF4,tmpL3,FPLIMB2);
    fp_mul_lazy(tmpL5,tmpF4.x0,tmpF4.x0);


    fp_sub_lazy(tmpL6,FPLIMB2,tmpL5,FPLIMB2,P1t.x.x0,FPLIMB);
    fp_sub_lazy(ANS->x.x0,FPLIMB2,tmpL6,FPLIMB2,P2t.x.x0,FPLIMB);
    fp_mod(&ANS->x,ANS->x.x0,FPLIMB2);

    fp_sub_lazy(tmp1,FPLIMB,P1t.x.x0,FPLIMB,ANS->x.x0,FPLIMB);
    fp_mul_lazy(tmpL3,tmpF4.x0,tmp1);
    fp_sub_lazy(bufL,FPLIMB2,tmpL3,FPLIMB2,P1t.y.x0,FPLIMB);
    fp_mod(&ANS->y,bufL,FPLIMB2);

}
void efp_eca_projective_lazy(efp_projective_t *ANS,efp_projective_t *P1,efp_projective_t *P2){
    static efp_projective_t Pt1,Pt2;
    static mp_limb_t Y1Z2[FPLIMB],X1Z2[FPLIMB],Z1Z2[FPLIMB];
    static mp_limb_t u[FPLIMB],uu[FPLIMB];
    static mp_limb_t v[FPLIMB],vv[FPLIMB],vvv[FPLIMB],R[FPLIMB],A[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp[FPLIMB];
    
    fp_t out;
    if(P1->infinity==1){
        efp_projective_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_projective_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_projective_lazy(ANS,P1);
            return;
        }
    }
    mpn_zero(bufL,FPLIMB2);
    mpn_zero(tmpL,FPLIMB2);
    
    efp_projective_set(&Pt1,P1);
    efp_projective_set(&Pt2,P2);

    //Y1Z2
    fp_mul_lazy(bufL,Pt1.y.x0,Pt2.z.x0);
    mpn_mod(Y1Z2,bufL,FPLIMB2);

    //X1Z2
    fp_mul_lazy(bufL,Pt1.x.x0,Pt2.z.x0);
    mpn_mod(X1Z2,bufL,FPLIMB2);
    
    //Z1Z2
    fp_mul_lazy(bufL,Pt1.z.x0,Pt2.z.x0);
    mpn_mod(Z1Z2,bufL,FPLIMB2);
    
    //u
    fp_mul_lazy(tmpL,Pt1.z.x0,Pt2.y.x0);//Y2*Z1
    mpn_mod(tmp,tmpL,FPLIMB2);
    fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,Y1Z2,FPLIMB);
    mpn_mod(u,buf,FPLIMB);
    
    //uu
    fp_sqr_lazy(bufL,u);
    mpn_mod(uu,bufL,FPLIMB2);
    
    //v
    fp_mul_lazy(tmpL,Pt1.z.x0,Pt2.x.x0);//X2*Z1
    mpn_mod(tmp,tmpL,FPLIMB2);
    fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,X1Z2,FPLIMB);
    mpn_mod(v,buf,FPLIMB);
    
    //vv
    fp_sqr_lazy(bufL,v);
    mpn_mod(vv,bufL,FPLIMB2);
    
    //vv
    fp_mul_lazy(bufL,v,vv);
    mpn_mod(vvv,bufL,FPLIMB2);
    
    //R
    fp_mul_lazy(bufL,X1Z2,vv);
    mpn_mod(R,bufL,FPLIMB2);

    //A
    fp_mul_lazy(tmpL,Z1Z2,uu);
    mpn_mod(tmp,tmpL,FPLIMB2);
    fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,vvv,FPLIMB);
    fp_sub_lazy(tmp,FPLIMB,buf,FPLIMB,R,FPLIMB);
    fp_sub_lazy(buf,FPLIMB,tmp,FPLIMB,R,FPLIMB);
    mpn_mod(A,buf,FPLIMB);
    //ANS->x
    fp_mul_lazy(bufL,A,v);
    fp_mod(&ANS->x,bufL,FPLIMB2);

    //ANS->y
    fp_sub_lazy(tmp,FPLIMB,R,FPLIMB,A,FPLIMB);
    fp_mul_lazy(tmpL,tmp,u);
    fp_mul_lazy(bufL,vvv,Y1Z2);

    fp_sub_lazy_mod(&ANS->y,tmpL,bufL);
    
    //ANS->z
    fp_mul_lazy(bufL,Z1Z2,vvv);
    fp_mod(&ANS->z,bufL,FPLIMB2);
}
void efp_eca_jacobian_lazy(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2){
    static efp_jacobian_t Pt1,Pt2;
    static mp_limb_t U1[FPLIMB],U2[FPLIMB],S1[FPLIMB],S2[FPLIMB],H[FPLIMB],r[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2],tmpL2[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB];
    static mp_limb_t tmpZ1[FPLIMB],tmpZ2[FPLIMB],tmpH2[FPLIMB],tmpH3[FPLIMB],tmpU1H2[FPLIMB];
    fp_t out;
    if(P1->infinity==1){
        efp_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_jacobian_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    efp_jacobian_set(&Pt1,P1);
    efp_jacobian_set(&Pt2,P2);

    //U1
    fp_sqr_lazy(bufL,Pt2.z.x0);
    //fp_mul_lazy(bufL,Pt2.z.x0,Pt2.z.x0);
    mpn_mod(tmpZ2,bufL,FPLIMB2);
    fp_mul_lazy(tmpL1,tmpZ2,Pt1.x.x0);
    mpn_mod(U1,tmpL1,FPLIMB2);

    //U2
    fp_sqr_lazy(bufL,Pt1.z.x0);
    //fp_mul_lazy(bufL,Pt1.z.x0,Pt1.z.x0);
    mpn_mod(tmpZ1,bufL,FPLIMB2);
    fp_mul_lazy(tmpL1,tmpZ1,Pt2.x.x0);
    mpn_mod(U2,tmpL1,FPLIMB2);

    //S1
    fp_mul_lazy(bufL,tmpZ2,Pt2.z.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    fp_mul_lazy(tmpL1,tmp1,Pt1.y.x0);
    mpn_mod(S1,tmpL1,FPLIMB2);
    //gmp_printf("S1=%Nu\n",S1,FPLIMB);

    //S2
    fp_mul_lazy(bufL,tmpZ1,Pt1.z.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    fp_mul_lazy(tmpL1,tmp1,Pt2.y.x0);
    mpn_mod(S2,tmpL1,FPLIMB2);
    //gmp_printf("S2=%Nu\n",S2,FPLIMB);
    
    //H
    fp_sub_lazy(buf,FPLIMB,U2,FPLIMB,U1,FPLIMB);
    mpn_mod(H,buf,FPLIMB);

    //r
    fp_sub_lazy(buf,FPLIMB,S2,FPLIMB,S1,FPLIMB);
    mpn_mod(r,buf,FPLIMB);
    //gmp_printf("r=%Nu\n",r,FPLIMB);

    //ANS->x
    fp_sqr_lazy(tmpL2,r);
    //fp_mul_lazy(tmpL2,r,r);

    fp_sqr_lazy(tmpL1,H);
    //fp_mul_lazy(tmpL1,H,H);
    mpn_mod(tmpH2,tmpL1,FPLIMB2);
    fp_mul_lazy(tmpL1,tmpH2,H);
    mpn_mod(tmpH3,tmpL1,FPLIMB2);
    fp_sub_lazy(bufL,FPLIMB2,tmpL2,FPLIMB2,tmpL1,FPLIMB2);

    fp_mul_lazy(tmpL1,tmpH2,U1);
    mpn_mod(tmpU1H2,tmpL1,FPLIMB2);
    fp_add_lazy(tmpL1,FPLIMB2,tmpL1,FPLIMB2,tmpL1,FPLIMB2);
    fp_sub_lazy_mod(&ANS->x,bufL,tmpL1);

    //ANS->y
    fp_sub_lazy(tmp1,FPLIMB,tmpU1H2,FPLIMB,ANS->x.x0,FPLIMB);
    fp_mul_lazy(bufL,tmp1,r);

    fp_mul_lazy(tmpL1,tmpH3,S1);
    fp_sub_lazy_mod(&ANS->y,bufL,tmpL1);
    
    //ANS->z
    fp_mul_lazy(tmpL1,Pt1.z.x0,Pt2.z.x0);
    mpn_mod(tmp1,tmpL1,FPLIMB2);
    fp_mul_lazy(bufL,tmp1,H);
    fp_mod(&ANS->z,bufL,FPLIMB2);
}
void efp_eca_jacobian_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2){
    static efp_jacobian_t Pt1,Pt2;
    static fp_t U1,U2,S1,S2,H,r;

    static fp_t buf,tmp1,tmp2;
    static fp_t tmpZ1,tmpZ2,tmpH2,tmpH3,tmpU1H2;
    
    if(P1->infinity==1){
        efp_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_jacobian_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    efp_jacobian_set(&Pt1,P1);
    efp_jacobian_set(&Pt2,P2);

    //U1
    fp_mulmod_montgomery(&tmpZ2,&Pt2.z,&Pt2.z);
    fp_mulmod_montgomery(&U1,&tmpZ2,&Pt1.x);
    //fp_printf("U1=",&U1);printf("\n");

    //U2
    fp_mulmod_montgomery(&tmpZ1,&Pt1.z,&Pt1.z);
    fp_mulmod_montgomery(&U2,&tmpZ1,&Pt2.x);
    //fp_printf("U2=",&U2);printf("\n");
    
    //S1
    fp_mulmod_montgomery(&tmp1,&tmpZ2,&Pt2.z);
    fp_mulmod_montgomery(&S1,&tmp1,&Pt1.y);
    //fp_printf("S1=",&S1);printf("\n");

    //S2
    fp_mulmod_montgomery(&tmp1,&tmpZ1,&Pt1.z);
    fp_mulmod_montgomery(&S2,&tmp1,&Pt2.y);
    //fp_printf("S2=",&S2);printf("\n");
    
    //H
    //fp_printf("U1=",&U1);printf("\n");
    fp_sub(&H,&U2,&U1);
    //fp_printf("H=",&H);printf("\n");

    //r
    fp_sub(&r,&S2,&S1);
    //fp_printf("r=",&r);printf("\n");

    //ANS->x
    fp_mulmod_montgomery(&tmp1,&r,&r);
    
    fp_mulmod_montgomery(&tmpH2,&H,&H);
    fp_mulmod_montgomery(&tmpH3,&tmpH2,&H);
    fp_sub(&tmp2,&tmp1,&tmpH3);

    fp_mulmod_montgomery(&tmpU1H2,&tmpH2,&U1);
    fp_add(&tmp1,&tmpU1H2,&tmpU1H2);
    fp_sub(&ANS->x,&tmp2,&tmp1);

    //ANS->y
    fp_sub_lazy(tmp1.x0,FPLIMB,tmpU1H2.x0,FPLIMB,ANS->x.x0,FPLIMB);
    fp_mulmod_montgomery(&tmp1,&tmp1,&r);

    fp_mulmod_montgomery(&tmp2,&tmpH3,&S1);
    fp_sub(&ANS->y,&tmp1,&tmp2);
    
    //ANS->z
    fp_mulmod_montgomery(&tmp1,&Pt1.z,&Pt2.z);
    fp_mulmod_montgomery(&ANS->z,&tmp1,&H);
    //getchar();
}
void efp_eca_mixture_lazy(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2){
    static efp_jacobian_t Pt1,Pt2;
    static mp_limb_t Z1Z1[FPLIMB],HH[FPLIMB],I[FPLIMB],J[FPLIMB],V[FPLIMB];
    static mp_limb_t U1[FPLIMB],U2[FPLIMB],S1[FPLIMB],S2[FPLIMB],H[FPLIMB],r[FPLIMB];

    static mp_limb_t bufL[FPLIMB2],tmpL1[FPLIMB2],tmpL2[FPLIMB2];
    static mp_limb_t buf[FPLIMB],tmp1[FPLIMB],tmp2[FPLIMB];
    fp_t out;
    if(P1->infinity==1){
        efp_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_jacobian_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    efp_jacobian_set(&Pt1,P1);
    efp_jacobian_set(&Pt2,P2);
    /*
    fp_printf("\nX1=",&Pt1.x);printf("\n");
    fp_printf("Y1=",&Pt1.y);printf("\n");
    fp_printf("Z1=",&Pt1.z);printf("\n");
    fp_printf("X2=",&Pt2.x);printf("\n");
    fp_printf("Y2=",&Pt2.y);printf("\n");
    fp_printf("Z2=",&Pt2.z);printf("\n");
*/
/*
    printf("mixture lazy\n");
    gmp_printf("Pt1.x=%Nu\n",Pt1.x.x0,FPLIMB);
    gmp_printf("Pt1.y=%Nu\n",Pt1.y.x0,FPLIMB);
    gmp_printf("Pt1.z=%Nu\n",Pt1.z.x0,FPLIMB);
    gmp_printf("Pt2.x=%Nu\n",Pt2.x.x0,FPLIMB);
    gmp_printf("Pt2.y=%Nu\n",Pt2.y.x0,FPLIMB);
    gmp_printf("Pt2.z=%Nu\n",Pt2.z.x0,FPLIMB);
    */
    //Z1Z1
    fp_mul_lazy(bufL,Pt1.z.x0,Pt1.z.x0);
    mpn_mod(Z1Z1,bufL,FPLIMB2);
    //gmp_printf("Z1Z1=%Nu\n",Z1Z1,FPLIMB);
    
    //U2
    fp_mul_lazy(bufL,Pt2.x.x0,Z1Z1);
    mpn_mod(U2,bufL,FPLIMB2);
    //gmp_printf("U2=%Nu\n",U2,FPLIMB);
    
    //S2
    fp_mul_lazy(bufL,Z1Z1,Pt1.z.x0);
    mpn_mod(tmp1,bufL,FPLIMB2);
    fp_mul_lazy(tmpL1,tmp1,Pt2.y.x0);
    mpn_mod(S2,tmpL1,FPLIMB2);
    //gmp_printf("S2=%Nu\n",S2,FPLIMB);
    
    //H
    fp_sub_lazy(buf,FPLIMB,U2,FPLIMB,Pt1.x.x0,FPLIMB);
    mpn_mod(H,buf,FPLIMB);
    //gmp_printf("H=%Nu\n",H,FPLIMB);
    
    //HH
    fp_mul_lazy(tmpL1,H,H);
    mpn_mod(HH,tmpL1,FPLIMB2);
    //gmp_printf("HH=%Nu\n",HH,FPLIMB);
    
    //I
    fp_add_lazy(I,FPLIMB,HH,FPLIMB,HH,FPLIMB);
    fp_add_lazy(I,FPLIMB,I,FPLIMB,I,FPLIMB);
    mpn_mod(I,I,FPLIMB);
    //gmp_printf("I=%Nu\n",I,FPLIMB);
    
    //J
    //fp_mul_lazy(bufL,I,H);
    fp_mul_lazy(bufL,HH,H);
    mpn_mod(J,bufL,FPLIMB2);
    //gmp_printf("J=%Nu\n",J,FPLIMB);
    
    //r
    //gmp_printf("Y1=%Nu\n",Pt1.y.x0,FPLIMB);
    fp_sub_lazy(buf,FPLIMB,S2,FPLIMB,Pt1.y.x0,FPLIMB);
    //fp_add_lazy(buf,FPLIMB,buf,FPLIMB,buf,FPLIMB);
    mpn_mod(r,buf,FPLIMB);
    //gmp_printf("r=%Nu\n",r,FPLIMB);
    
    //V
    //fp_mul_lazy(bufL,Pt1.x.x0,I);
    fp_mul_lazy(bufL,Pt1.x.x0,HH);
    mpn_mod(V,bufL,FPLIMB2);
    //gmp_printf("V=%Nu\n",V,FPLIMB);
    
    //X3
    fp_mul_lazy(tmpL1,r,r);
    fp_add_lazy(tmp1,FPLIMB,V,FPLIMB,V,FPLIMB);
    fp_sub_lazy(bufL,FPLIMB2,tmpL1,FPLIMB2,J,FPLIMB);
    fp_sub_lazy(bufL,FPLIMB2,bufL,FPLIMB2,tmp1,FPLIMB);
    fp_mod(&ANS->x,bufL,FPLIMB2);
    
    //Y3
    fp_sub_lazy(tmp1,FPLIMB,V,FPLIMB,ANS->x.x0,FPLIMB);
    fp_mul_lazy(tmpL1,tmp1,r);
    fp_mul_lazy(tmpL2,Pt1.y.x0,J);
    //fp_add_lazy(tmpL2,FPLIMB2,tmpL2,FPLIMB2,tmpL2,FPLIMB2);
    fp_sub_lazy_mod(&ANS->y,tmpL1,tmpL2);
    
    //Z3
    fp_mul_lazy(bufL,Pt1.z.x0,H);
    fp_mod(&ANS->z,bufL,FPLIMB2);
    //getchar();

}

void efp_eca_mixture_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P1,efp_jacobian_t *P2){
    static efp_jacobian_t Pt1,Pt2;
    static fp_t Z1Z1,HH,I,J,V;
    static fp_t U1,U2,S1,S2,H,r;
    static fp_t buf,tmp1,tmp2;
    
    if(P1->infinity==1){
        efp_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp_jacobian_set(ANS,P1);
        return;
    }else if(fp_cmp(&P1->x,&P2->x)==0){
        if(fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp_ecd_jacobian_lazy(ANS,P1);
            return;
        }
    }
    
    efp_jacobian_set(&Pt1,P1);
    efp_jacobian_set(&Pt2,P2);
    
    
    
    //fp_printf("Pt1.x=",&Pt1.x);printf("\n");
    //fp_printf("Pt1.y=",&Pt1.y);printf("\n");
    //fp_printf("Pt1.z=",&Pt1.z);printf("\n");
    //fp_printf("Pt2.x=",&Pt2.x);printf("\n");
    //fp_printf("Pt2.y=",&Pt2.y);printf("\n");
    //fp_printf("Pt2.z=",&Pt2.z);printf("\n");
    
    //fp_printf_montgomery("Pt1.x=",&Pt1.x);printf("\n");
    //fp_printf_montgomery("Pt1.y=",&Pt1.y);printf("\n");
    //fp_printf_montgomery("Pt1.z=",&Pt1.z);printf("\n");
    //fp_printf_montgomery("Pt2.x=",&Pt2.x);printf("\n");
    //fp_printf_montgomery("Pt2.y=",&Pt2.y);printf("\n");
    //fp_printf_montgomery("Pt2.z=",&Pt2.z);printf("\n");
    
    
    //Z1Z1
    fp_mulmod_montgomery(&Z1Z1,&Pt1.z,&Pt1.z);
    
    //U2
    fp_mulmod_montgomery(&U2,&Pt2.x,&Z1Z1);
    
    //S2
    fp_mulmod_montgomery(&tmp1,&Z1Z1,&Pt1.z);
    fp_mulmod_montgomery(&S2,&tmp1,&Pt2.y);
    
    //H
    fp_sub(&H,&U2,&Pt1.x);
    
    //HH
    fp_mulmod_montgomery(&HH,&H,&H);
    
    //I
    fp_add(&I,&HH,&HH);
    fp_add(&I,&I,&I);
    
    //J
    fp_mulmod_montgomery(&J,&HH,&H);
    
    //r
    fp_sub(&r,&S2,&Pt1.y);
    
    //V
    fp_mulmod_montgomery(&V,&Pt1.x,&HH);
    
    //X3
    fp_mulmod_montgomery(&tmp1,&r,&r);
    fp_add(&tmp2,&V,&V);
    fp_sub(&buf,&tmp1,&J);
    fp_sub(&ANS->x,&buf,&tmp2);
    
    //Y3
    fp_sub_lazy(tmp1.x0,FPLIMB,V.x0,FPLIMB,ANS->x.x0,FPLIMB);
    fp_mulmod_montgomery(&tmp2,&tmp1,&r);
    fp_mulmod_montgomery(&tmp1,&Pt1.y,&J);
    fp_sub(&ANS->y,&tmp2,&tmp1);
    
    
    //ANS->z
    fp_mulmod_montgomery(&ANS->z,&Pt1.z,&H);

}


void efp_scm(efp_t *ANS,efp_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp_set(ANS,P);
        return;
    }
    
    efp_t Tmp_P,Next_P;
    efp_init(&Tmp_P);
    efp_set(&Tmp_P,P);
    efp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    efp_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp_ecd(&Next_P,&Next_P);
        if(binary[i]=='1'){
            efp_eca(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    efp_set(ANS,&Next_P);

}
void efp_scm_jacobian(efp_t *ANS,efp_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp_set(ANS,P);
        return;
    }
    
    efp_jacobian_t Tmp_P,Next_P;
    efp_affine_to_jacobian(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    efp_jacobian_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp_ecd_jacobian(&Next_P,&Next_P);
        if(binary[i]=='1'){
	    efp_eca_jacobian(&Next_P,&Next_P,&Tmp_P);
        }
    }
    efp_jacobian_to_affine(ANS,&Next_P);
}
void efp_scm_lazy(efp_t *ANS,efp_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp_set(ANS,P);
        return;
    }
    
    efp_t Tmp_P,Next_P;
    efp_init(&Tmp_P);
    efp_set(&Tmp_P,P);
    efp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    efp_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp_ecd_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
            efp_eca_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    efp_set(ANS,&Next_P);

}
void efp_scm_jacobian_lazy(efp_t *ANS,efp_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp_set(ANS,P);
        return;
    }
    
    efp_jacobian_t Tmp_P,Next_P;
    efp_affine_to_jacobian(&Tmp_P,P);

    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    efp_jacobian_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp_ecd_jacobian_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
	    efp_eca_jacobian_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    efp_jacobian_to_affine(ANS,&Next_P);
}
//skew frobenius map
void efp_skew_frobenius_map_p2(efp_t *ANS,efp_t *A){
    fp_mul_mpn(&ANS->x,&A->x,epsilon1);
    fp_set_neg(&ANS->y,&A->y);
}
void efp_jacobian_skew_frobenius_map_p2(efp_jacobian_t *ANS,efp_jacobian_t *A){
    fp_mul_mpn(&ANS->x,&A->x,epsilon1);
    fp_set_neg(&ANS->y,&A->y);
    fp_set(&ANS->z,&A->z);
}
void efp_jacobian_skew_frobenius_map_p2_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *A){
    static fp_t buf;
    //TODO:global
    mpn_to_montgomery(buf.x0,epsilon1);
    fp_mulmod_montgomery(&ANS->x,&A->x,&buf);
    fp_set_neg(&ANS->y,&A->y);
    fp_set(&ANS->z,&A->z);
}
