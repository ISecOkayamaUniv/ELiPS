#include <ELiPS/efp2.h>

void efp2_init(efp2_t *P){
    fp2_init(&P->x);
    fp2_init(&P->y);
    P->infinity=0;
}

void efp2_projective_init(efp2_projective_t *P){
    fp2_init(&P->x);
    fp2_init(&P->y);
    fp2_init(&P->z);
    P->infinity=0;
}

void efp2_jacobian_init(efp2_jacobian_t *P){
    fp2_init(&P->x);
    fp2_init(&P->y);
    fp2_init(&P->z);
    P->infinity=0;
}

void efp2_printf(char *str,efp2_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf("",&P->x);
        printf(",");
        fp2_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_printf_montgomery(char *str,efp2_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf_montgomery("",&P->x);
        printf(",");
        fp2_printf_montgomery("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_println(char *str,efp2_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf("",&P->x);
        printf(",");
        fp2_printf("",&P->y);
        printf(")\n");
    }else{
        printf("0\n");
    }
}

void efp2_jacobian_printf(char *str,efp2_jacobian_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf("",&P->x);
        printf(",");
        fp2_printf("",&P->y);
        printf(",");
        fp2_printf("",&P->z);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_projective_printf(char *str,efp2_projective_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf("",&P->x);
        printf(",");
        fp2_printf("",&P->y);
        printf(",");
        fp2_printf("",&P->z);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_jacobian_printf_montgomery(char *str,efp2_jacobian_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf_montgomery("",&P->x);
        printf(",");
        fp2_printf_montgomery("",&P->y);
        printf(",");
        fp2_printf_montgomery("",&P->z);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_projective_printf_montgomery(char *str,efp2_projective_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf_montgomery("",&P->x);
        printf(",");
        fp2_printf_montgomery("",&P->y);
        printf(",");
        fp2_printf_montgomery("",&P->z);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_projective_printf_affine(char *str,efp2_projective_t *P){
    static efp2_t out;
    efp2_projective_to_affine(&out,P);
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf("",&out.x);
        printf(",");
        fp2_printf("",&out.y);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_projective_printf_affine_montgomery(char *str,efp2_projective_t *P){
    static efp2_t out;
    efp2_projective_to_affine_montgomery(&out,P);
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp2_printf_montgomery("",&out.x);
        printf(",");
        fp2_printf_montgomery("",&out.y);
        printf(")");
    }else{
        printf("0");
    }
}

void efp2_set(efp2_t *ANS,efp2_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void efp2_jacobian_set(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set(&ANS->y,&A->y);
    fp2_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}

void efp2_projective_set(efp2_projective_t *ANS,efp2_projective_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set(&ANS->y,&A->y);
    fp2_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}

void efp2_affine_to_jacobian(efp2_jacobian_t *ANS,efp2_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set(&ANS->y,&A->y);
    fp2_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}

void efp2_affine_to_projective(efp2_projective_t *ANS,efp2_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set(&ANS->y,&A->y);
    fp2_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}

void efp2_affine_to_projective_montgomery(efp2_projective_t *ANS,efp2_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set(&ANS->y,&A->y);
    fp2_set_mpn(&ANS->z,RmodP);
    ANS->infinity=A->infinity;
}

void efp2_affine_to_jacobian_montgomery(efp2_jacobian_t *ANS,efp2_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set(&ANS->y,&A->y);
    fp2_set_mpn(&ANS->z,RmodP);
    ANS->infinity=A->infinity;
}

void efp2_jacobian_to_affine(efp2_t *ANS,efp2_jacobian_t *A){
    static fp2_t Zi,Zt;
    //TODO:mul->mul_lazy
    fp2_inv(&Zi,&A->z);
    fp2_mul(&Zt,&Zi,&Zi);
    fp2_mul(&ANS->x,&A->x,&Zt);
    fp2_mul(&Zt,&Zt,&Zi);
    fp2_mul(&ANS->y,&A->y,&Zt);
    ANS->infinity=A->infinity;
}

void efp2_projective_to_affine(efp2_t *ANS,efp2_projective_t *A){
    static fp2_t Zi;
    //TODO:mul->mul_lazy
    fp2_inv_lazy(&Zi,&A->z);
    fp2_mul_lazy(&ANS->x,&A->x,&Zi);
    fp2_mul_lazy(&ANS->y,&A->y,&Zi);
    ANS->infinity=A->infinity;
}

void efp2_jacobian_to_affine_montgomery(efp2_t *ANS,efp2_jacobian_t *A){
    static fp2_t Zi,Zt;
    fp2_inv_lazy_montgomery(&Zi,&A->z);
    fp2_mul_lazy_montgomery(&Zt,&Zi,&Zi);
    fp2_mul_lazy_montgomery(&ANS->x,&A->x,&Zt);
    fp2_mul_lazy_montgomery(&Zt,&Zt,&Zi);
    fp2_mul_lazy_montgomery(&ANS->y,&A->y,&Zt);
    ANS->infinity=A->infinity;
}

void efp2_projective_to_affine_montgomery(efp2_t *ANS,efp2_projective_t *A){
    static fp2_t Zi;
    //TODO:mul->mul_lazy
    fp2_inv_lazy_montgomery(&Zi,&A->z);
    fp2_mul_lazy_montgomery(&ANS->x,&A->x,&Zi);
    fp2_mul_lazy_montgomery(&ANS->y,&A->y,&Zi);
    ANS->infinity=A->infinity;
}

void efp2_mix(efp2_jacobian_t *ANS,efp2_jacobian_t *A,fp2_t *Zi){
    static fp2_t Zt;
    //TODO:mul->mul_lazy
    fp2_mul(&Zt,Zi,Zi);
    fp2_mul(&ANS->x,&A->x,&Zt);
    fp2_mul(&Zt,&Zt,Zi);
    fp2_mul(&ANS->y,&A->y,&Zt);
    fp2_set_ui(&ANS->z,1);
    ANS->infinity=A->infinity;
}

void efp2_mix_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A,fp2_t *Zi){
    static fp2_t Zt;
    fp2_mul_lazy_montgomery(&Zt,Zi,Zi);
    fp2_mul_lazy_montgomery(&ANS->x,&A->x,&Zt);
    fp2_mul_lazy_montgomery(&Zt,&Zt,Zi);
    fp2_mul_lazy_montgomery(&ANS->y,&A->y,&Zt);
    fp2_set_mpn(&ANS->z,RmodP);
    ANS->infinity=A->infinity;
}

void efp2_set_ui(efp2_t *ANS,unsigned long int UI){
    fp2_set_ui(&ANS->x,UI);
    fp2_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}
void efp2_to_montgomery(efp2_t *ANS,efp2_t *A){
    fp2_to_montgomery(&ANS->x,&A->x);
    fp2_to_montgomery(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void efp2_projective_to_montgomery(efp2_projective_t *ANS,efp2_projective_t *A){
    fp2_to_montgomery(&ANS->x,&A->x);
    fp2_to_montgomery(&ANS->y,&A->y);
    fp2_to_montgomery(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void efp2_mod_montgomery(efp2_t *ANS,efp2_t *A){
    fp2_mod_montgomery(&ANS->x,&A->x);
    fp2_mod_montgomery(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void efp2_projective_mod_montgomery(efp2_projective_t *ANS,efp2_projective_t *A){
    fp2_mod_montgomery(&ANS->x,&A->x);
    fp2_mod_montgomery(&ANS->y,&A->y);
    fp2_mod_montgomery(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}
void efp2_set_mpn(efp2_t *ANS,mp_limb_t *A){
    fp2_set_mpn(&ANS->x,A);
    fp2_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}

void efp2_set_neg(efp2_t *ANS,efp2_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}
void efp2_jacobian_set_neg(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
    fp2_set(&ANS->x,&A->x);
    fp2_set_neg(&ANS->y,&A->y);
    fp2_set(&ANS->z,&A->z);
    ANS->infinity=A->infinity;
}

int  efp2_cmp(efp2_t *A,efp2_t *B){
    if(fp2_cmp(&A->x,&B->x)==0 && fp2_cmp(&A->y,&B->y)==0){
        return 0;
    }else if(A->infinity==1&&B->infinity==1){
	    return 0;
    }else{
        return 1;
    }
}

void efp2_rational_point(efp2_t *P){
    fp2_t tmp1,tmp2;
    fp2_init(&tmp1);
    fp2_init(&tmp2);

    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));

    while(1){
        fp2_set_random(&P->x,state);
        fp2_sqr(&tmp1,&P->x);
        fp2_mul(&tmp2,&tmp1,&P->x);
        fp_add_mpn(&tmp2.x0,&tmp2.x0,curve_b);
        if(fp2_legendre(&tmp2)==1){
            fp2_sqrt(&P->y,&tmp2);
            break;
        }
    }
    P->infinity=0;
}

void efp2_ecd(efp2_t *ANS,efp2_t *P){
    static efp2_t tmp1_efp2;
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    if(fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }

    efp2_set(&tmp1_efp2,P);

    fp2_add(&tmp1_fp2,&tmp1_efp2.y,&tmp1_efp2.y);

    fp2_inv(&tmp1_fp2,&tmp1_fp2);
    fp2_sqr(&tmp2_fp2,&tmp1_efp2.x);
    fp2_add(&tmp3_fp2,&tmp2_fp2,&tmp2_fp2);
    fp2_add(&tmp2_fp2,&tmp2_fp2,&tmp3_fp2);
    fp2_mul(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);

    fp2_sqr(&tmp1_fp2,&tmp3_fp2);
    fp2_add(&tmp2_fp2,&tmp1_efp2.x,&tmp1_efp2.x);
    fp2_sub(&ANS->x,&tmp1_fp2,&tmp2_fp2);

    fp2_sub(&tmp1_fp2,&tmp1_efp2.x,&ANS->x);
    fp2_mul(&tmp2_fp2,&tmp3_fp2,&tmp1_fp2);
    fp2_sub(&ANS->y,&tmp2_fp2,&tmp1_efp2.y);
}

void efp2_ecd_lazy_montgomery(efp2_t *ANS,efp2_t *P){
    static efp2_t tmp1_efp2;
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    if(fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }

    efp2_set(&tmp1_efp2,P);

    fp2_add(&tmp1_fp2,&tmp1_efp2.y,&tmp1_efp2.y);

    fp2_inv_lazy_montgomery(&tmp1_fp2,&tmp1_fp2);
    fp2_sqr_lazy_montgomery(&tmp2_fp2,&tmp1_efp2.x);
    fp2_add(&tmp3_fp2,&tmp2_fp2,&tmp2_fp2);
    fp2_add(&tmp2_fp2,&tmp2_fp2,&tmp3_fp2);
    fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);

    fp2_sqr_lazy_montgomery(&tmp1_fp2,&tmp3_fp2);
    fp2_add(&tmp2_fp2,&tmp1_efp2.x,&tmp1_efp2.x);
    fp2_sub(&ANS->x,&tmp1_fp2,&tmp2_fp2);

    fp2_sub(&tmp1_fp2,&tmp1_efp2.x,&ANS->x);
    fp2_mul_lazy_montgomery(&tmp2_fp2,&tmp3_fp2,&tmp1_fp2);
    fp2_sub(&ANS->y,&tmp2_fp2,&tmp1_efp2.y);
}
/*
void efp2_ecd_projective_lazy(efp2_projective_t *ANS,efp2_projective_t *P){
    static efp2_projective_t Pt;
    static fp2_t XX,ZZ;
    static fp2_t w,s,ss,sss;
    static fp2_t R,RR,B,h;

    static fp2_t bufL,tmpL;
    static fp2_t buf,tmp;

    if(fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }

    efp2_projective_set(&Pt,P);


    //XX
    fp2_sqr_lazy(&XX,&Pt.x);

    //ZZ
    fp2_sqr_lazy(&ZZ,&Pt.z);

    //w
    fp2_add(&tmp,&XX,&XX);
    fp2_add(&w,&tmp,&XX);

    //s
    fp2_mul_lazy(&tmp,&Pt.y,&Pt.z);//Y2*Z1
    fp2_add(&s,&tmp,&tmp);

    //ss
    fp2_sqr_lazy(&ss,&s);

    //sss
    fp2_mul_lazy(&sss,&s,&ss);

    //R
    fp2_mul_lazy(&R,&Pt.y,&s);//Y2*Z1

    //RR
    fp2_sqr_lazy(&RR,&R);

    //B
    fp2_add_nonmod_single(&tmp,&Pt.x,&R);
    fp2_sqr_lazy(&tmp,&tmp);
    fp2_sub(&tmp,&tmp,&XX);
    fp2_sub(&B,&tmp,&RR);

    //h
    fp2_sqr_lazy(&tmp,&w);
    fp2_sub(&tmp,&tmp,&B);
    fp2_sub(&h,&tmp,&B);

    //ANS->x
    fp2_mul_lazy(&ANS->x,&h,&s);

    //ANS->y
    fp2_sub_lazy(&tmp,&B,&h);
    fp2_mul_lazy(&tmp,&tmp,&w);
    fp2_sub(&tmp,&tmp,&RR);
    fp2_sub(&ANS->y,&tmp,&RR);

    //ANS->z
    fp2_set(&ANS->z,&sss);
}
*/
void efp2_ecd_jacobian_lazy_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *P){
    static fp2_t s,m,T;

    static fp2_t buf,tmp1;
    static fp2_t tmpY2;
    static efp2_jacobian_t Pt;
    if(fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }

    efp2_jacobian_set(&Pt,P);

    //s
    fp2_mul_lazy_montgomery(&tmpY2,&Pt.y,&Pt.y);
    fp2_mul_lazy_montgomery(&tmp1,&tmpY2,&Pt.x);
    fp2_add(&tmp1,&tmp1,&tmp1);
    fp2_add(&s,&tmp1,&tmp1);

    //m
    fp2_add_nonmod_single(&tmp1,&Pt.x,&Pt.x);
    fp2_add_nonmod_single(&tmp1,&tmp1,&Pt.x);
    fp2_mul_lazy_montgomery(&m,&tmp1,&Pt.x);

    //T
    fp2_mul_lazy_montgomery(&T,&m,&m);
    fp2_add(&tmp1,&s,&s);
    fp2_sub(&T,&T,&tmp1);

    //ANS->x
    fp2_set(&ANS->x,&T);

    //ANS->y
    fp2_sub_nonmod_single(&tmp1,&s,&T);
    fp2_mul_lazy_montgomery(&buf,&tmp1,&m);

    fp2_mul_lazy_montgomery(&tmp1,&tmpY2,&tmpY2);
    fp2_add(&tmp1,&tmp1,&tmp1);
    fp2_add(&tmp1,&tmp1,&tmp1);
    fp2_add(&tmp1,&tmp1,&tmp1);
    fp2_sub(&ANS->y,&buf,&tmp1);

    //ANS->z
    fp2_add_nonmod_single(&tmp1,&Pt.y,&Pt.y);
    fp2_mul_lazy_montgomery(&ANS->z,&tmp1,&Pt.z);
}
void efp2_eca(efp2_t *ANS,efp2_t *P1,efp2_t *P2){
    static efp2_t tmp1_efp2,tmp2_efp2;
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    if(P1->infinity==1){
        efp2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp2_set(ANS,P1);
        return;
    }else if(fp2_cmp(&P1->x,&P2->x)==0){
        if(fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp2_ecd(ANS,P1);
            return;
        }
    }

    efp2_set(&tmp1_efp2,P1);
    efp2_set(&tmp2_efp2,P2);

    fp2_sub(&tmp1_fp2,&tmp2_efp2.x,&tmp1_efp2.x);
    fp2_inv(&tmp1_fp2,&tmp1_fp2);
    fp2_sub(&tmp2_fp2,&tmp2_efp2.y,&tmp1_efp2.y);
    fp2_mul(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);
    fp2_sqr(&tmp1_fp2,&tmp3_fp2);
    fp2_sub(&tmp2_fp2,&tmp1_fp2,&tmp1_efp2.x);
    fp2_sub(&ANS->x,&tmp2_fp2,&tmp2_efp2.x);
    fp2_sub(&tmp1_fp2,&tmp1_efp2.x,&ANS->x);
    fp2_mul(&tmp2_fp2,&tmp3_fp2,&tmp1_fp2);
    fp2_sub(&ANS->y,&tmp2_fp2,&tmp1_efp2.y);
}

void efp2_eca_lazy_montgomery(efp2_t *ANS,efp2_t *P1,efp2_t *P2){
    static efp2_t tmp1_efp2,tmp2_efp2;
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    if(P1->infinity==1){
        efp2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp2_set(ANS,P1);
        return;
    }else if(fp2_cmp(&P1->x,&P2->x)==0){
        if(fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp2_ecd_lazy_montgomery(ANS,P1);
            return;
        }
    }

    efp2_set(&tmp1_efp2,P1);
    efp2_set(&tmp2_efp2,P2);

    fp2_sub(&tmp1_fp2,&tmp2_efp2.x,&tmp1_efp2.x);
    fp2_inv_lazy_montgomery(&tmp1_fp2,&tmp1_fp2);
    fp2_sub(&tmp2_fp2,&tmp2_efp2.y,&tmp1_efp2.y);
    fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp1_fp2,&tmp2_fp2);
    fp2_sqr_lazy_montgomery(&tmp1_fp2,&tmp3_fp2);
    fp2_sub(&tmp2_fp2,&tmp1_fp2,&tmp1_efp2.x);
    fp2_sub(&ANS->x,&tmp2_fp2,&tmp2_efp2.x);
    fp2_sub(&tmp1_fp2,&tmp1_efp2.x,&ANS->x);
    fp2_mul_lazy_montgomery(&tmp2_fp2,&tmp3_fp2,&tmp1_fp2);
    fp2_sub(&ANS->y,&tmp2_fp2,&tmp1_efp2.y);
}
/*
void efp2_eca_projective_lazy(efp2_projective_t *ANS,efp2_projective_t *P1,efp2_projective_t *P2){
    static efp2_projective_t Pt1,Pt2;
    static fp2_t Y1Z2,X1Z2,Z1Z2;
    static fp2_t u,uu;
    static fp2_t v,vv,vvv,R,A;

    static fp2_t bufL,tmpL;
    static fp2_t buf,tmp;

    if(P1->infinity==1){
        efp2_projective_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp2_projective_set(ANS,P1);
        return;
    }else if(fp2_cmp(&P1->x,&P2->x)==0){
        if(fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp2_ecd_projective_lazy(ANS,P1);
            return;
        }
    }

    efp2_projective_set(&Pt1,P1);
    efp2_projective_set(&Pt2,P2);

    //Y1Z2
    fp2_mul_lazy(&Y1Z2,&Pt1.y,&Pt2.z);

    //X1Z2
    fp2_mul_lazy(&X1Z2,&Pt1.x,&Pt2.z);

    //Z1Z2
    fp2_mul_lazy(&Z1Z2,&Pt1.z,&Pt2.z);

    //u
    fp2_mul_lazy(&tmp,&Pt1.z,&Pt2.y);//Y2*Z1
    fp2_sub(&u,&tmp,&Y1Z2);

    //uu
    fp2_sqr_lazy(&uu,&u);

    //v
    fp2_mul_lazy(&tmp,&Pt1.z,&Pt2.x);//X2*Z1
    fp2_sub(&v,&tmp,&X1Z2);

    //vv
    fp2_sqr_lazy(&vv,&v);

    //vv
    fp2_mul_lazy(&vvv,&v,&vv);

    //R
    fp2_mul_lazy(&R,&X1Z2,&vv);

    //A
    fp2_mul_lazy(&tmp,&Z1Z2,&uu);
    fp2_sub_lazy(&buf,&tmp,&vvv);
    fp2_sub_lazy(&tmp,&buf,&R);
    fp2_sub_lazy(&A,&tmp,&R);
    //ANS->x
    fp2_mul_lazy(&ANS->x,&A,&v);

    //ANS->y
    fp2_sub_lazy(&tmp,&R,&A);
    fp2_mul_lazy(&tmpL,&tmp,&u);
    fp2_mul_lazy(&buf,&vvv,&Y1Z2);

    fp2_sub(&ANS->y,&tmpL,&buf);

    //ANS->z
    fp2_mul_lazy(&ANS->z,&Z1Z2,&vvv);
}
*/
void efp2_eca_jacobian_lazy_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *P1,efp2_jacobian_t *P2){
    static efp2_jacobian_t Pt1,Pt2;
    static fp2_t U1,U2,S1,S2,H,r;

    static fp2_t buf,tmp1,tmp2;
    static fp2_t tmpZ1,tmpZ2,tmpH2,tmpH3,tmpU1H2;

    if(P1->infinity==1){
        efp2_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp2_jacobian_set(ANS,P1);
        return;
    }else if(fp2_cmp(&P1->x,&P2->x)==0){
        if(fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp2_ecd_jacobian_lazy_montgomery(ANS,P1);
            return;
        }
    }

    efp2_jacobian_set(&Pt1,P1);
    efp2_jacobian_set(&Pt2,P2);

    //U1
    fp2_mul_lazy_montgomery(&tmpZ2,&Pt2.z,&Pt2.z);
    fp2_mul_lazy_montgomery(&U1,&tmpZ2,&Pt1.x);
    //fp2_printf("U1=",&U1);printf("\n");

    //U2
    fp2_mul_lazy_montgomery(&tmpZ1,&Pt1.z,&Pt1.z);
    fp2_mul_lazy_montgomery(&U2,&tmpZ1,&Pt2.x);
    //fp2_printf("U2=",&U2);printf("\n");

    //S1
    fp2_mul_lazy_montgomery(&tmp1,&tmpZ2,&Pt2.z);
    fp2_mul_lazy_montgomery(&S1,&tmp1,&Pt1.y);
    //fp2_printf("S1=",&S1);printf("\n");

    //S2
    fp2_mul_lazy_montgomery(&tmp1,&tmpZ1,&Pt1.z);
    fp2_mul_lazy_montgomery(&S2,&tmp1,&Pt2.y);
    //fp2_printf("S2=",&S2);printf("\n");

    //H
    //fp2_printf("U1=",&U1);printf("\n");
    fp2_sub(&H,&U2,&U1);
    //fp2_printf("H=",&H);printf("\n");

    //r
    fp2_sub(&r,&S2,&S1);
    //fp2_printf("r=",&r);printf("\n");

    //ANS->x
    fp2_mul_lazy_montgomery(&tmp1,&r,&r);

    fp2_mul_lazy_montgomery(&tmpH2,&H,&H);
    fp2_mul_lazy_montgomery(&tmpH3,&tmpH2,&H);
    fp2_sub(&tmp2,&tmp1,&tmpH3);

    fp2_mul_lazy_montgomery(&tmpU1H2,&tmpH2,&U1);
    fp2_add(&tmp1,&tmpU1H2,&tmpU1H2);
    fp2_sub(&ANS->x,&tmp2,&tmp1);

    //ANS->y
    fp2_sub_nonmod_single(&tmp1,&tmpU1H2,&ANS->x);
    fp2_mul_lazy_montgomery(&tmp1,&tmp1,&r);

    fp2_mul_lazy_montgomery(&tmp2,&tmpH3,&S1);
    fp2_sub(&ANS->y,&tmp1,&tmp2);

    //ANS->z
    fp2_mul_lazy_montgomery(&tmp1,&Pt1.z,&Pt2.z);
    fp2_mul_lazy_montgomery(&ANS->z,&tmp1,&H);
    //getchar();
}
void efp2_eca_mixture_lazy_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *P1,efp2_jacobian_t *P2){
    static efp2_jacobian_t Pt1,Pt2;
    static fp2_t Z1Z1,HH,I,J,V;
    static fp2_t U1,U2,S1,S2,H,r;
    static fp2_t buf,tmp1,tmp2;

    if(P1->infinity==1){
        efp2_jacobian_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp2_jacobian_set(ANS,P1);
        return;
    }else if(fp2_cmp(&P1->x,&P2->x)==0){
        if(fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp2_ecd_jacobian_lazy_montgomery(ANS,P1);
            return;
        }
    }

    efp2_jacobian_set(&Pt1,P1);
    efp2_jacobian_set(&Pt2,P2);

    //Z1Z1
    fp2_mul_lazy_montgomery(&Z1Z1,&Pt1.z,&Pt1.z);

    //U2
    fp2_mul_lazy_montgomery(&U2,&Pt2.x,&Z1Z1);

    //S2
    fp2_mul_lazy_montgomery(&tmp1,&Z1Z1,&Pt1.z);
    fp2_mul_lazy_montgomery(&S2,&tmp1,&Pt2.y);

    //H
    fp2_sub(&H,&U2,&Pt1.x);

    //HH
    fp2_mul_lazy_montgomery(&HH,&H,&H);

    //I
    fp2_add(&I,&HH,&HH);
    fp2_add(&I,&I,&I);

    //J
    fp2_mul_lazy_montgomery(&J,&HH,&H);

    //r
    fp2_sub(&r,&S2,&Pt1.y);

    //V
    fp2_mul_lazy_montgomery(&V,&Pt1.x,&HH);

    //X3
    fp2_mul_lazy_montgomery(&tmp1,&r,&r);
    fp2_add(&tmp2,&V,&V);
    fp2_sub(&buf,&tmp1,&J);
    fp2_sub(&ANS->x,&buf,&tmp2);

    //Y3
    fp2_sub_nonmod_single(&tmp1,&V,&ANS->x);
    fp2_mul_lazy_montgomery(&tmp2,&tmp1,&r);
    fp2_mul_lazy_montgomery(&tmp1,&Pt1.y,&J);
    fp2_sub(&ANS->y,&tmp2,&tmp1);


    //ANS->z
    fp2_mul_lazy_montgomery(&ANS->z,&Pt1.z,&H);

}

void efp2_scm(efp2_t *ANS,efp2_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp2_set(ANS,P);
        return;
    }

    efp2_t Tmp_P,Next_P;
    efp2_init(&Tmp_P);
    efp2_set(&Tmp_P,P);
    efp2_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);

    efp2_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp2_ecd(&Next_P,&Next_P);
        if(binary[i]=='1'){
            efp2_eca(&Next_P,&Next_P,&Tmp_P);
        }
    }
    efp2_set(ANS,&Next_P);
}

//skew_frobenius_map
void efp2_skew_frobenius_map_p1(efp2_t *ANS,efp2_t *A){
#ifdef TWIST_PHI_INV
    //x(w2w3)
    fp_set(&ANS->x.x0,&A->x.x0);
    fp_set_neg(&ANS->x.x1,&A->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p1][1]);
    //y(w8w9)
    fp_set(&ANS->y.x0,&A->y.x0);
    fp_set_neg(&ANS->y.x1,&A->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p1][4]);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp_set(&ANS->x.x0,&ANS->x.x0);
    fp_set_neg(&ANS->x.x1,&ANS->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p1][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp_set(&ANS->y.x0,&ANS->y.x0);
    fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p1][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);
#endif
}

void efp2_skew_frobenius_map_p2(efp2_t *ANS,efp2_t *A){
#ifdef TWIST_PHI_INV
    //x
    fp2_mul(&ANS->x,&A->x,&frobenius_constant[f_p2][1]);
    //y
    fp2_mul(&ANS->y,&A->y,&frobenius_constant[f_p2][4]);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p2][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p2][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);
#endif

}

void efp2_skew_frobenius_map_p3(efp2_t *ANS,efp2_t *A){
#ifdef TWIST_PHI_INV
    //x
    fp_set(&ANS->x.x0,&A->x.x0);
    fp_set_neg(&ANS->x.x1,&A->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p3][1]);
    //y
    fp_set(&ANS->y.x0,&A->y.x0);
    fp_set_neg(&ANS->y.x1,&A->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p3][4]);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp_set(&ANS->x.x0,&ANS->x.x0);
    fp_set_neg(&ANS->x.x1,&ANS->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p3][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp_set(&ANS->y.x0,&ANS->y.x0);
    fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p3][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);
#endif
}

void efp2_jacobian_skew_frobenius_map_p1(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
#ifdef TWIST_PHI_INV
    //x
    fp_set(&ANS->x.x0,&A->x.x0);
    fp_set_neg(&ANS->x.x1,&A->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p1][1]);
    //y
    fp_set(&ANS->y.x0,&A->y.x0);
    fp_set_neg(&ANS->y.x1,&A->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p1][4]);
    //z
    fp2_set(&ANS->z,&A->z);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp_set(&ANS->x.x0,&ANS->x.x0);
    fp_set_neg(&ANS->x.x1,&ANS->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p1][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp_set(&ANS->y.x0,&ANS->y.x0);
    fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p1][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);

    fp2_set(&ANS->z,&A->z);
#endif
}

void efp2_jacobian_skew_frobenius_map_p2(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
#ifdef TWIST_PHI_INV
    //x
    fp2_mul(&ANS->x,&A->x,&frobenius_constant[f_p2][1]);
    //y
    fp2_mul(&ANS->y,&A->y,&frobenius_constant[f_p2][4]);
    //z
    fp2_set(&ANS->z,&A->z);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp_set(&ANS->x.x0,&ANS->x.x0);
    fp_set_neg(&ANS->x.x1,&ANS->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p1][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp_set(&ANS->y.x0,&ANS->y.x0);
    fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p1][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);

    fp2_set(&ANS->z,&A->z);
#endif
}

void efp2_jacobian_skew_frobenius_map_p3(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
    //x
    fp_set(&ANS->x.x0,&A->x.x0);
    fp_set_neg(&ANS->x.x1,&A->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant[f_p3][1]);
    //y
    fp_set(&ANS->y.x0,&A->y.x0);
    fp_set_neg(&ANS->y.x1,&A->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant[f_p3][4]);
    //z
    fp2_set(&ANS->z,&A->z);
}
void efp2_skew_frobenius_map_p1_montgomery(efp2_t *ANS,efp2_t *A){
#ifdef TWIST_PHI_INV
    //x(w2w3)
    fp_set(&ANS->x.x0,&A->x.x0);
    fp_set_neg(&ANS->x.x1,&A->x.x1);
    fp2_mul_lazy_montgomery(&ANS->x,&ANS->x,&frobenius_constant_montgomery[f_p1][1]);
    //y(w8w9)
    fp_set(&ANS->y.x0,&A->y.x0);
    fp_set_neg(&ANS->y.x1,&A->y.x1);
    fp2_mul_lazy_montgomery(&ANS->y,&ANS->y,&frobenius_constant_montgomery[f_p1][4]);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp_set(&ANS->x.x0,&ANS->x.x0);
    fp_set_neg(&ANS->x.x1,&ANS->x.x1);
    fp2_mul(&ANS->x,&ANS->x,&frobenius_constant_montgomery[f_p1][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp_set(&ANS->y.x0,&ANS->y.x0);
    fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    fp2_mul(&ANS->y,&ANS->y,&frobenius_constant_montgomery[f_p1][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);
#endif
}
void efp2_jacobian_skew_frobenius_map_p1_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
#ifdef TWIST_PHI_INV
    //x
    fp_set(&ANS->x.x0,&A->x.x0);
    fp_set_neg(&ANS->x.x1,&A->x.x1);
    fp2_mul_lazy_montgomery(&ANS->x,&ANS->x,&frobenius_constant_montgomery[f_p1][1]);
    //y
    fp_set(&ANS->y.x0,&A->y.x0);
    fp_set_neg(&ANS->y.x1,&A->y.x1);
    fp2_mul_lazy_montgomery(&ANS->y,&ANS->y,&frobenius_constant_montgomery[f_p1][4]);
    //z
    fp2_set(&ANS->z,&A->z);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp_set(&ANS->x.x0,&ANS->x.x0);
    fp_set_neg(&ANS->x.x1,&ANS->x.x1);
    fp2_mul_lazy_montgomery(&ANS->x,&ANS->x,&frobenius_constant_montgomery[f_p1][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp_set(&ANS->y.x0,&ANS->y.x0);
    fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    fp2_mul_lazy_montgomery(&ANS->y,&ANS->y,&frobenius_constant_montgomery[f_p1][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);

    fp2_set(&ANS->z,&A->z);
#endif
}

void efp2_jacobian_skew_frobenius_map_p2_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
#ifdef TWIST_PHI_INV
    //x
    fp2_mul_lazy_montgomery(&ANS->x,&A->x,&frobenius_constant_montgomery[f_p2][1]);
    //y
    fp2_mul_lazy_montgomery(&ANS->y,&A->y,&frobenius_constant_montgomery[f_p2][4]);
    //z
    fp2_set(&ANS->z,&A->z);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp2_mul_lazy_montgomery(&ANS->x,&ANS->x,&frobenius_constant_montgomery[f_p2][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp2_mul_lazy_montgomery(&ANS->y,&ANS->y,&frobenius_constant_montgomery[f_p2][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);

    fp2_set(&ANS->z,&A->z);
#endif
}

void efp2_jacobian_skew_frobenius_map_p3_montgomery(efp2_jacobian_t *ANS,efp2_jacobian_t *A){
#ifdef TWIST_PHI_INV
    //x
    fp_set(&ANS->x.x0,&A->x.x0);
    fp_set_neg(&ANS->x.x1,&A->x.x1);
    fp2_mul_lazy_montgomery(&ANS->x,&ANS->x,&frobenius_constant_montgomery[f_p3][1]);
    //y
    fp_set(&ANS->y.x0,&A->y.x0);
    fp_set_neg(&ANS->y.x1,&A->y.x1);
    fp2_mul_lazy_montgomery(&ANS->y,&ANS->y,&frobenius_constant_montgomery[f_p3][4]);
    //z
    fp2_set(&ANS->z,&A->z);
#endif
#ifdef TWIST_PHI
    //x(w4w5)
    fp2_inv_basis(&ANS->x,&A->x);
    fp_set(&ANS->x.x0,&ANS->x.x0);
    fp_set_neg(&ANS->x.x1,&ANS->x.x1);
    fp2_mul_lazy_montgomery(&ANS->x,&ANS->x,&frobenius_constant_montgomery[f_p3][2]);
    fp2_mul_basis(&ANS->x,&ANS->x);

    //y
    fp2_inv_basis(&ANS->y,&A->y);
    fp_set(&ANS->y.x0,&ANS->y.x0);
    fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    fp2_mul_lazy_montgomery(&ANS->y,&ANS->y,&frobenius_constant_montgomery[f_p3][4]);
    fp2_mul_basis(&ANS->y,&ANS->y);

    fp2_set(&ANS->z,&A->z);
#endif
}
void efp2_skew_frobenius_map_p10(efp2_t *ANS,efp2_t *A){
    //x
    fp2_mul(&ANS->x,&A->x,&frobenius_constant[f_p10][1]);
    //y
    fp2_mul(&ANS->y,&A->y,&frobenius_constant[f_p10][4]);
}
