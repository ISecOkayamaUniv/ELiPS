#include <ELiPS/efp.h>

void efp_init(efp_t *P) {
  fp_init(&P->x);
  fp_init(&P->y);
  P->infinity = 1;
}

void efp_projective_init(efp_projective_t *P) {
  fp_init(&P->x);
  fp_init(&P->y);
  fp_init(&P->z);
  P->infinity = 1;
}

void efp_jacobian_init(efp_jacobian_t *P) {
  fp_init(&P->x);
  fp_init(&P->y);
  fp_init(&P->z);
  //P->infinity=0;
  P->infinity = 1;
}

void efp_printf(char *str, efp_t *P) {
  printf("%s", str);
  if (P->infinity == 0) {
    printf("(");
    fp_printf("", &P->x);
    printf(",");
    fp_printf("", &P->y);
    printf(")");
  } else {
    printf("0");
  }
}

void efp_println(char *str, efp_t *P) {
  printf("%s", str);
  if (P->infinity == 0) {
    printf("(");
    fp_printf("", &P->x);
    printf(",");
    fp_printf("", &P->y);
    printf(")\n");
  } else {
    printf("0\n");
  }
}

void efp_projective_printf(char *str, efp_projective_t *P) {
  printf("%s", str);
  if (P->infinity == 0) {
    printf("(");
    fp_printf("", &P->x);
    printf(",");
    fp_printf("", &P->y);
    printf(",");
    fp_printf("", &P->z);
    printf(")");
  } else {
    printf("Infinity");
  }
}

void efp_jacobian_printf(char *str, efp_jacobian_t *P) {
  printf("%s", str);
  if (P->infinity == 0) {
    printf("(");
    fp_printf("", &P->x);
    printf(",");
    fp_printf("", &P->y);
    printf(",");
    fp_printf("", &P->z);
    printf(")");
  } else {
    printf("Infinity");
  }
}

void efp_set(efp_t *ANS, efp_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}

void efp_projective_set(efp_projective_t *ANS, efp_projective_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set(&ANS->y, &A->y);
  fp_set(&ANS->z, &A->z);
  ANS->infinity = A->infinity;
}

void efp_jacobian_set(efp_jacobian_t *ANS, efp_jacobian_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set(&ANS->y, &A->y);
  fp_set(&ANS->z, &A->z);
  ANS->infinity = A->infinity;
}

void efp_affine_to_projective(efp_projective_t *ANS, efp_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set(&ANS->y, &A->y);
  fp_set_ui(&ANS->z, 1);
  ANS->infinity = A->infinity;
}

void efp_affine_to_jacobian(efp_jacobian_t *ANS, efp_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set(&ANS->y, &A->y);
  fp_set_ui(&ANS->z, 1);
  ANS->infinity = A->infinity;
}

void efp_affine_to_jacobian_montgomery(efp_jacobian_t *ANS, efp_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set(&ANS->y, &A->y);
  mpn_copyd(ANS->z.x0, RmodP, FPLIMB);
  ANS->infinity = A->infinity;
}

void efp_projective_to_affine(efp_t *ANS, efp_projective_t *A) {
  fp_t Zi;
  fp_inv(&Zi, &A->z);
  fp_mul(&ANS->x, &A->x, &Zi);
  fp_mul(&ANS->y, &A->y, &Zi);
  ANS->infinity = A->infinity;
}

void efp_jacobian_to_affine(efp_t *ANS, efp_jacobian_t *A) {
  static fp_t Zi, Zt;
  fp_inv(&Zi, &A->z);
  fp_sqr(&Zt, &Zi);
  fp_mul(&ANS->x, &A->x, &Zt);
  fp_mul(&Zt, &Zt, &Zi);
  fp_mul(&ANS->y, &A->y, &Zt);
  ANS->infinity = A->infinity;
}

void efp_jacobian_to_affine_montgomery(efp_t *ANS, efp_jacobian_t *A) {
  fp_t Zi, Zt;
  fp_inv_montgomery(&Zi, &A->z);
  fp_sqrmod_montgomery(&Zt, &Zi);
  fp_mulmod_montgomery(&ANS->x, &A->x, &Zt);
  fp_mulmod_montgomery(&Zt, &Zt, &Zi);
  fp_mulmod_montgomery(&ANS->y, &A->y, &Zt);
  ANS->infinity = A->infinity;
}

void efp_jacobian_to_mixture_noninv(efp_jacobian_t *ANS, efp_jacobian_t *A, fp_t *Zi) {
  static fp_t Zt;
  fp_sqr(&Zt, Zi);
  fp_mul(&ANS->x, &A->x, &Zt);
  fp_mul(&Zt, &Zt, Zi);
  fp_mul(&ANS->y, &A->y, &Zt);
  fp_set_ui(&ANS->z, 1);
  ANS->infinity = A->infinity;
}

void efp_jacobian_to_mixture_noninv_montgomery(efp_jacobian_t *ANS, efp_jacobian_t *A, fp_t *Zi) {
  static fp_t Zt;
  fp_sqrmod_montgomery(&Zt, Zi);
  fp_mulmod_montgomery(&ANS->x, &A->x, &Zt);
  fp_mulmod_montgomery(&Zt, &Zt, Zi);
  fp_mulmod_montgomery(&ANS->y, &A->y, &Zt);
  //fp_set_ui(&ANS->z,1);
  mpn_copyd(ANS->z.x0, RmodP, FPLIMB);
  ANS->infinity = A->infinity;
}

void efp_to_montgomery(efp_t *ANS, efp_t *A) {
  fp_to_montgomery(&ANS->x, &A->x);
  fp_to_montgomery(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}

void efp_mod_montgomery(efp_t *ANS, efp_t *A) {
  fp_mod_montgomery(&ANS->x, &A->x);
  fp_mod_montgomery(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}

void efp_set_ui(efp_t *ANS, unsigned long int UI) {
  fp_set_ui(&ANS->x, UI);
  fp_set_ui(&ANS->y, UI);
  ANS->infinity = 0;
}

void efp_set_mpn(efp_t *ANS, mp_limb_t *A) {
  fp_set_mpn(&ANS->x, A);
  fp_set_mpn(&ANS->y, A);
  ANS->infinity = 0;
}

void efp_set_neg(efp_t *ANS, efp_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set_neg(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}

void efp_projective_set_neg(efp_projective_t *ANS, efp_projective_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set_neg(&ANS->y, &A->y);
  fp_set(&ANS->z, &A->z);
  ANS->infinity = A->infinity;
}

void efp_jacobian_set_neg(efp_jacobian_t *ANS, efp_jacobian_t *A) {
  fp_set(&ANS->x, &A->x);
  fp_set_neg(&ANS->y, &A->y);
  fp_set(&ANS->z, &A->z);
  ANS->infinity = A->infinity;
}

int efp_cmp(efp_t *A, efp_t *B) {
  if (!fp_cmp(&A->x, &B->x) && !fp_cmp(&A->y, &B->y))
    return 0;
  else if (A->infinity && B->infinity)
    return 0;
  else
    return 1;
}

void efp_set_random(efp_t *P) {
  fp_t tmp1, tmp2, tmp_x;
  fp_init(&tmp1);
  fp_init(&tmp2);
  fp_init(&tmp_x);

  while (1) {
    fp_set_random(&P->x, state);
    fp_mul(&tmp1, &P->x, &P->x);
    fp_mul(&tmp2, &tmp1, &P->x);
    fp_add_mpn(&tmp_x, &tmp2, curve_b);
    if (fp_legendre(&tmp_x) == 1) {
      fp_sqrt(&P->y, &tmp_x);
      break;
    }
  }
  P->infinity = 0;
}

void efp_ecd(efp_t *ANS, efp_t *P) {
  static efp_t tmp1_efp;
  static fp_t tmp1_fp, tmp2_fp, tmp3_fp;
  if (fp_cmp_zero(&P->y) == 0 || P->infinity == 1) {
    ANS->infinity = 1;
  } else {
    efp_set(&tmp1_efp, P);

    fp_add(&tmp1_fp, &tmp1_efp.y, &tmp1_efp.y);
    fp_inv(&tmp1_fp, &tmp1_fp);

    fp_sqr(&tmp2_fp, &tmp1_efp.x);
    fp_add(&tmp3_fp, &tmp2_fp, &tmp2_fp);
    fp_add(&tmp2_fp, &tmp2_fp, &tmp3_fp);

    fp_mul(&tmp3_fp, &tmp1_fp, &tmp2_fp);
    fp_sqr(&tmp1_fp, &tmp3_fp);

    fp_add(&tmp2_fp, &tmp1_efp.x, &tmp1_efp.x);
    fp_sub(&ANS->x, &tmp1_fp, &tmp2_fp);

    fp_sub(&tmp1_fp, &tmp1_efp.x, &ANS->x);
    fp_mul(&tmp2_fp, &tmp3_fp, &tmp1_fp);
    fp_sub(&ANS->y, &tmp2_fp, &tmp1_efp.y);
    ANS->infinity = 0;
  }
}

void efp_ecd_lazy_montgomery(efp_t *ANS, efp_t *P) {
  static efp_t tmp1_efp;
  static fp_t tmp1_fp, tmp2_fp, tmp3_fp;
  if (fp_cmp_zero(&P->y) == 0 || P->infinity == 1) {
    ANS->infinity = 1;
  } else {
    efp_set(&tmp1_efp, P);

    fp_add(&tmp1_fp, &tmp1_efp.y, &tmp1_efp.y);
    fp_inv_montgomery(&tmp1_fp, &tmp1_fp);

    fp_sqrmod_montgomery(&tmp2_fp, &tmp1_efp.x);
    fp_add(&tmp3_fp, &tmp2_fp, &tmp2_fp);
    fp_add(&tmp2_fp, &tmp2_fp, &tmp3_fp);

    fp_mulmod_montgomery(&tmp3_fp, &tmp1_fp, &tmp2_fp);
    fp_sqrmod_montgomery(&tmp1_fp, &tmp3_fp);

    fp_add(&tmp2_fp, &tmp1_efp.x, &tmp1_efp.x);
    fp_sub(&ANS->x, &tmp1_fp, &tmp2_fp);

    fp_sub(&tmp1_fp, &tmp1_efp.x, &ANS->x);
    fp_mulmod_montgomery(&tmp2_fp, &tmp3_fp, &tmp1_fp);
    fp_sub(&ANS->y, &tmp2_fp, &tmp1_efp.y);
    ANS->infinity = 0;
  }
}
void efp_ecd_jacobian_lazy_montgomery(efp_jacobian_t *ANS, efp_jacobian_t *P) {
  static fp_t s, m, T;

  static fp_t buf, tmp1;
  static fp_t tmpY2;
  static efp_jacobian_t Pt;
  if (fp_cmp_zero(&P->y) == 0 || P->infinity == 1) {
    ANS->infinity = 1;
  } else {
    efp_jacobian_set(&Pt, P);

    //s
    fp_sqrmod_montgomery(&tmpY2, &Pt.y);

    fp_mulmod_montgomery(&tmp1, &tmpY2, &Pt.x);
    fp_add(&tmp1, &tmp1, &tmp1);
    fp_add(&s, &tmp1, &tmp1);
    //m
    fp_add_nonmod_single(&tmp1, &Pt.x, &Pt.x);
    fp_add_nonmod_single(&tmp1, &tmp1, &Pt.x);
    fp_mulmod_montgomery(&m, &tmp1, &Pt.x);

    //T
    fp_sqrmod_montgomery(&T, &m);
    fp_add(&tmp1, &s, &s);
    fp_sub(&T, &T, &tmp1);

    //ANS->x
    fp_set(&ANS->x, &T);

    //ANS->y
    fp_sub_nonmod_single(&tmp1, &s, &T);
    fp_mulmod_montgomery(&buf, &tmp1, &m);

    fp_sqrmod_montgomery(&tmp1, &tmpY2);
    fp_add(&tmp1, &tmp1, &tmp1);
    fp_add(&tmp1, &tmp1, &tmp1);
    fp_add(&tmp1, &tmp1, &tmp1);
    fp_sub(&ANS->y, &buf, &tmp1);

    //ANS->z
    fp_add_nonmod_single(&tmp1, &Pt.y, &Pt.y);
    fp_mulmod_montgomery(&ANS->z, &tmp1, &Pt.z);
    ANS->infinity = 0;
  }
}

// void efp_ecd_jacobian_lazy_montgomery(efp_jacobian_t *ANS,efp_jacobian_t *P){
//     static fp_t s,m,T;

//     static fp_t buf,tmp1;
//     static fp_t tmpY2;
//     static efp_jacobian_t Pt;

//     efp_jacobian_set(&Pt,P);

//     //s
//     fp_sqrmod_montgomery(&tmpY2,&Pt.y);

//     fp_mulmod_montgomery(&tmp1,&tmpY2,&Pt.x);
//     fp_add(&tmp1,&tmp1,&tmp1);
//     fp_add(&s,&tmp1,&tmp1);
//     //m
//     fp_add_nonmod_single(&tmp1,&Pt.x,&Pt.x);
//     fp_add_nonmod_single(&tmp1,&tmp1,&Pt.x);
//     fp_mulmod_montgomery(&m,&tmp1,&Pt.x);

//     //T
//     fp_sqrmod_montgomery(&T,&m);
//     fp_add(&tmp1,&s,&s);
//     fp_sub(&T,&T,&tmp1);

//     //ANS->x
//     fp_set(&ANS->x,&T);

//     //ANS->y
//     fp_sub_nonmod_single(&tmp1,&s,&T);
//     fp_mulmod_montgomery(&buf,&tmp1,&m);

//     fp_sqrmod_montgomery(&tmp1,&tmpY2);
//     fp_add(&tmp1,&tmp1,&tmp1);
//     fp_add(&tmp1,&tmp1,&tmp1);
//     fp_add(&tmp1,&tmp1,&tmp1);
//     fp_sub(&ANS->y,&buf,&tmp1);

//     //ANS->z
//     fp_add_nonmod_single(&tmp1,&Pt.y,&Pt.y);
//     fp_mulmod_montgomery(&ANS->z,&tmp1,&Pt.z);
//     ANS->infinity=0;
// }
void efp_eca(efp_t *ANS, efp_t *P1, efp_t *P2) {
  static efp_t tmp1_efp, tmp2_efp;
  static fp_t tmp1_fp, tmp2_fp, tmp3_fp;
  if (P1->infinity == 1) {
    efp_set(ANS, P2);
    return;
  } else if (P2->infinity == 1) {
    efp_set(ANS, P1);
    return;
  } else if (fp_cmp(&P1->x, &P2->x) == 0) {
    if (fp_cmp(&P1->y, &P2->y) != 0) {
      ANS->infinity = 1;
      return;
    } else {
      efp_ecd(ANS, P1);
      return;
    }
  } else {
    efp_set(&tmp1_efp, P1);
    efp_set(&tmp2_efp, P2);

    fp_sub(&tmp1_fp, &tmp2_efp.x, &tmp1_efp.x);
    fp_inv(&tmp1_fp, &tmp1_fp);
    fp_sub(&tmp2_fp, &tmp2_efp.y, &tmp1_efp.y);
    fp_mul(&tmp3_fp, &tmp1_fp, &tmp2_fp);
    fp_mul(&tmp1_fp, &tmp3_fp, &tmp3_fp);

    fp_sub(&tmp2_fp, &tmp1_fp, &tmp1_efp.x);
    fp_sub(&ANS->x, &tmp2_fp, &tmp2_efp.x);

    fp_sub(&tmp1_fp, &tmp1_efp.x, &ANS->x);
    fp_mul(&tmp2_fp, &tmp3_fp, &tmp1_fp);
    fp_sub(&ANS->y, &tmp2_fp, &tmp1_efp.y);
    ANS->infinity = 0;
  }
}
void efp_eca_lazy_montgomery(efp_t *ANS, efp_t *P1, efp_t *P2) {
  static efp_t tmp1_efp, tmp2_efp;
  static fp_t tmp1_fp, tmp2_fp, tmp3_fp;
  if (P1->infinity == 1) {
    efp_set(ANS, P2);
    return;
  } else if (P2->infinity == 1) {
    efp_set(ANS, P1);
    return;
  } else if (fp_cmp(&P1->x, &P2->x) == 0) {
    if (fp_cmp(&P1->y, &P2->y) != 0) {
      ANS->infinity = 1;
      return;
    } else {
      efp_ecd_lazy_montgomery(ANS, P1);
      return;
    }
  } else {
    efp_set(&tmp1_efp, P1);
    efp_set(&tmp2_efp, P2);

    fp_sub(&tmp1_fp, &tmp2_efp.x, &tmp1_efp.x);
    fp_inv_montgomery(&tmp1_fp, &tmp1_fp);
    fp_sub(&tmp2_fp, &tmp2_efp.y, &tmp1_efp.y);
    fp_mulmod_montgomery(&tmp3_fp, &tmp1_fp, &tmp2_fp);
    fp_mulmod_montgomery(&tmp1_fp, &tmp3_fp, &tmp3_fp);

    fp_sub(&tmp2_fp, &tmp1_fp, &tmp1_efp.x);
    fp_sub(&ANS->x, &tmp2_fp, &tmp2_efp.x);

    fp_sub(&tmp1_fp, &tmp1_efp.x, &ANS->x);
    fp_mulmod_montgomery(&tmp2_fp, &tmp3_fp, &tmp1_fp);
    fp_sub(&ANS->y, &tmp2_fp, &tmp1_efp.y);
    ANS->infinity = 0;
  }
}
/*
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
*/
void efp_eca_jacobian_lazy_montgomery(efp_jacobian_t *ANS, efp_jacobian_t *P1, efp_jacobian_t *P2) {
  static efp_jacobian_t Pt1, Pt2;
  static fp_t U1, U2, S1, S2, H, r;

  static fp_t buf, tmp1, tmp2;
  static fp_t tmpZ1, tmpZ2, tmpH2, tmpH3, tmpU1H2;

  if (P1->infinity == 1) {
    efp_jacobian_set(ANS, P2);
  } else if (P2->infinity == 1) {
    efp_jacobian_set(ANS, P1);
  } else if (fp_cmp(&P1->x, &P2->x) == 0) {
    if (fp_cmp(&P1->y, &P2->y) != 0) {
      ANS->infinity = 1;
    } else {
      efp_ecd_jacobian_lazy_montgomery(ANS, P1);
    }
  } else {
    efp_jacobian_set(&Pt1, P1);
    efp_jacobian_set(&Pt2, P2);

    //U1
    fp_sqrmod_montgomery(&tmpZ2, &Pt2.z);
    fp_mulmod_montgomery(&U1, &tmpZ2, &Pt1.x);
    //fp_printf("U1=",&U1);printf("\n");

    //U2
    fp_sqrmod_montgomery(&tmpZ1, &Pt1.z);
    fp_mulmod_montgomery(&U2, &tmpZ1, &Pt2.x);
    //fp_printf("U2=",&U2);printf("\n");

    //S1
    fp_mulmod_montgomery(&tmp1, &tmpZ2, &Pt2.z);
    fp_mulmod_montgomery(&S1, &tmp1, &Pt1.y);
    //fp_printf("S1=",&S1);printf("\n");

    //S2
    fp_mulmod_montgomery(&tmp1, &tmpZ1, &Pt1.z);
    fp_mulmod_montgomery(&S2, &tmp1, &Pt2.y);
    //fp_printf("S2=",&S2);printf("\n");

    //H
    //fp_printf("U1=",&U1);printf("\n");
    fp_sub(&H, &U2, &U1);
    //fp_printf("H=",&H);printf("\n");

    //r
    fp_sub(&r, &S2, &S1);
    //fp_printf("r=",&r);printf("\n");

    //ANS->x
    fp_sqrmod_montgomery(&tmp1, &r);

    fp_sqrmod_montgomery(&tmpH2, &H);
    fp_mulmod_montgomery(&tmpH3, &tmpH2, &H);
    fp_sub(&tmp2, &tmp1, &tmpH3);

    fp_mulmod_montgomery(&tmpU1H2, &tmpH2, &U1);
    fp_add(&tmp1, &tmpU1H2, &tmpU1H2);
    fp_sub(&ANS->x, &tmp2, &tmp1);

    //ANS->y
    //fp_sub_lazy(tmp1.x0,FPLIMB,tmpU1H2.x0,FPLIMB,ANS->x.x0,FPLIMB);
    fp_sub_nonmod_single(&tmp1, &tmpU1H2, &ANS->x);
    fp_mulmod_montgomery(&tmp1, &tmp1, &r);

    fp_mulmod_montgomery(&tmp2, &tmpH3, &S1);
    fp_sub(&ANS->y, &tmp1, &tmp2);

    //ANS->z
    fp_mulmod_montgomery(&tmp1, &Pt1.z, &Pt2.z);
    fp_mulmod_montgomery(&ANS->z, &tmp1, &H);
    //getchar();

    ANS->infinity = 0;
  }
}

void efp_eca_mixture_lazy_montgomery(efp_jacobian_t *ANS, efp_jacobian_t *P1, efp_jacobian_t *P2) {
  static efp_jacobian_t Pt1, Pt2;
  static fp_t Z1Z1, HH, I, J, V;
  static fp_t U1, U2, S1, S2, H, r;
  static fp_t buf, tmp1, tmp2;

  if (P1->infinity == 1) {
    efp_jacobian_set(ANS, P2);
  } else if (P2->infinity == 1) {
    efp_jacobian_set(ANS, P1);
  } else if (fp_cmp(&P1->x, &P2->x) == 0) {
    if (fp_cmp(&P1->y, &P2->y) != 0) {
      ANS->infinity = 1;
    } else {
      efp_ecd_jacobian_lazy_montgomery(ANS, P1);
    }
  } else {
    efp_jacobian_set(&Pt1, P1);
    efp_jacobian_set(&Pt2, P2);

    //Z1Z1
    fp_sqrmod_montgomery(&Z1Z1, &Pt1.z);

    //U2
    fp_mulmod_montgomery(&U2, &Pt2.x, &Z1Z1);

    //S2
    fp_mulmod_montgomery(&tmp1, &Z1Z1, &Pt1.z);
    fp_mulmod_montgomery(&S2, &tmp1, &Pt2.y);

    //H
    fp_sub(&H, &U2, &Pt1.x);

    //HH
    fp_sqrmod_montgomery(&HH, &H);

    //I
    fp_add(&I, &HH, &HH);
    fp_add(&I, &I, &I);

    //J
    fp_mulmod_montgomery(&J, &HH, &H);

    //r
    fp_sub(&r, &S2, &Pt1.y);

    //V
    fp_mulmod_montgomery(&V, &Pt1.x, &HH);

    //X3
    fp_sqrmod_montgomery(&tmp1, &r);
    fp_add(&tmp2, &V, &V);
    fp_sub(&buf, &tmp1, &J);
    fp_sub(&ANS->x, &buf, &tmp2);

    //Y3
    //fp_sub_lazy(tmp1.x0,FPLIMB,V.x0,FPLIMB,ANS->x.x0,FPLIMB);
    fp_sub_nonmod_single(&tmp1, &V, &ANS->x);
    fp_mulmod_montgomery(&tmp2, &tmp1, &r);
    fp_mulmod_montgomery(&tmp1, &Pt1.y, &J);
    fp_sub(&ANS->y, &tmp2, &tmp1);

    //ANS->z
    fp_mulmod_montgomery(&ANS->z, &Pt1.z, &H);

    ANS->infinity = 0;
  }
}

void efp_eca_mixture_lazy_montgomery_ignore_inf(efp_jacobian_t *ANS, efp_jacobian_t *P1, efp_jacobian_t *P2) {
  static efp_jacobian_t Pt1, Pt2;
  static fp_t Z1Z1, HH, I, J, V;
  static fp_t U1, U2, S1, S2, H, r;
  static fp_t buf, tmp1, tmp2;

  // if(P1->infinity==1){
  //     efp_jacobian_set(ANS,P2);
  //     return;
  // }else if(P2->infinity==1){
  //     efp_jacobian_set(ANS,P1);
  //     return;
  // }else if(fp_cmp(&P1->x,&P2->x)==0){
  //     if(fp_cmp(&P1->y,&P2->y)!=0){
  //         ANS->infinity=1;
  //         return;
  //     }else{
  //         efp_ecd_jacobian_lazy_montgomery(ANS,P1);
  //         return;
  //     }
  // }

  efp_jacobian_set(&Pt1, P1);
  efp_jacobian_set(&Pt2, P2);

  //Z1Z1
  fp_sqrmod_montgomery(&Z1Z1, &Pt1.z);

  //U2
  fp_mulmod_montgomery(&U2, &Pt2.x, &Z1Z1);

  //S2
  fp_mulmod_montgomery(&tmp1, &Z1Z1, &Pt1.z);
  fp_mulmod_montgomery(&S2, &tmp1, &Pt2.y);

  //H
  fp_sub(&H, &U2, &Pt1.x);

  //HH
  fp_sqrmod_montgomery(&HH, &H);

  //I
  fp_add(&I, &HH, &HH);
  fp_add(&I, &I, &I);

  //J
  fp_mulmod_montgomery(&J, &HH, &H);

  //r
  fp_sub(&r, &S2, &Pt1.y);

  //V
  fp_mulmod_montgomery(&V, &Pt1.x, &HH);

  //X3
  fp_sqrmod_montgomery(&tmp1, &r);
  fp_add(&tmp2, &V, &V);
  fp_sub(&buf, &tmp1, &J);
  fp_sub(&ANS->x, &buf, &tmp2);

  //Y3
  //fp_sub_lazy(tmp1.x0,FPLIMB,V.x0,FPLIMB,ANS->x.x0,FPLIMB);
  fp_sub_nonmod_single(&tmp1, &V, &ANS->x);
  fp_mulmod_montgomery(&tmp2, &tmp1, &r);
  fp_mulmod_montgomery(&tmp1, &Pt1.y, &J);
  fp_sub(&ANS->y, &tmp2, &tmp1);

  //ANS->z
  fp_mulmod_montgomery(&ANS->z, &Pt1.z, &H);
  ANS->infinity = 0;
}

void efp_scm(efp_t *ANS, efp_t *P, mpz_t scalar) {
  if (mpz_cmp_ui(scalar, 0) == 0) {
    ANS->infinity = 1;
    return;
  } else if (mpz_cmp_ui(scalar, 1) == 0) {
    efp_set(ANS, P);
    return;
  }

  efp_t Tmp_P, Next_P;
  efp_init(&Tmp_P);
  efp_set(&Tmp_P, P);
  efp_init(&Next_P);
  int i, length;
  length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);

  efp_set(&Next_P, &Tmp_P);
  for (i = 1; i < length; i++) {
    efp_ecd(&Next_P, &Next_P);
    if (binary[i] == '1') {
      efp_eca(&Next_P, &Next_P, &Tmp_P);
    }
  }

  efp_set(ANS, &Next_P);
}

void efp_scm_jacobian_lazy_montgomery(efp_jacobian_t *ANS, efp_jacobian_t *P, mpz_t scalar) {
  if (mpz_cmp_ui(scalar, 0) == 0) {
    ANS->infinity = 1;
    return;
  } else if (mpz_cmp_ui(scalar, 1) == 0) {
    efp_jacobian_set(ANS, P);
    return;
  }

  efp_jacobian_t tmp_ANS, tmp_P;
  efp_jacobian_init(&tmp_ANS);
  efp_jacobian_init(&tmp_P);
  int length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);

  efp_jacobian_set(&tmp_P, P);
  efp_jacobian_set(&tmp_ANS, &tmp_P);
  for (int i = 1; i < length; i++) {
    efp_ecd_jacobian_lazy_montgomery(&tmp_ANS, &tmp_ANS);
    if (binary[i] == '1') {
      efp_eca_jacobian_lazy_montgomery(&tmp_ANS, &tmp_ANS, &tmp_P);
    }
  }

  efp_jacobian_set(ANS, &tmp_ANS);
}

//skew frobenius map
void efp_skew_frobenius_map_p2(efp_t *ANS, efp_t *A) {
  fp_mul_mpn(&ANS->x, &A->x, epsilon1);
  fp_set_neg(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}
void efp_skew_frobenius_map_p2_montgomery(efp_t *ANS, efp_t *A) {
  static fp_t buf;
  fp_mulmod_montgomery(&ANS->x, &A->x, &epsilon1_montgomery);
  fp_set_neg(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}
void efp_jacobian_skew_frobenius_map_p2(efp_jacobian_t *ANS, efp_jacobian_t *A) {
  fp_mul_mpn(&ANS->x, &A->x, epsilon1);
  fp_set_neg(&ANS->y, &A->y);
  fp_set(&ANS->z, &A->z);
  ANS->infinity = A->infinity;
}
void efp_jacobian_skew_frobenius_map_p2_montgomery(efp_jacobian_t *ANS, efp_jacobian_t *A) {
  static fp_t buf;
  fp_mulmod_montgomery(&ANS->x, &A->x, &epsilon1_montgomery);
  fp_set_neg(&ANS->y, &A->y);
  fp_set(&ANS->z, &A->z);
  ANS->infinity = A->infinity;
}

void efp_set_random_fast(efp_t *P, gmp_randstate_t state) {
  fp_t tmp1, tmp2, tmp_x;
  fp_init(&tmp1);
  fp_init(&tmp2);
  fp_init(&tmp_x);

  while (1) {
    fp_set_random(&P->x, state);
    fp_mul(&tmp1, &P->x, &P->x);
    fp_mul(&tmp2, &tmp1, &P->x);
    fp_add_mpn(&tmp_x, &tmp2, curve_b);
    // if(fp_legendre_montgomery(&tmp_x)==1){
    //     fp_sqrt_montgomery(&P->y,&tmp_x);
    //     break;
    // }
    if (fp_legendre_sqrt(&P->y, &tmp_x) == 1) break;
  }
  P->infinity = 0;
}
