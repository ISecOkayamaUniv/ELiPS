#include <ELiPS/efp12.h>
//efp12
void efp12_init(efp12_t *P) {
  fp12_init(&P->x);
  fp12_init(&P->y);
  P->infinity = 0;
}
void sym_init(sym_t *A) {
  efp12_init(&A->p);
  efp12_init(&A->q);
}
void efp12_printf(char *str, efp12_t *P) {
  printf("%s", str);
  if (P->infinity == 0) {
    printf("(");
    fp12_printf("", &P->x);
    printf(",");
    fp12_printf("", &P->y);
    printf(")");
  } else {
    printf("0");
  }
}
void efp12_println(char *str, efp12_t *P) {
  printf("%s", str);
  if (P->infinity == 0) {
    printf("(");
    fp12_printf("", &P->x);
    printf(",");
    fp12_printf("", &P->y);
    printf(")\n");
  } else {
    printf("0\n");
  }
}
void efp12_set(efp12_t *ANS, efp12_t *A) {
  fp12_set(&ANS->x, &A->x);
  fp12_set(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}

void efp12_set_ui(efp12_t *ANS, unsigned long int UI1, unsigned long int UI2) {
  fp12_set_ui(&ANS->x, UI1);
  fp12_set_ui(&ANS->y, UI2);
  ANS->infinity = 0;
}

void efp12_set_mpn(efp12_t *ANS, mp_limb_t *A) {
  fp12_set_mpn(&ANS->x, A);
  fp12_set_mpn(&ANS->y, A);
  ANS->infinity = 0;
}

void efp12_set_neg(efp12_t *ANS, efp12_t *A) {
  fp12_set(&ANS->x, &A->x);
  fp12_set_neg(&ANS->y, &A->y);
  ANS->infinity = A->infinity;
}

int efp12_cmp(efp12_t *A, efp12_t *B) {
  if (fp12_cmp(&A->x, &B->x) == 0 && fp12_cmp(&A->y, &B->y) == 0) {
    return 0;
  } else if (A->infinity == 1 && B->infinity == 1) {
    return 0;
  } else {
    return 1;
  }
}

void efp12_rational_point(efp12_t *P) {
  fp12_t tmp1, tmp2;
  fp12_init(&tmp1);
  fp12_init(&tmp2);

  //gmp_randinit_default (state);
  //gmp_randseed_ui(state,(unsigned long)time(NULL));
  //gmp_randseed_ui(state,1);

  while (1) {
    fp12_set_random(&P->x, state);
    fp12_sqr(&tmp1, &P->x);
    fp12_mul(&tmp2, &tmp1, &P->x);
    fp_add_mpn(&tmp2.x0.x0.x0, &tmp2.x0.x0.x0, curve_b);
    if (fp12_legendre(&tmp2) == 1) {
      fp12_sqrt(&P->y, &tmp2);
      break;
    }
  }
  P->infinity = 0;
}

void bls12_generate_g1(efp12_t *P) {
  efp_t tmp_P;
  efp_init(&tmp_P);
  mpz_t exp;
  mpz_init(exp);

  efp_set_random(&tmp_P);
  efp12_set_ui(P, 0, 0);
  mpz_tdiv_q(exp, efp_total, order_z);
  efp_scm(&tmp_P, &tmp_P, exp);
  fp_set(&P->x.x0.x0.x0, &tmp_P.x);
  fp_set(&P->y.x0.x0.x0, &tmp_P.y);
  P->infinity = tmp_P.infinity;

  mpz_clear(exp);
}

void bls12_generate_g2(efp12_t *Q) {
  efp12_t random_P, P, frobenius_P;
  efp12_init(&random_P);
  efp12_init(&P);
  efp12_init(&frobenius_P);
  mpz_t exp;
  mpz_init(exp);

  efp12_rational_point(&random_P);
  mpz_pow_ui(exp, order_z, 2);
  mpz_tdiv_q(exp, efp12_total, exp);
  efp12_scm(&P, &random_P, exp);
  fp12_frobenius_map_p1(&frobenius_P.x, &P.x);
  fp12_frobenius_map_p1(&frobenius_P.y, &P.y);
  efp12_set_neg(&P, &P);
  efp12_eca(Q, &P, &frobenius_P);

  mpz_clear(exp);
}

void bls12_generate_symmetric_point(sym_t *A, mpz_t a) {
  efp12_t tmp_p;
  efp12_t tmp_q;
  mpz_t s;

  mpz_init(s);

  mpz_urandomm(a, state, order_z);
  bls12_generate_g1(&tmp_p);
  bls12_generate_g2(&tmp_q);
  efp12_scm(&tmp_q, &tmp_q, a);

  mpz_urandomm(s, state, order_z);

  efp12_scm(&A->p, &tmp_p, s);
  efp12_scm(&A->q, &tmp_q, s);

  mpz_clear(s);
}

void efp12_ecd(efp12_t *ANS, efp12_t *P) {
  static efp12_t tmp1_efp12;
  static fp12_t tmp1_fp12, tmp2_fp12, tmp3_fp12;
  if (fp12_cmp_zero(&P->y) == 0) {
    ANS->infinity = 1;
    return;
  }

  efp12_set(&tmp1_efp12, P);

  fp12_add(&tmp1_fp12, &tmp1_efp12.y, &tmp1_efp12.y);

  fp12_inv(&tmp1_fp12, &tmp1_fp12);
  fp12_sqr(&tmp2_fp12, &tmp1_efp12.x);
  fp12_add(&tmp3_fp12, &tmp2_fp12, &tmp2_fp12);
  fp12_add(&tmp2_fp12, &tmp2_fp12, &tmp3_fp12);
  fp12_mul(&tmp3_fp12, &tmp1_fp12, &tmp2_fp12);

  fp12_sqr(&tmp1_fp12, &tmp3_fp12);
  fp12_add(&tmp2_fp12, &tmp1_efp12.x, &tmp1_efp12.x);
  fp12_sub(&ANS->x, &tmp1_fp12, &tmp2_fp12);

  fp12_sub(&tmp1_fp12, &tmp1_efp12.x, &ANS->x);
  fp12_mul(&tmp2_fp12, &tmp3_fp12, &tmp1_fp12);
  fp12_sub(&ANS->y, &tmp2_fp12, &tmp1_efp12.y);
}

void efp12_ecd_lazy(efp12_t *ANS, efp12_t *P) {
  static efp12_t tmp1_efp12;
  static fp12_t tmp1_fp12, tmp2_fp12, tmp3_fp12;
  if (fp12_cmp_zero(&P->y) == 0) {
    ANS->infinity = 1;
    return;
  }

  efp12_set(&tmp1_efp12, P);

  fp12_add(&tmp1_fp12, &tmp1_efp12.y, &tmp1_efp12.y);

  fp12_inv(&tmp1_fp12, &tmp1_fp12);
  fp12_sqr_lazy(&tmp2_fp12, &tmp1_efp12.x);
  fp12_add_nonmod_single(&tmp3_fp12, &tmp2_fp12, &tmp2_fp12);
  fp12_add_nonmod_single(&tmp2_fp12, &tmp2_fp12, &tmp3_fp12);
  fp12_mul_lazy(&tmp3_fp12, &tmp1_fp12, &tmp2_fp12);

  fp12_sqr_lazy(&tmp1_fp12, &tmp3_fp12);
  fp12_add(&tmp2_fp12, &tmp1_efp12.x, &tmp1_efp12.x);
  fp12_sub(&ANS->x, &tmp1_fp12, &tmp2_fp12);

  fp12_sub_nonmod_single(&tmp1_fp12, &tmp1_efp12.x, &ANS->x);
  fp12_mul_lazy(&tmp2_fp12, &tmp3_fp12, &tmp1_fp12);
  fp12_sub(&ANS->y, &tmp2_fp12, &tmp1_efp12.y);
}

void efp12_eca(efp12_t *ANS, efp12_t *P1, efp12_t *P2) {
  static efp12_t tmp1_efp12, tmp2_efp12;
  static fp12_t tmp1_fp12, tmp2_fp12, tmp3_fp12;
  if (P1->infinity == 1) {
    efp12_set(ANS, P2);
    return;
  } else if (P2->infinity == 1) {
    efp12_set(ANS, P1);
    return;
  } else if (fp12_cmp(&P1->x, &P2->x) == 0) {
    if (fp12_cmp(&P1->y, &P2->y) != 0) {
      ANS->infinity = 1;
      return;
    } else {
      efp12_ecd(ANS, P1);
      return;
    }
  }

  efp12_set(&tmp1_efp12, P1);
  efp12_set(&tmp2_efp12, P2);

  fp12_sub(&tmp1_fp12, &tmp2_efp12.x, &tmp1_efp12.x);
  fp12_inv(&tmp1_fp12, &tmp1_fp12);
  fp12_sub(&tmp2_fp12, &tmp2_efp12.y, &tmp1_efp12.y);
  fp12_mul(&tmp3_fp12, &tmp1_fp12, &tmp2_fp12);
  fp12_sqr(&tmp1_fp12, &tmp3_fp12);
  fp12_sub(&tmp2_fp12, &tmp1_fp12, &tmp1_efp12.x);
  fp12_sub(&ANS->x, &tmp2_fp12, &tmp2_efp12.x);
  fp12_sub(&tmp1_fp12, &tmp1_efp12.x, &ANS->x);
  fp12_mul(&tmp2_fp12, &tmp3_fp12, &tmp1_fp12);
  fp12_sub(&ANS->y, &tmp2_fp12, &tmp1_efp12.y);
}
void efp12_eca_lazy(efp12_t *ANS, efp12_t *P1, efp12_t *P2) {
  static efp12_t tmp1_efp12, tmp2_efp12;
  static fp12_t tmp1_fp12, tmp2_fp12, tmp3_fp12;
  if (P1->infinity == 1) {
    efp12_set(ANS, P2);
    return;
  } else if (P2->infinity == 1) {
    efp12_set(ANS, P1);
    return;
  } else if (fp12_cmp(&P1->x, &P2->x) == 0) {
    if (fp12_cmp(&P1->y, &P2->y) != 0) {
      ANS->infinity = 1;
      return;
    } else {
      efp12_ecd_lazy(ANS, P1);
      return;
    }
  }

  efp12_set(&tmp1_efp12, P1);
  efp12_set(&tmp2_efp12, P2);

  fp12_sub(&tmp1_fp12, &tmp2_efp12.x, &tmp1_efp12.x);
  fp12_inv(&tmp1_fp12, &tmp1_fp12);
  fp12_sub_nonmod_single(&tmp2_fp12, &tmp2_efp12.y, &tmp1_efp12.y);
  fp12_mul_lazy(&tmp3_fp12, &tmp1_fp12, &tmp2_fp12);
  fp12_sqr_lazy(&tmp1_fp12, &tmp3_fp12);
  fp12_sub_nonmod_single(&tmp2_fp12, &tmp1_fp12, &tmp1_efp12.x);
  fp12_sub(&ANS->x, &tmp2_fp12, &tmp2_efp12.x);
  fp12_sub_nonmod_single(&tmp1_fp12, &tmp1_efp12.x, &ANS->x);
  fp12_mul_lazy(&tmp2_fp12, &tmp3_fp12, &tmp1_fp12);
  fp12_sub(&ANS->y, &tmp2_fp12, &tmp1_efp12.y);
}

void efp12_scm(efp12_t *ANS, efp12_t *P, mpz_t scalar) {
  if (mpz_cmp_ui(scalar, 0) == 0) {
    ANS->infinity = 1;
    return;
  } else if (mpz_cmp_ui(scalar, 1) == 0) {
    efp12_set(ANS, P);
    return;
  }

  efp12_t Tmp_P, Next_P;
  efp12_init(&Tmp_P);
  efp12_set(&Tmp_P, P);
  efp12_init(&Next_P);
  int i, length;
  length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);

  efp12_set(&Next_P, &Tmp_P);
  for (i = 1; i < length; i++) {
    efp12_ecd(&Next_P, &Next_P);
    if (binary[i] == '1') {
      efp12_eca(&Next_P, &Next_P, &Tmp_P);
    }
  }
  efp12_set(ANS, &Next_P);
}

void efp12_scm_lazy(efp12_t *ANS, efp12_t *P, mpz_t scalar) {
  if (mpz_cmp_ui(scalar, 0) == 0) {
    ANS->infinity = 1;
    return;
  } else if (mpz_cmp_ui(scalar, 1) == 0) {
    efp12_set(ANS, P);
    return;
  }

  efp12_t Tmp_P, Next_P;
  efp12_init(&Tmp_P);
  efp12_set(&Tmp_P, P);
  efp12_init(&Next_P);
  int i, length;
  length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);

  efp12_set(&Next_P, &Tmp_P);
  for (i = 1; i < length; i++) {
    efp12_ecd_lazy(&Next_P, &Next_P);
    if (binary[i] == '1') {
      efp12_eca_lazy(&Next_P, &Next_P, &Tmp_P);
    }
  }
  efp12_set(ANS, &Next_P);
}

void bls12_generate_g1_fast(efp12_t *P) {
  efp_t tmp_P;
  efp_init(&tmp_P);
  mpz_t exp;
  mpz_init(exp);

  efp_set_random_fast(&tmp_P, state);
  efp12_set_ui(P, 0, 0);
  mpz_tdiv_q(exp, efp_total, order_z);
  efp_scm(&tmp_P, &tmp_P, exp);
  fp_set(&P->x.x0.x0.x0, &tmp_P.x);
  fp_set(&P->y.x0.x0.x0, &tmp_P.y);
  P->infinity = tmp_P.infinity;

  mpz_clear(exp);
}