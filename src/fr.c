#include <ELiPS/fr.h>
/**************set up*******************/
void mpn_set_mpz_size(mp_limb_t *ans, mpz_t a, mp_size_t size) {
  char *str;

  str = (char *)malloc(mpz_sizeinbase(a, 10) + 2);

  //gmp_printf("a=%Zd\n",a);
  str = mpz_get_str(str, 10, a);
  //printf("str1=%s",str);
  mpn_set_char(ans, size, str);

  free(str);
}

void mpz_set_mpn_size(mpz_t ans, mp_limb_t *a, int bits, mp_size_t size) {
  unsigned char str[bits / 4];  //16進数表記(4bit)でするので、割る4している
  mp_size_t str_size;
  mp_size_t a_size = -1;
  for (int i = size - 1; i >= 0; i--) {
    if (a[i] == 0) continue;
    a_size = i;
    break;
  }
  a_size++;
  str_size = mpn_get_str(str, 16, a, a_size);
  mpz_set_ui(ans, 0);
  for (int i = 0; i < str_size; i++) {
    mpz_mul_ui(ans, ans, 16);
    mpz_add_ui(ans, ans, str[i]);
  }
}

void fr_order_init() {
  mpn_set_mpz_size(order, order_z, FRLIMB);
#ifdef X_PLUS
  mpn_set_mpz_size(X_abs, X_z, FXLIMB);
#endif
#ifdef X_MINUS
  mpz_neg(X_z, X_z);
  mpn_set_mpz_size(X_abs, X_z, FXLIMB);
  mpz_neg(X_z, X_z);
#endif
  mp_limb_t tmp[FXLIMB * 2];  //TODO
  mpn_mul(tmp, X_abs, FXLIMB, X_abs, FXLIMB);
  mpn_copyd(X2, tmp, FXLIMB2);
}

void to_g1_expo_init() {
  mpz_init(to_g1_expo);
  mpz_tdiv_q(to_g1_expo, efp_total, order_z);
}

void to_g2_expo_init() {
  mpz_init(to_g2_expo);

  // (x^8 -4x^7 +5x^6 -4x^4 +6x^3 -4x^2 -4x +13) / 9
  mpz_t tmp;
  mpz_init(tmp);

  mpz_pow_ui(tmp, X_z, 8);
  mpz_set(to_g2_expo, tmp);
  mpz_pow_ui(tmp, X_z, 7);
  mpz_mul_ui(tmp, tmp, 4);
  mpz_sub(to_g2_expo, to_g2_expo, tmp);
  mpz_pow_ui(tmp, X_z, 6);
  mpz_mul_ui(tmp, tmp, 5);
  mpz_add(to_g2_expo, to_g2_expo, tmp);
  mpz_pow_ui(tmp, X_z, 4);
  mpz_mul_ui(tmp, tmp, 4);
  mpz_sub(to_g2_expo, to_g2_expo, tmp);
  mpz_pow_ui(tmp, X_z, 3);
  mpz_mul_ui(tmp, tmp, 6);
  mpz_add(to_g2_expo, to_g2_expo, tmp);
  mpz_pow_ui(tmp, X_z, 2);
  mpz_mul_ui(tmp, tmp, 4);
  mpz_sub(to_g2_expo, to_g2_expo, tmp);
  mpz_mul_ui(tmp, X_z, 4);
  mpz_sub(to_g2_expo, to_g2_expo, tmp);
  mpz_add_ui(to_g2_expo, to_g2_expo, 13);

  if (mpz_divisible_ui_p(to_g2_expo, 9) == 0) {
    printf("error to_g2_expo_init()\n");
    fflush(stdout);
  }
  mpz_divexact_ui(to_g2_expo, to_g2_expo, 9);

  mpz_clear(tmp);
}

void curve_b_montgomery_init() {
  fp_init(&curve_b_montgomery);
  fp_set_mpn(&curve_b_montgomery, curve_b);
  fp_to_montgomery(&curve_b_montgomery, &curve_b_montgomery);
}

void twist_curve_b_montgomery_init() {
  fp2_init(&twist_curve_b_montgomery);
  fp_set_mpn(&twist_curve_b_montgomery.x0, curve_b);
  fp_set_mpn(&twist_curve_b_montgomery.x1, curve_b);
  fp2_to_montgomery(&twist_curve_b_montgomery, &twist_curve_b_montgomery);
}

void inv2_montgomery_init() {
  fp_init(&inv2_montgomery);
  fp_set_ui(&inv2_montgomery, 2);
  fp_inv(&inv2_montgomery, &inv2_montgomery);
  fp_to_montgomery(&inv2_montgomery, &inv2_montgomery);
}
int w_naf_frt(int *dw, mpz_t d, int w) {
  int i = 0;
  mp_limb_t tmp_d[FRLIMB];
  int n_int;

  mpn_set_mpz_size(tmp_d, d, FRLIMB);

  if (w == 5) {
    n_int = 32;
    while (mpn_zero_p(tmp_d, FRLIMB) != 1) {
      if (tmp_d[0] % 2 == 1) {
        dw[i] = tmp_d[0] % n_int;
        if (dw[i] <= 16) {
          mpn_sub_ui(tmp_d, tmp_d, FRLIMB, dw[i]);  //
        } else if (mpn_zero_p(tmp_d, FRLIMB) != 1) {
          dw[i] = dw[i] - n_int;
          mpn_add_ui(tmp_d, tmp_d, FRLIMB, -dw[i]);  //
        }
      } else
        dw[i] = 0;
      mpn_rshift(tmp_d, tmp_d, FRLIMB, 1);  //
      i++;
    }
  }
  return i - 1;
}
// int w_naf_frt_mpz(int *dw,mp_limb_t *d,mp_size_t size,int w){
//     mpz_t tmp;
//     mpz_init(&tmp);
//     mpz_set_mpn_size(tmp,d,,size);
//     mpn_set_mpz_size(dw,tmp,w);
//     mpz_clear(&tmp);
// }

// int w_naf_frt_mpn(int *dw,mp_limb_t d[FXLIMB2],int w){
// 	int i=0;
// 	//mpz_t dw_t,buf,n;
// 	mp_limb_t tmp_d[FXLIMB2];
//     mpn_copyd(tmp_d,d,FXLIMB2);
//     //mpz_init(dw_t);
//     //mpz_init(buf);
//     //mpz_init(n);
//     int n_int;

//     //mpn_set_mpz_size(tmp_d,d,FXLIMB);

// 	if(w==5){
//         n_int=32;
//         //mpz_set_ui(n,32);
//         //while(mpz_cmp_ui(d,0)!=0){
//         while(mpn_zero_p(tmp_d,FXLIMB2)!=1){
//             //if(mpz_odd_p(d)!=0){
//             if(tmp_d[0]%2==1){
//             //if(mpn_scan1(tmp_d,0)==0){
//                 //mpz_mod(dw_t,d,n);
//                 dw[i]=tmp_d[0]%n_int;
//                 //if(mpz_cmp_ui(dw_t,16)<=0)	{
//                 if(dw[i]<=16)	{
//                     //dw[i]=mpz_get_ui(dw_t);
//                     //mpz_sub_ui(d,d,dw[i]);
//                     mpn_sub_ui(tmp_d,tmp_d,FXLIMB2,dw[i]);//
//                 }
//                 //else if(mpz_cmp_ui(d,0)!=0) {
//                 else if(mpn_zero_p(tmp_d,FXLIMB2)!=1) {
//                     //mpz_sub(buf,n,dw_t);
//                     //dw[i]=-(32-dw[i]);
//                     dw[i]=dw[i]-n_int;
//                     //dw[i]=-mpz_get_ui(buf);
//                     //mpz_add(d,d,buf);
//                     //mpz_add_ui(d,d,-dw[i]);
//                     //dw[i]=32-dw[i];
//                     mpn_add_ui(tmp_d,tmp_d,FXLIMB2,-dw[i]);//
//                 }
//             }else dw[i]=0;
//             // gmp_printf("A=%Nu",tmp_d,FRLIMB);
//             // gmp_printf("A=%Zd",d);
//             // printf("\n");
//             // getchar();
//             //mpz_tdiv_q_2exp(d,d,1);
//             mpn_rshift(tmp_d,tmp_d,FXLIMB2,1);//
//             //gmp_printf("A=%Nu",tmp_d,FRLIMB);
//             //gmp_printf("A=%Zd",d);
//             //getchar();
//             i++;
//         }
//     }
// 	return i-1;
// }
int w_naf_frt_mpn_size(int *dw, mp_limb_t *d, mp_size_t size, int w) {
  int i = 0;
  mp_limb_t tmp_d[size];
  mpn_copyd(tmp_d, d, size);
  int n_int;
  if (w == 5) {
    n_int = 32;
    while (mpn_zero_p(tmp_d, size) != 1) {
      if (tmp_d[0] % 2 == 1) {
        dw[i] = tmp_d[0] % n_int;
        if (dw[i] <= 16)
          mpn_sub_ui(tmp_d, tmp_d, size, dw[i]);
        else if (mpn_zero_p(tmp_d, size) != 1) {
          dw[i] = dw[i] - n_int;
          mpn_add_ui(tmp_d, tmp_d, size, -dw[i]);  //
        }
      } else
        dw[i] = 0;
      mpn_rshift(tmp_d, tmp_d, size, 1);
      i++;
    }
  }
  // printf("2: ok! %d\n",i);
  // mpn_printf("d=",d,size);
  // //getchar();
  // debug_get("");
  return i - 1;
}
/************FR_T**********************/
void fr_init(fr_t *A) {
  mpn_zero(A->x0, FRLIMB);
}

void mpz_set_fr(mpz_t ans, fr_t *a) {
  unsigned char str[FRLIMB_BITS / 4];  //16進数表記(4bit)でするので、割る4している
  mp_size_t str_size;
  mp_size_t a_size = -1;
  for (int i = FRLIMB - 1; i >= 0; i--) {
    if (a->x0[i] == 0) continue;
    a_size = i;
    break;
  }
  a_size++;
  str_size = mpn_get_str(str, 16, a->x0, a_size);
  mpz_set_ui(ans, 0);
  for (int i = 0; i < str_size; i++) {
    mpz_mul_ui(ans, ans, 16);
    mpz_add_ui(ans, ans, str[i]);
  }
}

size_t fr_sizeinbase_2(fr_t *A) {
  return mpn_sizeinbase(A->x0, FRLIMB, 2);
}
void fr_printf(char *str, fr_t *A) {
  gmp_printf("%s%Nu", str, A->x0, FRLIMB);
}

void fr_println(char *str, fr_t *A) {
  gmp_printf("%s%Nu\n", str, A->x0, FRLIMB);
}

void fr_set(fr_t *ANS, fr_t *A) {
  mpn_copyd(ANS->x0, A->x0, FRLIMB);
}

void fr_set_ui(fr_t *ANS, unsigned long int UI) {
  mpn_set_ui(ANS->x0, FRLIMB, UI);
}

void fr_set_mpn(fr_t *ANS, mp_limb_t *A) {
  mpn_copyd(ANS->x0, A, FRLIMB);
}

void fr_mod(fr_t *ans, mp_limb_t *a, mp_size_t size_a) {
  mp_limb_t dumy[size_a];
  mpn_tdiv_qr(dumy, ans->x0, 0, a, size_a, order, FRLIMB);
}

void fr_add(fr_t *ANS, fr_t *A, fr_t *B) {
  mpn_add_n(ANS->x0, A->x0, B->x0, FRLIMB);
  if (mpn_cmp(ANS->x0, order, FRLIMB) >= 0) mpn_sub_n(ANS->x0, ANS->x0, order, FRLIMB);
}

void fr_neg(fr_t *ANS, fr_t *A) {
  if (mpn_cmp_ui(A->x0, FRLIMB, 0) == 0)
    fr_set(ANS, A);
  else
    mpn_sub_n(ANS->x0, order, A->x0, FPLIMB);
}

void fr_mul(fr_t *ANS, fr_t *A, fr_t *B) {
  static mp_limb_t tmp_mul[FRLIMB2];
  mpn_mul_n(tmp_mul, A->x0, B->x0, FRLIMB);
  fr_mod(ANS, tmp_mul, FRLIMB2);
}

void fr_inv(fr_t *ANS, fr_t *A) {
  static mp_limb_t order_tmp[FRLIMB], gp[FRLIMB], sp[FRLIMB], buf[FRLIMB];
  static mp_size_t buf_size;

  mpn_init(sp, FRLIMB);

  mpn_copyd(order_tmp, order, FRLIMB);

  mpn_add_n(buf, A->x0, order, FRLIMB);
  mpn_gcdext(gp, sp, &buf_size, buf, FRLIMB, order_tmp, FRLIMB);

  if (buf_size < 0) {
    mpn_sub_n(ANS->x0, order, sp, FRLIMB);
  } else {
    mpn_copyd(ANS->x0, sp, FRLIMB);
  }
}
void fr_pow(fr_t *ANS, fr_t *A, mpz_t scalar) {
  int i, length;
  length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);
  fr_t tmp;
  fr_init(&tmp);

  fr_set(&tmp, A);

  for (i = 1; i < length; i++) {
    fr_mul(&tmp, &tmp, &tmp);
    if (binary[i] == '1') {
      fr_mul(&tmp, A, &tmp);
    }
  }
  fr_set(ANS, &tmp);
}
void fr_sub(fr_t *ANS, fr_t *A, fr_t *B) {
  static mp_limb_t buf[FRLIMB];
  if (mpn_cmp(A->x0, B->x0, FRLIMB) < 0) {
    mpn_sub_n(buf, A->x0, B->x0, FRLIMB);
    mpn_add_n(ANS->x0, order, buf, FRLIMB);
  } else {
    mpn_sub_n(ANS->x0, A->x0, B->x0, FRLIMB);
  }
}

void fr_set_random(fr_t *ANS, gmp_randstate_t state) {
  mpz_t tmp;
  mpz_init(tmp);
  mpz_urandomm(tmp, state, order_z);
  mpn_set_mpz_size(ANS->x0, tmp, FRLIMB);
  mpz_clear(tmp);
}

/************g1_t**************/
void g1_init(g1_t *A) {
  efp_init(A);
  //A->infinity=1;
}

void g1_printf(char *s, g1_t *A) {
  g1_t tmp;
  efp_mod_montgomery(&tmp, A);
  efp_printf(s, &tmp);
}

void g1_println(char *s, g1_t *A) {
  g1_t tmp;
  efp_mod_montgomery(&tmp, A);
  efp_println(s, &tmp);
}

void g1_set(g1_t *ANS, g1_t *A) {
  efp_set(ANS, A);
}

int g1_cmp(g1_t *ANS, g1_t *A) {
  return efp_cmp(ANS, A);
}

void g1_set_random(g1_t *ANS, gmp_randstate_t state) {
#if 0
  efp12_t P;
  g1_t P_twisted;
  //bls12_generate_g1(&P);
  bls12_generate_g1_fast(&P);

  efp12_to_efp(&P_twisted, &P);
  efp_to_montgomery(ANS, &P_twisted);
  ANS->infinity = 0;
#else
  efp_t P;
  efp_init(&P);
  P.infinity = 0;
  fp_set_random(&P.x, state);
  fp_to_montgomery(&P.x, &P.x);

  fp_t tmp1, tmp2, tmp_x;
  fp_init(&tmp1);
  fp_init(&tmp2);
  fp_init(&tmp_x);
  while (1) {
    fp_sqrmod_montgomery(&tmp1, &P.x);
    fp_mulmod_montgomery(&tmp2, &tmp1, &P.x);
    fp_add(&tmp_x, &tmp2, &curve_b_montgomery);
    if (fp_legendre_sqrt_montgomery(&P.y, &tmp_x) != -1) {
      break;
    }
    fp_add_mpn(&P.x, &P.x, RmodP);
  }

  efp_jacobian_t P_jacobian;
  efp_affine_to_jacobian_montgomery(&P_jacobian, &P);
  efp_scm_jacobian_lazy_montgomery(&P_jacobian, &P_jacobian, to_g1_expo);
  efp_jacobian_to_affine_montgomery(ANS, &P_jacobian);
#endif
}

//check function
int g1_cmp_efp12(g1_t *A, efp12_t *B) {
  efp12_t A_efp12;
  efp_t A_non_monty;
  efp_init(&A_non_monty);
  efp12_init(&A_efp12);
  efp_mod_montgomery(&A_non_monty, A);
  efp_to_efp12(&A_efp12, &A_non_monty);
  return efp12_cmp(&A_efp12, B);
}

void g1_ecd(g1_t *ANS, g1_t *P) {
  efp_ecd_lazy_montgomery(ANS, P);
}

void g1_eca(g1_t *ANS, g1_t *P, g1_t *Q) {
  efp_eca_lazy_montgomery(ANS, P, Q);
}

void g1_neg(g1_t *ANS, g1_t *P) {
  fp_set(&ANS->x, &P->x);
  fp_set_neg(&ANS->y, &P->y);
  ANS->infinity = P->infinity;
}

void g1_scm(g1_t *ANS, g1_t *P, fr_t *sca) {
  if (P->infinity == 1) {
    g1_set(ANS, P);
  } else {
    //s=s0+s1[x^4]
    int i, length_s[2], loop_length;
    efp_t next_tmp_P, tmp_P;
    efp_jacobian_t next_tmpJ_P, tmpJ_P[8], tmpJ_P_neg[8], tmpJ_P_4x[8], tmpJ_P_4x_neg[8], tmpJ_2P;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_jacobian_init(&tmpJ_2P);
    for (i = 0; i < 8; i++) {
      efp_jacobian_init(&tmpJ_P[i]);
      efp_jacobian_init(&tmpJ_P_neg[i]);
      efp_jacobian_init(&tmpJ_P_4x[i]);
      efp_jacobian_init(&tmpJ_P_4x_neg[i]);
    }
    efp_jacobian_init(&next_tmpJ_P);

    mp_limb_t s1[FRLIMB - FXLIMB2 + 1], s0[FXLIMB2];
    //table
    efp_jacobian_t table0[17], table1[17];
    for (i = 0; i < 17; i++) {
      efp_jacobian_init(&table0[i]);
      efp_jacobian_init(&table1[i]);
    }

    //set
    efp_set(&tmp_P, P);
    efp_affine_to_jacobian_montgomery(&tmpJ_P[0], &tmp_P);
    efp_ecd_jacobian_lazy_montgomery(&tmpJ_2P, &tmpJ_P[0]);
    for (i = 1; i < 8; i++) {
      efp_eca_jacobian_lazy_montgomery(&tmpJ_P[i], &tmpJ_P[i - 1], &tmpJ_2P);
    }

    fp_t point_table[8], inv_table[8];
    for (i = 0; i < 8; i++) fp_set(&point_table[i], &tmpJ_P[i].z);
    fp_montgomery_trick_montgomery(inv_table, point_table, 8);
    for (i = 0; i < 8; i++) efp_jacobian_to_mixture_noninv_montgomery(&tmpJ_P[i], &tmpJ_P[i], &inv_table[i]);

    for (i = 0; i < 8; i++) {
      efp_jacobian_set_neg(&tmpJ_P_neg[i], &tmpJ_P[i]);                          //tmp_P_neg
      efp_jacobian_skew_frobenius_map_p2_montgomery(&tmpJ_P_4x[i], &tmpJ_P[i]);  //tmp_P_4x
      efp_jacobian_set_neg(&tmpJ_P_4x_neg[i], &tmpJ_P_4x[i]);                    //tmp_P_4x_neg
    }
    //set table
    table0[0].infinity = 1;  //0
    table1[0].infinity = 1;  //0

    for (i = 0; i < 8; i++) {
      efp_jacobian_set(&table0[i + 1], &tmpJ_P[i]);         //[1]P
      efp_jacobian_set(&table0[i + 9], &tmpJ_P_neg[i]);     //[-1]P
      efp_jacobian_set(&table1[i + 1], &tmpJ_P_4x[i]);      //[1]P'
      efp_jacobian_set(&table1[i + 9], &tmpJ_P_4x_neg[i]);  //[-1]P'
    }
    mpn_tdiv_qr(s1, s0, 0, sca->x0, FRLIMB, X2, FXLIMB2);

    //get loop_length
    loop_length = 0;
    length_s[0] = (int)mpn_sizeinbase(s0, FXLIMB2, 2);               //mpn
    length_s[1] = (int)mpn_sizeinbase(s1, FRLIMB - FXLIMB2 + 1, 2);  //mpn
    if (length_s[0] > length_s[1]) {
      loop_length = length_s[0];
    } else {
      loop_length = length_s[1];
    }

    //naf
    int naf_length, naf_length0, naf_length1;
    int naf_binary[2][loop_length + 1];
    for (i = 0; i < loop_length + 1; i++) {
      naf_binary[0][i] = 0;
      naf_binary[1][i] = 0;
    }
    int *naf_pointer[2];
    naf_pointer[0] = naf_binary[0];
    naf_pointer[1] = naf_binary[1];

    naf_length0 = w_naf_frt_mpn_size(naf_binary[0], s0, FXLIMB2, 5);
    naf_length1 = w_naf_frt_mpn_size(naf_binary[1], s1, FRLIMB - FXLIMB2 + 1, 5);
    if (naf_length0 < naf_length1)
      naf_length = naf_length1;
    else
      naf_length = naf_length0;

    //naf_length=loop_length-1;
    int binary0[naf_length + 1], binary1[naf_length + 1];

    for (i = naf_length; i >= 0; i--) {
      if (naf_binary[0][i] == 0)
        binary0[i] = 0;
      else if (naf_binary[0][i] > 0)
        binary0[i] = (naf_binary[0][i] + 1) >> 1;
      else
        binary0[i] = ((17 - (naf_binary[0][i] + 16)) >> 1) + 8;

      if (naf_binary[1][i] == 0)
        binary1[i] = 0;
      else if (naf_binary[1][i] > 0)
        binary1[i] = (naf_binary[1][i] + 1) >> 1;
      else
        binary1[i] = ((17 - (naf_binary[1][i] + 16)) >> 1) + 8;
    }
    if (naf_length0 == naf_length1)
      efp_eca_jacobian_lazy_montgomery(&next_tmpJ_P, &table0[binary0[naf_length]], &table1[binary1[naf_length]]);
    else if (naf_length0 < naf_length1)
      efp_jacobian_set(&next_tmpJ_P, &table1[binary1[naf_length]]);
    else
      efp_jacobian_set(&next_tmpJ_P, &table0[binary0[naf_length]]);

    //scm
    for (i = naf_length - 1; i >= 0; i--) {
      efp_ecd_jacobian_lazy_montgomery(&next_tmpJ_P, &next_tmpJ_P);
      if (binary0[i] != 0) efp_eca_mixture_lazy_montgomery_ignore_inf(&next_tmpJ_P, &next_tmpJ_P, &table0[binary0[i]]);
      if (binary1[i] != 0) efp_eca_mixture_lazy_montgomery_ignore_inf(&next_tmpJ_P, &next_tmpJ_P, &table1[binary1[i]]);
    }

    efp_jacobian_to_affine_montgomery(&next_tmp_P, &next_tmpJ_P);
    efp_set(ANS, &next_tmp_P);
    ANS->infinity = next_tmp_P.infinity;
  }
}

void g1_set_random_with_basepoint(g1_t *ANS, g1_t *basepoint, gmp_randstate_t state) {
  fr_t s;
  fr_set_random(&s, state);
  g1_scm(ANS, basepoint, &s);
  ANS->infinity = 0;
}

void g1_test(int scm) {
  g1_t A_g1;
  efp12_t A_efp12;
  fr_t sca_fr;
  mpz_t sca;
  //initialize
  mpz_init(sca);
  fr_init(&sca_fr);
  g1_init(&A_g1);
  efp12_init(&A_efp12);
  //set
  bls12_generate_g1(&A_efp12);
  efp12_to_efp(&A_g1, &A_efp12);
  efp_to_montgomery(&A_g1, &A_g1);
  fr_set_random(&sca_fr, state);
  //fr_set_ui(&sca_fr,4);
  //mpz_set_ui(sca,10);
  mpz_set_fr(sca, &sca_fr);
  //calc
  printf("\n===========================\n");
  g1_scm(&A_g1, &A_g1, &sca_fr);
  //g1_scm_mpz(&A_g1,&A_g1,sca);
  //efp12_scm(&A_efp12,&A_efp12,order_z);
  efp12_scm(&A_efp12, &A_efp12, sca);
  efp12_printf("A=", &A_efp12);
  printf("\n\n");
  g1_printf("A=", &A_g1);
  printf("\n");
  printf("g1 test:");
  if (g1_cmp_efp12(&A_g1, &A_efp12) == 0)
    printf("ok!\n");
  else
    printf("ng\n");
  //bench
  int i = 0;
  float scm_time = 0;
  cost tmp, scm_cost;
  struct timeval tv_A, tv_B;
  for (i = 0; i < scm; i++) {
    fr_set_random(&sca_fr, state);
    //faster type
    cost_zero();
    gettimeofday(&tv_A, NULL);
    g1_scm(&A_g1, &A_g1, &sca_fr);
    gettimeofday(&tv_B, NULL);
    scm_time += timedifference_msec(tv_A, tv_B);
    cost_check(&tmp);
    cost_addition(&scm_cost, &tmp);
  }
  printf("bls12 g1 scm.     : %.4f[ms]\n", scm_time / scm);

#ifdef DEBUG_COST_A
  printf("*********bls12 g1 scm fp COST.********         \n");
  cost_printf("bls12 g1 scm", &scm_cost, scm);
  printf("***************************************         \n");
#endif
  mpz_clear(sca);
}

void g1_set_random_schoolbook(g1_t *ANS, gmp_randstate_t state) {
  efp12_t P;
  g1_t P_twisted;
  bls12_generate_g1(&P);
  efp12_to_efp(&P_twisted, &P);
  efp_to_montgomery(ANS, &P_twisted);
  ANS->infinity = 0;
}

void g1_set_random_test(int scm) {
  g1_t P_test;
  int i;
  float random_time = 0, random_fast_time = 0;
  struct timeval tv_A, tv_B;
  for (i = 0; i < scm; i++) {
    //normal type
    gettimeofday(&tv_A, NULL);
    g1_set_random_schoolbook(&P_test, state);
    gettimeofday(&tv_B, NULL);
    random_time += timedifference_msec(tv_A, tv_B);

    //faster type
    gettimeofday(&tv_A, NULL);
    g1_set_random(&P_test, state);
    gettimeofday(&tv_B, NULL);
    random_fast_time += timedifference_msec(tv_A, tv_B);
  }
  printf("g1_set_random_schoolbook.     : %.4f[ms]\n", random_time / scm);
  printf("g1_set_random.          : %.4f[ms]\n", random_fast_time / scm);
}

/************g2_t**************/
void g2_init(g2_t *A) {
  efp2_init(A);
  A->infinity = 1;
}

void g2_printf(char *s, g2_t *A) {
  g2_t tmp;
  efp2_mod_montgomery(&tmp, A);
  efp2_printf(s, &tmp);
}

void g2_println(char *s, g2_t *A) {
  g2_t tmp;
  efp2_mod_montgomery(&tmp, A);
  efp2_println(s, &tmp);
}

void g2_set(g2_t *ANS, g2_t *A) {
  efp2_set(ANS, A);
}

int g2_cmp(g2_t *ANS, g2_t *A) {
  return efp2_cmp(ANS, A);
}

void g2_set_random(g2_t *ANS, gmp_randstate_t state) {
#if 0
  efp12_t Q;
  g2_t Q_twisted;
  bls12_generate_g2(&Q);
  efp12_to_efp2(&Q_twisted, &Q);
  efp2_to_montgomery(ANS, &Q_twisted);
  ANS->infinity = 0;
#else
  efp2_t P;
  efp2_init(&P);
  P.infinity = 0;
  fp2_set_random(&P.x, state);
  fp2_to_montgomery(&P.x, &P.x);

  fp2_t tmp1, tmp2, tmp_x;
  fp2_init(&tmp1);
  fp2_init(&tmp2);
  fp2_init(&tmp_x);
  while (1) {
    fp2_sqr_lazy_montgomery(&tmp1, &P.x);
    fp2_mul_lazy_montgomery(&tmp2, &tmp1, &P.x);
#ifdef TWIST_PHI
    fp2_add(&tmp_x, &tmp2, &twist_curve_b_montgomery);
#endif
#ifdef TWIST_PHI_INV
    printf("TODO: peks_hash1_g2(), EP_TYPE2\n");
#endif
    if (fp2_sqrt_complex_method_montgomery(&P.y, &tmp_x) == 1) {
      break;
    }
    fp_add_mpn(&P.x.x0, &P.x.x0, RmodP);
  }

  map_to_g2_montgomery(ANS, &P);
#endif
}

//check function
int g2_cmp_efp12(g2_t *A, efp12_t *B) {
  efp12_t A_efp12;
  efp2_t A_non_monty;
  efp2_init(&A_non_monty);
  efp12_init(&A_efp12);
  efp2_mod_montgomery(&A_non_monty, A);
  efp2_to_efp12(&A_efp12, &A_non_monty);
  return efp12_cmp(&A_efp12, B);
}
void g2_ecd(g2_t *ANS, g2_t *Q) {
  efp2_ecd_lazy_montgomery(ANS, Q);
}
void g2_eca(g2_t *ANS, g2_t *P, g2_t *Q) {
  efp2_eca_lazy_montgomery(ANS, P, Q);
}
void g2_neg(g2_t *ANS, g2_t *P) {
  fp2_set(&ANS->x, &P->x);
  fp2_set_neg(&ANS->y, &P->y);
  ANS->infinity = P->infinity;
}
void map_to_g2(g2_t *ANS, efp2_t *A) {
  // (χ<0, z = α + 1)
  // It is diffrent output between Scott et al. method and Fuentes et al. method. But both method is map to g2.
  // Fuentes et al. method is faster than the other

#if 0  // Scott et al. method
  efp2_t A0, A1, A2;
  efp2_t tmp1, tmp2;
  efp2_init(&A0);
  efp2_init(&A1);
  efp2_init(&A2);
  efp2_init(&tmp1);
  efp2_init(&tmp2);

  efp2_set(&A0, A);                      // A0 = A
  efp2_skew_frobenius_map_p1(&A1, &A0);  // A1 = ψ(A)
  efp2_skew_frobenius_map_p2(&A2, &A0);  // A2 = ψ^2(A)

  efp2_eca(&tmp1, &A0, &A1);
  efp2_set_neg(&tmp1, &tmp1);  // tmp1 = -(A + ψ(A))
  efp2_scm(ANS, &tmp1, X_abs_z);

  efp2_set_neg(&tmp2, &A2);
  efp2_eca(&tmp2, &tmp1, &tmp2);  // tmp2 = -(A + ψ(A) + ψ^2(A))
  efp2_eca(ANS, ANS, &tmp2);
  efp2_scm(ANS, ANS, X_abs_z);
  efp2_set_neg(ANS, ANS);

  efp2_eca(&tmp2, &tmp1, &A2);  // tmp2 = -(A + ψ(A)) + ψ^2(A)
  efp2_eca(&tmp1, &tmp2, &A2);  // tmp1 = -(A + ψ(A)) + [2]ψ^2(A)
  efp2_eca(ANS, ANS, &tmp1);
  efp2_scm(ANS, ANS, X_abs_z);
  efp2_set_neg(ANS, ANS);

  efp2_ecd(&tmp1, &A0);
  efp2_eca(&tmp1, &tmp1, &A0);    // tmp1 = [3](A)
  efp2_set_neg(&tmp2, &tmp2);     // tmp2 = A + ψ(A) - ψ^2(A)
  efp2_eca(&tmp2, &tmp2, &tmp1);  // tmp2 = [4]A + ψ(A) - ψ^2(A)
  efp2_eca(ANS, ANS, &tmp2);
#else  // Fuentes et al. method
  efp2_t A0, A1, A2;
  efp2_init(&A0);
  efp2_init(&A1);
  efp2_init(&A2);

  efp2_set_neg(&A0, A);                // A0 = -A
  efp2_skew_frobenius_map_p1(&A1, A);  // A1 = ψ(A)
  efp2_ecd(&A2, A);
  efp2_skew_frobenius_map_p2(&A2, &A2);  // A2 = ψ^2(2A)

  efp2_scm(ANS, &A0, X_abs_z);

  efp2_eca(ANS, ANS, &A0);
  efp2_eca(ANS, ANS, &A1);
  efp2_scm(ANS, ANS, X_abs_z);
  efp2_set_neg(ANS, ANS);

  efp2_set_neg(&A1, &A1);  // A1 = -ψ(A)
  efp2_eca(ANS, ANS, &A0);
  efp2_eca(ANS, ANS, &A1);
  efp2_eca(ANS, ANS, &A2);
#endif
  efp2_to_montgomery(ANS, ANS);
}
void map_to_g2_montgomery(g2_t *ANS, efp2_t *A) {
  // (χ<0, z = α + 1)
  // It is diffrent output between Scott et al. method and Fuentes et al. method. But both method is map to g2.

  // Fuentes et al. method
  efp2_t tmp, ans;
  efp2_init(&tmp);
  efp2_init(&ans);
  efp2_jacobian_t A_jacobian, ANS_jacobian, A0, A1, A2;
  efp2_jacobian_init(&A_jacobian);
  efp2_jacobian_init(&A0);
  efp2_jacobian_init(&A1);
  efp2_jacobian_init(&A2);

  efp2_affine_to_jacobian_montgomery(&A_jacobian, A);
  efp2_jacobian_set_neg(&A0, &A_jacobian);                           // A0 = -A
  efp2_jacobian_skew_frobenius_map_p1_montgomery(&A1, &A_jacobian);  // A1 = ψ(A)
  efp2_ecd_jacobian_lazy_montgomery(&A2, &A_jacobian);
  efp2_jacobian_skew_frobenius_map_p2_montgomery(&A2, &A2);  // A2 = ψ^2(2A)

  efp2_scm_X_jacobian_lazy_montgomery(&ANS_jacobian, &A0);
  efp2_jacobian_set_neg(&ANS_jacobian, &ANS_jacobian);

  efp2_eca_jacobian_lazy_montgomery(&ANS_jacobian, &ANS_jacobian, &A0);
  efp2_eca_jacobian_lazy_montgomery(&ANS_jacobian, &ANS_jacobian, &A1);
  efp2_scm_X_jacobian_lazy_montgomery(&ANS_jacobian, &ANS_jacobian);

  efp2_jacobian_set_neg(&A1, &A1);  // A1 = -ψ(A)
  efp2_eca_jacobian_lazy_montgomery(&ANS_jacobian, &ANS_jacobian, &A0);
  efp2_eca_jacobian_lazy_montgomery(&ANS_jacobian, &ANS_jacobian, &A1);
  efp2_eca_jacobian_lazy_montgomery(&ANS_jacobian, &ANS_jacobian, &A2);

  efp2_jacobian_to_affine_montgomery(ANS, &ANS_jacobian);
}


void g2_scm(g2_t *ANS, g2_t *Q, fr_t *sca) {
#if ARCBIT == 64
  if (Q->infinity == 1) {
    g2_set(ANS, Q);
  } else {
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i, j, length_s[4], loop_length;
    efp2_t next_twisted_Q, twisted_Q, twisted_2Q;
    efp2_jacobian_t next_twistedJ_Q, twistedJ_2Q;
    efp2_jacobian_t twistedJ_Q[8], twistedJ_Q_x[8], twistedJ_Q_2x[8], twistedJ_Q_3x[8];
    efp2_jacobian_t twistedJ_Q_neg[8], twistedJ_Q_x_neg[8], twistedJ_Q_2x_neg[8], twistedJ_Q_3x_neg[8];

    efp2_init(&twisted_Q);
    efp2_init(&twisted_2Q);
    efp2_init(&next_twisted_Q);
    efp2_jacobian_init(&next_twistedJ_Q);
    efp2_jacobian_init(&twistedJ_2Q);
    for (i = 0; i < 8; i++) {
      efp2_jacobian_init(&twistedJ_Q[i]);
      efp2_jacobian_init(&twistedJ_Q_x[i]);
      efp2_jacobian_init(&twistedJ_Q_2x[i]);
      efp2_jacobian_init(&twistedJ_Q_3x[i]);
      efp2_jacobian_init(&twistedJ_Q_neg[i]);
      efp2_jacobian_init(&twistedJ_Q_x_neg[i]);
      efp2_jacobian_init(&twistedJ_Q_2x_neg[i]);
      efp2_jacobian_init(&twistedJ_Q_3x_neg[i]);
    }

    mp_limb_t A_mpn[FXLIMB2], B_mpn[FRLIMB - FXLIMB2 + 1], s0[FXLIMB], s1[FXLIMB2 - FXLIMB + 1], s2[FXLIMB], s3[FRLIMB - FXLIMB2 + 1 - FXLIMB + 1];
    //table
    efp2_jacobian_t table[4][17];
    for (i = 0; i < 17; i++) {
      efp2_jacobian_init(&table[0][i]);
      efp2_jacobian_init(&table[1][i]);
      efp2_jacobian_init(&table[2][i]);
      efp2_jacobian_init(&table[3][i]);
    }

    //set
    efp2_affine_to_jacobian_montgomery(&twistedJ_Q[0], Q);
    efp2_ecd_jacobian_lazy_montgomery(&twistedJ_2Q, &twistedJ_Q[0]);
    for (i = 1; i < 8; i++) {
      efp2_eca_jacobian_lazy_montgomery(&twistedJ_Q[i], &twistedJ_Q[i - 1], &twistedJ_2Q);
    }

    fp2_t point_table[8], inv_table[8];
    for (i = 0; i < 8; i++) fp2_set(&point_table[i], &twistedJ_Q[i].z);
    fp2_montgomery_trick_montgomery(inv_table, point_table, 8);
    for (i = 0; i < 8; i++) efp2_mix_montgomery(&twistedJ_Q[i], &twistedJ_Q[i], &inv_table[i]);

    for (i = 0; i < 8; i++) {
#ifdef X_PLUS
      efp2_jacobian_skew_frobenius_map_p1_montgomery(&twistedJ_Q_x[i], &twistedJ_Q[i]);   //twisted_Q_x
      efp2_jacobian_skew_frobenius_map_p2_montgomery(&twistedJ_Q_2x[i], &twistedJ_Q[i]);  //twisted_Q_2x
      efp2_jacobian_skew_frobenius_map_p3_montgomery(&twistedJ_Q_3x[i], &twistedJ_Q[i]);  //twisted_Q_3x

      efp2_jacobian_set_neg(&twistedJ_Q_neg[i], &twistedJ_Q[i]);        //twisted_P_neg
      efp2_jacobian_set_neg(&twistedJ_Q_x_neg[i], &twistedJ_Q_x[i]);    //twisted_P_4x_neg
      efp2_jacobian_set_neg(&twistedJ_Q_2x_neg[i], &twistedJ_Q_2x[i]);  //twisted_P_4x_neg
      efp2_jacobian_set_neg(&twistedJ_Q_3x_neg[i], &twistedJ_Q_3x[i]);  //twisted_P_4x_neg
#endif
#ifdef X_MINUS
      efp2_jacobian_skew_frobenius_map_p1_montgomery(&twistedJ_Q_x_neg[i], &twistedJ_Q[i]);   //twisted_Q_x
      efp2_jacobian_skew_frobenius_map_p2_montgomery(&twistedJ_Q_2x[i], &twistedJ_Q[i]);      //twisted_Q_2x
      efp2_jacobian_skew_frobenius_map_p3_montgomery(&twistedJ_Q_3x_neg[i], &twistedJ_Q[i]);  //twisted_Q_3x

      efp2_jacobian_set_neg(&twistedJ_Q_neg[i], &twistedJ_Q[i]);        //twisted_P_neg
      efp2_jacobian_set_neg(&twistedJ_Q_x[i], &twistedJ_Q_x_neg[i]);    //twisted_P_4x_neg
      efp2_jacobian_set_neg(&twistedJ_Q_2x_neg[i], &twistedJ_Q_2x[i]);  //twisted_P_4x_neg
      efp2_jacobian_set_neg(&twistedJ_Q_3x[i], &twistedJ_Q_3x_neg[i]);  //twisted_P_4x_neg
#endif
    }

    //set table
    table[0][0].infinity = 1;  //0
    table[1][0].infinity = 1;  //0
    table[2][0].infinity = 1;  //0
    table[3][0].infinity = 1;  //0

    for (i = 0; i < 8; i++) {
      efp2_jacobian_set(&table[0][i + 1], &twistedJ_Q[i]);
      efp2_jacobian_set(&table[0][i + 9], &twistedJ_Q_neg[i]);
      efp2_jacobian_set(&table[1][i + 1], &twistedJ_Q_x[i]);      //
      efp2_jacobian_set(&table[1][i + 9], &twistedJ_Q_x_neg[i]);  //
      efp2_jacobian_set(&table[2][i + 1], &twistedJ_Q_2x[i]);
      efp2_jacobian_set(&table[2][i + 9], &twistedJ_Q_2x_neg[i]);
      efp2_jacobian_set(&table[3][i + 1], &twistedJ_Q_3x[i]);
      efp2_jacobian_set(&table[3][i + 9], &twistedJ_Q_3x_neg[i]);
    }

    //set
    //s0,s1,s2,s3
    mpn_tdiv_qr(B_mpn, A_mpn, 0, sca->x0, FRLIMB, X2, FXLIMB2);
    mpn_tdiv_qr(s1, s0, 0, A_mpn, FXLIMB2, X_abs, FXLIMB);
    mpn_tdiv_qr(s3, s2, 0, B_mpn, FRLIMB - FXLIMB2 + 1, X_abs, FXLIMB);

    //s0[FXLIMB],s1[FXLIMB2-FXLIMB+1],s2[FXLIMB],s3[FRLIMB-FXLIMB2+1-FXLIMB+1]
    // mpn_println("s0=",s0,FXLIMB);
    // mpn_println("s1=",s1,FXLIMB2-FXLIMB+1);
    // mpn_println("s2=",s2,FXLIMB);
    // mpn_println("s3=",s3,FRLIMB-FXLIMB2+1-FXLIMB+1);

    // mpn_println("s=",sca->x0,FRLIMB);
    // getchar();

    //binary
    loop_length = 0;
    length_s[0] = (int)mpn_sizeinbase(s0, 2, FXLIMB);
    if (loop_length < length_s[0]) {
      loop_length = length_s[0];
    }
    length_s[1] = (int)mpn_sizeinbase(s1, 2, FXLIMB2 - FXLIMB + 1);
    if (loop_length < length_s[1]) {
      loop_length = length_s[1];
    }
    length_s[2] = (int)mpn_sizeinbase(s2, 2, FXLIMB);
    if (loop_length < length_s[2]) {
      loop_length = length_s[2];
    }
    length_s[3] = (int)mpn_sizeinbase(s3, 2, FRLIMB - FXLIMB2 + 1 - FXLIMB + 1);
    if (loop_length < length_s[3]) {
      loop_length = length_s[3];
    }
    // printf("loop length %d\n",loop_length);
    // getchar();

    //naf
    int naf_length[5];
    int naf_binary[4][loop_length + 1];
    for (i = 0; i < loop_length + 1; i++) {
      naf_binary[0][i] = 0;
      naf_binary[1][i] = 0;
      naf_binary[2][i] = 0;
      naf_binary[3][i] = 0;
    }
    int *naf_pointer[4];
    naf_pointer[0] = naf_binary[0];
    naf_pointer[1] = naf_binary[1];
    naf_pointer[2] = naf_binary[2];
    naf_pointer[3] = naf_binary[3];

    naf_length[1] = w_naf_frt_mpn_size(naf_binary[0], s0, FXLIMB, 5);
    naf_length[2] = w_naf_frt_mpn_size(naf_binary[1], s1, FXLIMB2 - FXLIMB + 1, 5);
    naf_length[3] = w_naf_frt_mpn_size(naf_binary[2], s2, FXLIMB, 5);
    naf_length[4] = w_naf_frt_mpn_size(naf_binary[3], s3, FRLIMB - FXLIMB2 + 1 - FXLIMB + 1, 5);
    naf_length[0] = naf_length[1];
    for (i = 2; i < 5; i++) {
      if (naf_length[0] < naf_length[i]) naf_length[0] = naf_length[i];
    }

    //naf_length=loop_length-1;
    int binary[4][naf_length[0] + 1];

    for (i = naf_length[0]; i >= 0; i--) {
      if (naf_binary[0][i] == 0)
        binary[0][i] = 0;
      else if (naf_binary[0][i] > 0)
        binary[0][i] = (naf_binary[0][i] + 1) >> 1;
      else
        binary[0][i] = ((17 - (naf_binary[0][i] + 16)) >> 1) + 8;

      if (naf_binary[1][i] == 0)
        binary[1][i] = 0;
      else if (naf_binary[1][i] > 0)
        binary[1][i] = (naf_binary[1][i] + 1) >> 1;
      else
        binary[1][i] = ((17 - (naf_binary[1][i] + 16)) >> 1) + 8;

      if (naf_binary[2][i] == 0)
        binary[2][i] = 0;
      else if (naf_binary[2][i] > 0)
        binary[2][i] = (naf_binary[2][i] + 1) >> 1;
      else
        binary[2][i] = ((17 - (naf_binary[2][i] + 16)) >> 1) + 8;

      if (naf_binary[3][i] == 0)
        binary[3][i] = 0;
      else if (naf_binary[3][i] > 0)
        binary[3][i] = (naf_binary[3][i] + 1) >> 1;
      else
        binary[3][i] = ((17 - (naf_binary[3][i] + 16)) >> 1) + 8;
    }

    next_twistedJ_Q.infinity = 1;
    for (i = 1; i < 5; i++) {
      if (naf_length[0] == naf_length[i]) {
        efp2_eca_jacobian_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[i - 1][binary[i - 1][naf_length[0]]]);
      }
    }

    //scm
    for (i = naf_length[0] - 1; i >= 0; i--) {
      efp2_ecd_jacobian_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q);
      if (binary[0][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[0][binary[0][i]]);
      if (binary[1][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[1][binary[1][i]]);
      if (binary[2][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[2][binary[2][i]]);
      if (binary[3][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[3][binary[3][i]]);
    }

    efp2_jacobian_to_affine_montgomery(ANS, &next_twistedJ_Q);
    ANS->infinity = next_twisted_Q.infinity;
  }
#else
  if (Q->infinity == 1) {
    g2_set(ANS, Q);
  } else {
    mpz_t scalar;
    mpz_init(scalar);
    mpz_set_fr(scalar, sca);
    // fr_printf("sca=",sca);
    // printf("\n");
    // mpz_printf("scalar=",scalar);
    // getchar();
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i, j, length_s[4], loop_length;
    efp2_t next_twisted_Q, twisted_Q, twisted_2Q;
    efp2_jacobian_t next_twistedJ_Q, twistedJ_2Q;
    efp2_jacobian_t twistedJ_Q[8], twistedJ_Q_x[8], twistedJ_Q_2x[8], twistedJ_Q_3x[8];
    efp2_jacobian_t twistedJ_Q_neg[8], twistedJ_Q_x_neg[8], twistedJ_Q_2x_neg[8], twistedJ_Q_3x_neg[8];

    efp2_init(&twisted_Q);
    efp2_init(&twisted_2Q);
    efp2_init(&next_twisted_Q);
    efp2_jacobian_init(&next_twistedJ_Q);
    efp2_jacobian_init(&twistedJ_2Q);
    for (i = 0; i < 8; i++) {
      efp2_jacobian_init(&twistedJ_Q[i]);
      efp2_jacobian_init(&twistedJ_Q_x[i]);
      efp2_jacobian_init(&twistedJ_Q_2x[i]);
      efp2_jacobian_init(&twistedJ_Q_3x[i]);
      efp2_jacobian_init(&twistedJ_Q_neg[i]);
      efp2_jacobian_init(&twistedJ_Q_x_neg[i]);
      efp2_jacobian_init(&twistedJ_Q_2x_neg[i]);
      efp2_jacobian_init(&twistedJ_Q_3x_neg[i]);
    }

    mpz_t A, B, s[4], x_2, x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for (i = 0; i < 4; i++) {
      mpz_init(s[i]);
    }
    //table
    efp2_jacobian_t table[4][17];
    for (i = 0; i < 17; i++) {
      efp2_jacobian_init(&table[0][i]);
      efp2_jacobian_init(&table[1][i]);
      efp2_jacobian_init(&table[2][i]);
      efp2_jacobian_init(&table[3][i]);
    }

    //set
    //efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    //TODO:lazy
    // efp2_ecd(&twisted_2Q,&twisted_Q);
    // efp2_to_montgomery(&twisted_Q,&twisted_Q);
    // efp2_to_montgomery(&twisted_2Q,&twisted_2Q);
    // efp2_affine_to_jacobian_montgomery(&twistedJ_Q[0],&twisted_Q);
    // efp2_affine_to_jacobian_montgomery(&twistedJ_2Q,&twisted_2Q);
    efp2_affine_to_jacobian_montgomery(&twistedJ_Q[0], Q);
    efp2_ecd_jacobian_lazy_montgomery(&twistedJ_2Q, &twistedJ_Q[0]);
    for (i = 1; i < 8; i++) {
      efp2_eca_jacobian_lazy_montgomery(&twistedJ_Q[i], &twistedJ_Q[i - 1], &twistedJ_2Q);
    }

    fp2_t point_table[8], inv_table[8];
    for (i = 0; i < 8; i++) fp2_set(&point_table[i], &twistedJ_Q[i].z);
    fp2_montgomery_trick_montgomery(inv_table, point_table, 8);
    for (i = 0; i < 8; i++) efp2_mix_montgomery(&twistedJ_Q[i], &twistedJ_Q[i], &inv_table[i]);

    for (i = 0; i < 8; i++) {
      efp2_jacobian_set_neg(&twistedJ_Q_neg[i], &twistedJ_Q[i]);                         //twisted_P_neg
#ifdef X_PLUS
      efp2_jacobian_skew_frobenius_map_p1_montgomery(&twistedJ_Q_x[i], &twistedJ_Q[i]);  //twisted_Q_x
#endif
#ifdef X_MINUS
      efp2_jacobian_skew_frobenius_map_p1_montgomery(&twistedJ_Q_x[i], &twistedJ_Q[i]);  //twisted_Q_x
      efp2_jacobian_set_neg(&twistedJ_Q_x[i], &twistedJ_Q_x[i]);
#endif
      efp2_jacobian_skew_frobenius_map_p2_montgomery(&twistedJ_Q_2x[i], &twistedJ_Q[i]);  //twisted_Q_2x
#ifdef X_PLUS
      efp2_jacobian_skew_frobenius_map_p3_montgomery(&twistedJ_Q_3x[i], &twistedJ_Q[i]);  //twisted_Q_3x
#endif
#ifdef X_MINUS
      efp2_jacobian_skew_frobenius_map_p3_montgomery(&twistedJ_Q_3x[i], &twistedJ_Q[i]);  //twisted_Q_3x
      efp2_jacobian_set_neg(&twistedJ_Q_3x[i], &twistedJ_Q_3x[i]);
#endif
      efp2_jacobian_set_neg(&twistedJ_Q_x_neg[i], &twistedJ_Q_x[i]);    //twisted_P_4x_neg
      efp2_jacobian_set_neg(&twistedJ_Q_2x_neg[i], &twistedJ_Q_2x[i]);  //twisted_P_4x_neg
      efp2_jacobian_set_neg(&twistedJ_Q_3x_neg[i], &twistedJ_Q_3x[i]);  //twisted_P_4x_neg
    }

    //set table
    table[0][0].infinity = 1;  //0
    table[1][0].infinity = 1;  //0
    table[2][0].infinity = 1;  //0
    table[3][0].infinity = 1;  //0

    for (i = 0; i < 8; i++) {
      efp2_jacobian_set(&table[0][i + 1], &twistedJ_Q[i]);
      efp2_jacobian_set(&table[0][i + 9], &twistedJ_Q_neg[i]);
      efp2_jacobian_set(&table[1][i + 1], &twistedJ_Q_x[i]);
      efp2_jacobian_set(&table[1][i + 9], &twistedJ_Q_x_neg[i]);
      efp2_jacobian_set(&table[2][i + 1], &twistedJ_Q_2x[i]);
      efp2_jacobian_set(&table[2][i + 9], &twistedJ_Q_2x_neg[i]);
      efp2_jacobian_set(&table[3][i + 1], &twistedJ_Q_3x[i]);
      efp2_jacobian_set(&table[3][i + 9], &twistedJ_Q_3x_neg[i]);
    }

    //set
    //s0,s1,s2,s3
    mpz_set(x_1, X_abs_z);
    mpz_mul(x_2, x_1, x_1);
    mpz_tdiv_qr(B, A, scalar, x_2);
    mpz_tdiv_qr(s[1], s[0], A, x_1);
    mpz_tdiv_qr(s[3], s[2], B, x_1);

    //binary
    loop_length = 0;
    for (i = 0; i < 4; i++) {
      length_s[i] = (int)mpz_sizeinbase(s[i], 2);
      if (loop_length < length_s[i]) {
        loop_length = length_s[i];
      }
    }

    //naf
    int naf_length[5];
    int naf_binary[4][loop_length + 1];
    for (i = 0; i < loop_length + 1; i++) {
      naf_binary[0][i] = 0;
      naf_binary[1][i] = 0;
      naf_binary[2][i] = 0;
      naf_binary[3][i] = 0;
    }
    int *naf_pointer[4];
    naf_pointer[0] = naf_binary[0];
    naf_pointer[1] = naf_binary[1];
    naf_pointer[2] = naf_binary[2];
    naf_pointer[3] = naf_binary[3];

    naf_length[1] = w_naf(naf_binary[0], s[0], 5);
    naf_length[2] = w_naf(naf_binary[1], s[1], 5);
    naf_length[3] = w_naf(naf_binary[2], s[2], 5);
    naf_length[4] = w_naf(naf_binary[3], s[3], 5);

    naf_length[0] = naf_length[1];
    for (i = 2; i < 5; i++) {
      if (naf_length[0] < naf_length[i]) naf_length[0] = naf_length[i];
    }
    //naf_length=loop_length-1;
    int binary[4][naf_length[0] + 1];

    for (i = naf_length[0]; i >= 0; i--) {
      if (naf_binary[0][i] == 0)
        binary[0][i] = 0;
      else if (naf_binary[0][i] > 0)
        binary[0][i] = (naf_binary[0][i] + 1) >> 1;
      else
        binary[0][i] = ((17 - (naf_binary[0][i] + 16)) >> 1) + 8;

      if (naf_binary[1][i] == 0)
        binary[1][i] = 0;
      else if (naf_binary[1][i] > 0)
        binary[1][i] = (naf_binary[1][i] + 1) >> 1;
      else
        binary[1][i] = ((17 - (naf_binary[1][i] + 16)) >> 1) + 8;

      if (naf_binary[2][i] == 0)
        binary[2][i] = 0;
      else if (naf_binary[2][i] > 0)
        binary[2][i] = (naf_binary[2][i] + 1) >> 1;
      else
        binary[2][i] = ((17 - (naf_binary[2][i] + 16)) >> 1) + 8;

      if (naf_binary[3][i] == 0)
        binary[3][i] = 0;
      else if (naf_binary[3][i] > 0)
        binary[3][i] = (naf_binary[3][i] + 1) >> 1;
      else
        binary[3][i] = ((17 - (naf_binary[3][i] + 16)) >> 1) + 8;
    }

    next_twistedJ_Q.infinity = 1;
    for (i = 1; i < 5; i++) {
      if (naf_length[0] == naf_length[i]) {
        efp2_eca_jacobian_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[i - 1][binary[i - 1][naf_length[0]]]);
      }
    }

    //scm
    for (i = naf_length[0] - 1; i >= 0; i--) {
      efp2_ecd_jacobian_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q);
      if (binary[0][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[0][binary[0][i]]);
      if (binary[1][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[1][binary[1][i]]);
      if (binary[2][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[2][binary[2][i]]);
      if (binary[3][i] != 0) efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q, &next_twistedJ_Q, &table[3][binary[3][i]]);
    }

    efp2_jacobian_to_affine_montgomery(ANS, &next_twistedJ_Q);
    // efp2_jacobian_to_affine_montgomery(&next_twisted_Q,&next_twistedJ_Q);
    // efp2_mod_montgomery(&next_twisted_Q,&next_twisted_Q);
    // efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity = next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for (i = 0; i < 4; i++) {
      mpz_clear(s[i]);
    }
    mpz_clear(scalar);
  }
#endif
}

void g2_set_random_with_basepoint(g2_t *ANS, g2_t *basepoint, gmp_randstate_t state) {
  fr_t s;
  fr_set_random(&s, state);
  g2_scm(ANS, basepoint, &s);
}

void g2_test(int scm) {
  g2_t A_g2;
  efp12_t A_efp12;
  fr_t sca_fr;
  mpz_t sca;
  //initialize
  mpz_init(sca);
  fr_init(&sca_fr);
  g2_init(&A_g2);
  efp12_init(&A_efp12);
  //set
  bls12_generate_g2(&A_efp12);
  efp12_printf("A=", &A_efp12);

  efp12_to_efp2(&A_g2, &A_efp12);
  efp2_to_montgomery(&A_g2, &A_g2);
  fr_set_random(&sca_fr, state);
  mpz_set_fr(sca, &sca_fr);
  //calc
  printf("\n===========================\n");
  g2_scm(&A_g2, &A_g2, &sca_fr);
  efp12_scm(&A_efp12, &A_efp12, sca);
  efp12_printf("A=", &A_efp12);
  printf("\n\n");
  g2_printf("A=", &A_g2);
  printf("\n");
  printf("g2 test:");
  if (g2_cmp_efp12(&A_g2, &A_efp12) == 0)
    printf("ok!\n");
  else
    printf("ng\n");

  //bench
  int i = 0;
  float scm_time = 0;
  cost tmp, scm_cost;
  struct timeval tv_A, tv_B;
  for (i = 0; i < scm; i++) {
    fr_set_random(&sca_fr, state);
    //faster type
    cost_zero();
    gettimeofday(&tv_A, NULL);
    g2_scm(&A_g2, &A_g2, &sca_fr);
    gettimeofday(&tv_B, NULL);
    scm_time += timedifference_msec(tv_A, tv_B);
    cost_check(&tmp);
    cost_addition(&scm_cost, &tmp);
  }
  printf("bls12 g2 scm.     : %.4f[ms]\n", scm_time / scm);

#ifdef DEBUG_COST_A
  printf("*********bls12 g2 scm fp COST.********         \n");
  cost_printf("bls12 g2 scm", &scm_cost, scm);
  printf("***************************************         \n");
#endif
  mpz_clear(sca);
}

void g2_set_random_test(int scm) {
  g2_t P_test;
  int i;
  float random_time = 0, random_fast_time = 0;
  struct timeval tv_A, tv_B;
  for (i = 0; i < scm; i++) {
    //normal type
    // gettimeofday(&tv_A, NULL);
    // g2_set_random_schoolbook(&P_test, state);
    // gettimeofday(&tv_B, NULL);
    // random_time += timedifference_msec(tv_A, tv_B);

    //faster type
    gettimeofday(&tv_A, NULL);
    g2_set_random(&P_test, state);
    gettimeofday(&tv_B, NULL);
    random_fast_time += timedifference_msec(tv_A, tv_B);
  }
  //printf("g2_set_random_schoolbook.     : %.4f[ms]\n", random_time / scm);
  printf("g2_set_random.          : %.4f[ms]\n", random_fast_time / scm);

}
/*******************g3_t******************************/
void g3_init(g3_t *A) {
  fp12_init(A);
}
void g3_set(g3_t *A, g3_t *B) {
  fp12_set(A, B);
}
int g3_cmp(g3_t *A, g3_t *B) {
  return fp12_cmp(A, B);
}
void g3_mul(g3_t *ANS, g3_t *A, g3_t *B) {
  fp12_mul_lazy_montgomery(ANS, A, B);
}
void g3_sqr(g3_t *ANS, g3_t *A) {
  //fp12_sqr_lazy_montgomery(ANS,A);
  fp12_sqr_GS_lazy_montgomery(ANS, A);
}
void g3_inv(g3_t *ANS, g3_t *A) {
  fp12_inv_lazy_montgomery(ANS, A);
}
void g3_printf(char *s, g3_t *A) {
  g3_t tmp;
  fp12_mod_montgomery(&tmp, A);
  fp12_printf(s, &tmp);
}

void g3_println(char *s, g3_t *A) {
  g3_t tmp;
  fp12_mod_montgomery(&tmp, A);
  fp12_println(s, &tmp);
}
#if ARCBIT == 64
void g3_exp(g3_t *ANS, g3_t *A, fr_t *sca) {
  int i, length_s[4], loop_length;
  fp12_t Buf;
  fp12_init(&Buf);
  fp12_t next_f, f2, f[8], f_2x[8], f_4x[8], f_6x[8];
  fp12_t f_neg[8], f_2x_neg[8], f_4x_neg[8], f_6x_neg[8];

  fp12_init(&next_f);
  fp12_init(&f2);
  for (i = 0; i < 8; i++) {
    fp12_init(&f[i]);
    fp12_init(&f_2x[i]);
    fp12_init(&f_4x[i]);
    fp12_init(&f_6x[i]);
    fp12_init(&f_neg[i]);
    fp12_init(&f_2x_neg[i]);
    fp12_init(&f_4x_neg[i]);
    fp12_init(&f_6x_neg[i]);
  }
  mp_limb_t A_mpn[FXLIMB2], B_mpn[FRLIMB - FXLIMB2 + 1], s0[FXLIMB], s1[FXLIMB2 - FXLIMB + 1], s2[FXLIMB], s3[FRLIMB - FXLIMB2 + 1 - FXLIMB + 1];

  //table
  fp12_t table[4][17];
  for (i = 0; i < 17; i++) {
    fp12_init(&table[0][i]);
    fp12_init(&table[1][i]);
    fp12_init(&table[2][i]);
    fp12_init(&table[3][i]);
  }

  //set
  fp12_set(&f[0], A);

  fp12_sqr_lazy_montgomery(&f2, &f[0]);
  for (i = 1; i < 8; i++) {
    fp12_mul_lazy_montgomery(&f[i], &f[i - 1], &f2);
  }

  for (i = 0; i < 8; i++) {
#ifdef X_PLUS
    fp12_frobenius_map_p6_montgomery(&f_neg[i], &f[i]);        //twisted_P_neg
    fp12_frobenius_map_p1_lazy_montgomery(&f_2x[i], &f[i]);    //f_2x
    fp12_frobenius_map_p2_montgomery(&f_4x[i], &f[i]);         //f_4x
    fp12_frobenius_map_p3_lazy_montgomery(&f_6x[i], &f[i]);    //f_6x
    fp12_frobenius_map_p6_montgomery(&f_2x_neg[i], &f_2x[i]);  //f_2x
    fp12_frobenius_map_p6_montgomery(&f_4x_neg[i], &f_4x[i]);  //f_4x
    fp12_frobenius_map_p6_montgomery(&f_6x_neg[i], &f_6x[i]);  //f_6x
#endif
#ifdef X_MINUS
    fp12_frobenius_map_p6_montgomery(&f_neg[i], &f[i]);          //twisted_P_neg
    fp12_frobenius_map_p1_lazy_montgomery(&f_2x_neg[i], &f[i]);  //f_2x
    fp12_frobenius_map_p2_montgomery(&f_4x[i], &f[i]);           //f_4x
    fp12_frobenius_map_p3_lazy_montgomery(&f_6x_neg[i], &f[i]);  //f_6x
    fp12_frobenius_map_p6_montgomery(&f_2x[i], &f_2x_neg[i]);    //f_2x
    fp12_frobenius_map_p6_montgomery(&f_4x_neg[i], &f_4x[i]);    //f_4x
    fp12_frobenius_map_p6_montgomery(&f_6x[i], &f_6x_neg[i]);    //f_6x
#endif
  }

  //set table
  fp12_set_mpn(&table[0][0], RmodP);
  fp12_set_mpn(&table[1][0], RmodP);
  fp12_set_mpn(&table[2][0], RmodP);
  fp12_set_mpn(&table[3][0], RmodP);

  for (i = 0; i < 8; i++) {
    fp12_set(&table[0][i + 1], &f[i]);
    fp12_set(&table[0][i + 9], &f_neg[i]);
    fp12_set(&table[1][i + 1], &f_2x[i]);
    fp12_set(&table[1][i + 9], &f_2x_neg[i]);
    fp12_set(&table[2][i + 1], &f_4x[i]);
    fp12_set(&table[2][i + 9], &f_4x_neg[i]);
    fp12_set(&table[3][i + 1], &f_6x[i]);
    fp12_set(&table[3][i + 9], &f_6x_neg[i]);
  }

  //set
  //s0,s1,s2,s3
  mpn_tdiv_qr(B_mpn, A_mpn, 0, sca->x0, FRLIMB, X2, FXLIMB2);
  mpn_tdiv_qr(s1, s0, 0, A_mpn, FXLIMB2, X_abs, FXLIMB);
  mpn_tdiv_qr(s3, s2, 0, B_mpn, FRLIMB - FXLIMB2 + 1, X_abs, FXLIMB);

  length_s[0] = (int)mpn_sizeinbase(s0, 2, FXLIMB);
  if (loop_length < length_s[0]) {
    loop_length = length_s[0];
  }
  length_s[1] = (int)mpn_sizeinbase(s1, 2, FXLIMB2 - FXLIMB + 1);
  if (loop_length < length_s[1]) {
    loop_length = length_s[1];
  }
  length_s[2] = (int)mpn_sizeinbase(s2, 2, FXLIMB);
  if (loop_length < length_s[2]) {
    loop_length = length_s[2];
  }
  length_s[3] = (int)mpn_sizeinbase(s3, 2, FRLIMB - FXLIMB2 + 1 - FXLIMB + 1);
  if (loop_length < length_s[3]) {
    loop_length = length_s[3];
  }

  //naf
  int naf_length[5];
  int naf_binary[4][loop_length + 1];
  for (i = 0; i < loop_length + 1; i++) {
    naf_binary[0][i] = 0;
    naf_binary[1][i] = 0;
    naf_binary[2][i] = 0;
    naf_binary[3][i] = 0;
  }
  int *naf_pointer[4];
  naf_pointer[0] = naf_binary[0];
  naf_pointer[1] = naf_binary[1];
  naf_pointer[2] = naf_binary[2];
  naf_pointer[3] = naf_binary[3];

  naf_length[1] = w_naf_frt_mpn_size(naf_binary[0], s0, FXLIMB, 5);
  naf_length[2] = w_naf_frt_mpn_size(naf_binary[1], s1, FXLIMB2 - FXLIMB + 1, 5);
  naf_length[3] = w_naf_frt_mpn_size(naf_binary[2], s2, FXLIMB, 5);
  naf_length[4] = w_naf_frt_mpn_size(naf_binary[3], s3, FRLIMB - FXLIMB2 + 1 - FXLIMB + 1, 5);

  naf_length[0] = naf_length[1];
  for (i = 2; i < 5; i++) {
    if (naf_length[0] < naf_length[i]) naf_length[0] = naf_length[i];
  }
  //naf_length=loop_length-1;
  int binary[4][naf_length[0] + 1];

  for (i = naf_length[0]; i >= 0; i--) {
    if (naf_binary[0][i] == 0)
      binary[0][i] = 0;
    else if (naf_binary[0][i] > 0)
      binary[0][i] = (naf_binary[0][i] + 1) >> 1;
    else
      binary[0][i] = ((17 - (naf_binary[0][i] + 16)) >> 1) + 8;

    if (naf_binary[1][i] == 0)
      binary[1][i] = 0;
    else if (naf_binary[1][i] > 0)
      binary[1][i] = (naf_binary[1][i] + 1) >> 1;
    else
      binary[1][i] = ((17 - (naf_binary[1][i] + 16)) >> 1) + 8;

    if (naf_binary[2][i] == 0)
      binary[2][i] = 0;
    else if (naf_binary[2][i] > 0)
      binary[2][i] = (naf_binary[2][i] + 1) >> 1;
    else
      binary[2][i] = ((17 - (naf_binary[2][i] + 16)) >> 1) + 8;

    if (naf_binary[3][i] == 0)
      binary[3][i] = 0;
    else if (naf_binary[3][i] > 0)
      binary[3][i] = (naf_binary[3][i] + 1) >> 1;
    else
      binary[3][i] = ((17 - (naf_binary[3][i] + 16)) >> 1) + 8;
  }
  fp12_set_mpn(&next_f, RmodP);
  for (i = 1; i < 5; i++) {
    if (naf_length[0] == naf_length[i]) {
      fp12_mul_lazy_montgomery(&next_f, &next_f, &table[i - 1][binary[i - 1][naf_length[0]]]);
    }
  }

  //scm
  for (i = naf_length[0] - 1; i >= 0; i--) {
    fp12_sqr_GS_lazy_montgomery(&next_f, &next_f);
    if (binary[0][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[0][binary[0][i]]);
    if (binary[1][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[1][binary[1][i]]);
    if (binary[2][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[2][binary[2][i]]);
    if (binary[3][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[3][binary[3][i]]);
  }
  fp12_set(ANS, &next_f);
}
#else
void g3_exp(fp12_t *ANS, fp12_t *A, fr_t *sca) {
  mpz_t scalar;
  mpz_init(scalar);
  mpz_set_fr(scalar, sca);
  //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
  int i, length_s[4], loop_length;
  fp12_t Buf;
  fp12_init(&Buf);
  fp12_t next_f, f2, f[8], f_2x[8], f_4x[8], f_6x[8];
  fp12_t f_neg[8], f_2x_neg[8], f_4x_neg[8], f_6x_neg[8];

  fp12_init(&next_f);
  fp12_init(&f2);
  for (i = 0; i < 8; i++) {
    fp12_init(&f[i]);
    fp12_init(&f_2x[i]);
    fp12_init(&f_4x[i]);
    fp12_init(&f_6x[i]);
    fp12_init(&f_neg[i]);
    fp12_init(&f_2x_neg[i]);
    fp12_init(&f_4x_neg[i]);
    fp12_init(&f_6x_neg[i]);
  }

  mpz_t A_s, B_s, s[4], x_4, x_2;
  mpz_init(A_s);
  mpz_init(B_s);
  mpz_init(x_2);
  mpz_init(x_4);
  for (i = 0; i < 4; i++) {
    mpz_init(s[i]);
  }

  //table
  fp12_t table[4][17];
  for (i = 0; i < 17; i++) {
    fp12_init(&table[0][i]);
    fp12_init(&table[1][i]);
    fp12_init(&table[2][i]);
    fp12_init(&table[3][i]);
  }

  //set
  //fp12_to_montgomery(&f[0],A);
  fp12_set(&f[0], A);

  fp12_sqr_lazy_montgomery(&f2, &f[0]);
  for (i = 1; i < 8; i++) {
    fp12_mul_lazy_montgomery(&f[i], &f[i - 1], &f2);
  }

  for (i = 0; i < 8; i++) {
#ifdef X_PLUS
    fp12_frobenius_map_p6_montgomery(&f_neg[i], &f[i]);          //twisted_P_neg
    fp12_frobenius_map_p1_lazy_montgomery(&f_2x[i], &f[i]);      //f_2x
    fp12_frobenius_map_p2_montgomery(&f_4x[i], &f[i]);           //f_4x
    fp12_frobenius_map_p3_lazy_montgomery(&f_6x[i], &f[i]);      //f_6x
    fp12_frobenius_map_p6_montgomery(&f_2x_neg[i], &f_2x[i]);    //f_2x
    fp12_frobenius_map_p6_montgomery(&f_4x_neg[i], &f_4x[i]);    //f_4x
    fp12_frobenius_map_p6_montgomery(&f_6x_neg[i], &f_6x[i]);    //f_6x
#endif
#ifdef X_MINUS
    fp12_frobenius_map_p6_montgomery(&f_neg[i], &f[i]);          //twisted_P_neg
    fp12_frobenius_map_p1_lazy_montgomery(&f_2x_neg[i], &f[i]);  //f_2x
    fp12_frobenius_map_p2_montgomery(&f_4x[i], &f[i]);           //f_4x
    fp12_frobenius_map_p3_lazy_montgomery(&f_6x_neg[i], &f[i]);  //f_6x
    fp12_frobenius_map_p6_montgomery(&f_2x[i], &f_2x_neg[i]);    //f_2x
    fp12_frobenius_map_p6_montgomery(&f_4x_neg[i], &f_4x[i]);    //f_4x
    fp12_frobenius_map_p6_montgomery(&f_6x[i], &f_6x_neg[i]);    //f_6x
#endif
  }

  //set table
  fp12_set_mpn(&table[0][0], RmodP);
  fp12_set_mpn(&table[1][0], RmodP);
  fp12_set_mpn(&table[2][0], RmodP);
  fp12_set_mpn(&table[3][0], RmodP);

  for (i = 0; i < 8; i++) {
    fp12_set(&table[0][i + 1], &f[i]);
    fp12_set(&table[0][i + 9], &f_neg[i]);
    fp12_set(&table[1][i + 1], &f_2x[i]);
    fp12_set(&table[1][i + 9], &f_2x_neg[i]);
    fp12_set(&table[2][i + 1], &f_4x[i]);
    fp12_set(&table[2][i + 9], &f_4x_neg[i]);
    fp12_set(&table[3][i + 1], &f_6x[i]);
    fp12_set(&table[3][i + 9], &f_6x_neg[i]);
  }

  //set
  //s0,s1,s2,s3
  mpz_set(x_2, X_abs_z);
  mpz_mul(x_4, x_2, x_2);
  mpz_tdiv_qr(B_s, A_s, scalar, x_4);
  mpz_tdiv_qr(s[1], s[0], A_s, x_2);
  mpz_tdiv_qr(s[3], s[2], B_s, x_2);

  //binary
  loop_length = 0;
  for (i = 0; i < 4; i++) {
    length_s[i] = (int)mpz_sizeinbase(s[i], 2);
    if (loop_length < length_s[i]) {
      loop_length = length_s[i];
    }
  }

  //naf
  int naf_length[5];
  int naf_binary[4][loop_length + 1];
  for (i = 0; i < loop_length + 1; i++) {
    naf_binary[0][i] = 0;
    naf_binary[1][i] = 0;
    naf_binary[2][i] = 0;
    naf_binary[3][i] = 0;
  }
  int *naf_pointer[4];
  naf_pointer[0] = naf_binary[0];
  naf_pointer[1] = naf_binary[1];
  naf_pointer[2] = naf_binary[2];
  naf_pointer[3] = naf_binary[3];

  naf_length[1] = w_naf(naf_binary[0], s[0], 5);
  naf_length[2] = w_naf(naf_binary[1], s[1], 5);
  naf_length[3] = w_naf(naf_binary[2], s[2], 5);
  naf_length[4] = w_naf(naf_binary[3], s[3], 5);

  naf_length[0] = naf_length[1];
  for (i = 2; i < 5; i++) {
    if (naf_length[0] < naf_length[i]) naf_length[0] = naf_length[i];
  }
  //naf_length=loop_length-1;
  int binary[4][naf_length[0] + 1];

  for (i = naf_length[0]; i >= 0; i--) {
    if (naf_binary[0][i] == 0)
      binary[0][i] = 0;
    else if (naf_binary[0][i] > 0)
      binary[0][i] = (naf_binary[0][i] + 1) >> 1;
    else
      binary[0][i] = ((17 - (naf_binary[0][i] + 16)) >> 1) + 8;

    if (naf_binary[1][i] == 0)
      binary[1][i] = 0;
    else if (naf_binary[1][i] > 0)
      binary[1][i] = (naf_binary[1][i] + 1) >> 1;
    else
      binary[1][i] = ((17 - (naf_binary[1][i] + 16)) >> 1) + 8;

    if (naf_binary[2][i] == 0)
      binary[2][i] = 0;
    else if (naf_binary[2][i] > 0)
      binary[2][i] = (naf_binary[2][i] + 1) >> 1;
    else
      binary[2][i] = ((17 - (naf_binary[2][i] + 16)) >> 1) + 8;

    if (naf_binary[3][i] == 0)
      binary[3][i] = 0;
    else if (naf_binary[3][i] > 0)
      binary[3][i] = (naf_binary[3][i] + 1) >> 1;
    else
      binary[3][i] = ((17 - (naf_binary[3][i] + 16)) >> 1) + 8;
  }
  fp12_set_mpn(&next_f, RmodP);
  for (i = 1; i < 5; i++) {
    if (naf_length[0] == naf_length[i]) {
      fp12_mul_lazy_montgomery(&next_f, &next_f, &table[i - 1][binary[i - 1][naf_length[0]]]);
    }
  }

  //scm
  for (i = naf_length[0] - 1; i >= 0; i--) {
    fp12_sqr_GS_lazy_montgomery(&next_f, &next_f);
    if (binary[0][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[0][binary[0][i]]);
    if (binary[1][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[1][binary[1][i]]);
    if (binary[2][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[2][binary[2][i]]);
    if (binary[3][i] != 0) fp12_mul_lazy_montgomery(&next_f, &next_f, &table[3][binary[3][i]]);
  }
  //fp12_set(ANS,&next_f);
  //fp12_mod_montgomery(ANS,&next_f);
  fp12_set(ANS, &next_f);

  mpz_clear(x_2);
  mpz_clear(x_4);
  mpz_clear(scalar);
  for (i = 0; i < 4; i++) {
    mpz_clear(s[i]);
  }
}
#endif
/*******************pairing***************************/
void g1g2_to_g3_miller_algo(g3_t *ANS, g1_t *P, g2_t *Q) {
  efp2_projective_t T;
  efp2_projective_t mapped_Q, mapped_Q_neg;
  efp_t mapped_P_ltt, mapped_P_ltq;
  fp12_t f;
  int i;

  efp2_projective_init(&T);
  efp2_projective_init(&mapped_Q);
  efp2_projective_init(&mapped_Q_neg);
  efp_init(&mapped_P_ltt);
  efp_init(&mapped_P_ltq);
  fp12_init(&f);

  //set
  efp_set(&mapped_P_ltq, P);  //set mapped_P
  efp2_affine_to_projective_montgomery(&mapped_Q, Q);
  efp2_projective_set(&mapped_Q_neg, &mapped_Q);  //set mapped_Q_neg
  fp2_set_neg(&mapped_Q_neg.y, &mapped_Q_neg.y);
  efp2_projective_set(&T, &mapped_Q);  //set T
  fp_set_mpn(&f.x0.x0.x0, RmodP);

  //precompute P
  fp_add_nonmod_single(&mapped_P_ltt.x, &mapped_P_ltq.x, &mapped_P_ltq.x);
  fp_add_nonmod_single(&mapped_P_ltt.x, &mapped_P_ltt.x, &mapped_P_ltq.x);
  fp_set_neg(&mapped_P_ltt.y, &mapped_P_ltq.y);
  fp_set_neg(&mapped_P_ltq.x, &mapped_P_ltq.x);

#ifdef X_PLUS
  //miller
  for (i = bls12_X_length - 1; i >= 0; i--) {
    switch (bls12_X_binary[i]) {
      case 0:
        ff_ltt_projective_lazy_montgomery(&f, &T, &mapped_P_ltt);
        break;
      case 1:
        ff_ltt_projective_lazy_montgomery(&f, &T, &mapped_P_ltt);
        f_ltq_projective_lazy_montgomery(&f, &T, &mapped_Q, &mapped_P_ltq);
        break;
      case -1:
        ff_ltt_projective_lazy_montgomery(&f, &T, &mapped_P_ltt);
        f_ltq_projective_lazy_montgomery(&f, &T, &mapped_Q_neg, &mapped_P_ltq);
        break;
      default:
        break;
    }
  }
#endif
#ifdef X_MINUS
  //miller
  for (i = bls12_X_length - 1; i >= 0; i--) {
    switch (bls12_X_binary[i]) {
      case 0:
        ff_ltt_projective_lazy_montgomery(&f, &T, &mapped_P_ltt);
        break;
      case -1:
        ff_ltt_projective_lazy_montgomery(&f, &T, &mapped_P_ltt);
        f_ltq_projective_lazy_montgomery(&f, &T, &mapped_Q, &mapped_P_ltq);
        break;
      case 1:
        ff_ltt_projective_lazy_montgomery(&f, &T, &mapped_P_ltt);
        f_ltq_projective_lazy_montgomery(&f, &T, &mapped_Q_neg, &mapped_P_ltq);
        break;
      default:
        break;
    }
  }
#endif
  fp12_set(ANS, &f);
}
void g1g2_to_g3_miller_algo_affine(g3_t *ANS, g1_t *P, g2_t *Q) {
  efp2_t T;
  efp2_t mapped_Q, mapped_Q_neg;
  efp_t mapped_P;
  fp12_t f;
  fp_t L;
  int i;

  efp2_init(&T);
  efp2_init(&mapped_Q);
  efp2_init(&mapped_Q_neg);
  efp_init(&mapped_P);
  fp12_init(&f);
  fp_init(&L);

  //set
  efp_set(&mapped_P, P);   //set mapped_P
  efp2_set(&mapped_Q, Q);  //set mapped_Q
  Pseudo_8_sparse_mapping_montgomery(&mapped_P, &mapped_Q, &L);
  efp2_set(&mapped_Q_neg, &mapped_Q);  //set mapped_Q_neg
  fp2_set_neg(&mapped_Q_neg.y, &mapped_Q_neg.y);
  efp2_set(&T, &mapped_Q);  //set T
  //TODO:1->RmodP?
  //fp_set_ui(&f.x0.x0.x0,1);
  fp_set_mpn(&f.x0.x0.x0, RmodP);

//miller
#ifdef X_PLUS
  for (i = bls12_X_length - 1; i >= 0; i--) {
    switch (bls12_X_binary[i]) {
      case 0:
        ff_ltt_lazy_montgomery(&f, &T, &mapped_P, &L);
        break;
      case 1:
        ff_ltt_lazy_montgomery(&f, &T, &mapped_P, &L);
        f_ltq_lazy_montgomery(&f, &T, &mapped_Q, &mapped_P, &L);
        break;
      case -1:
        ff_ltt_lazy_montgomery(&f, &T, &mapped_P, &L);
        f_ltq_lazy_montgomery(&f, &T, &mapped_Q_neg, &mapped_P, &L);
        break;
      default:
        break;
    }
  }
#endif

#ifdef X_MINUS
  //miller
  for (i = bls12_X_length - 1; i >= 0; i--) {
    switch (bls12_X_binary[i]) {
      case 0:
        ff_ltt_lazy_montgomery(&f, &T, &mapped_P, &L);
        break;
      case -1:
        ff_ltt_lazy_montgomery(&f, &T, &mapped_P, &L);
        f_ltq_lazy_montgomery(&f, &T, &mapped_Q, &mapped_P, &L);
        break;
      case 1:
        ff_ltt_lazy_montgomery(&f, &T, &mapped_P, &L);
        f_ltq_lazy_montgomery(&f, &T, &mapped_Q_neg, &mapped_P, &L);
        break;
      default:
        break;
    }
  }
#endif

  fp12_set(ANS, &f);
}
void g1g2_to_g3_hardpart_compress_montrick(fp12_t *ANS, fp12_t *A) {
  // fp12_t tmp,t0,t1,t2,t3,t4,t5, test,At;
  // fp12_init(&tmp);
  // fp12_init(&t0);
  // fp12_init(&t1);
  // fp12_init(&t2);
  // fp12_init(&t3);
  // fp12_init(&t5);
  // fp12_init(&t4);
  // fp12_init(&test);
  // fp12_init(&At);
  // //HARDPART
  // fp12_sqr_GS_lazy_montgomery(&t0, A);
  // bls12_fp12_pow_X_compress_montrick(&t1, &t0);

  // bls12_fp12_pow_X2_compress_montrick(&t2,&t1);//t2:=t1^(u2);
  // fp12_frobenius_map_p6_montgomery(&t3,A);//t3:=f^(-1);
  // fp12_mul_lazy_montgomery(&t1,&t3,&t1);//t1:=t3*t1;
  // fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);
  // fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;

  // bls12_fp12_pow_X_compress_montrick(&t2,&t1);//t2:=t1^(u);
  // bls12_fp12_pow_X_compress_montrick(&t3,&t2);//t3:=t2^(u);
  // fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);

  // fp12_mul_lazy_montgomery(&t3,&t1,&t3);//t3:=t1*t3;
  // fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);
  // fp12_frobenius_map_p3_lazy_montgomery(&t1,&t1);//t1:=t1^(p^3);
  // fp12_frobenius_map_p2_montgomery(&t2,&t2);//t2:=t2^(p^2);

  // fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;
  // bls12_fp12_pow_X_compress_montrick(&t2,&t3);//t2:=t3^(u);
  // fp12_mul_lazy_montgomery(&t2,&t2,&t0);//t2:=t2*t0;
  // fp12_mul_lazy_montgomery(&t2,&t2,A);//t2:=t2*f;
  // fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;

  // fp12_frobenius_map_p1_lazy_montgomery(&t2,&t3);//t2:=t3^p;
  // fp12_mul_lazy_montgomery(ANS,&t1,&t2);//t1:=t1*t2;
  //fp12_mod_montgomery(ANS,ANS);
  ////hayashida et al. in 2021
  fp12_t t1, t2, t3, t4;
  fp12_init(&t1);
  fp12_init(&t2);
  fp12_init(&t3);
  fp12_init(&t4);
  //HARDPART
  fp12_sqr_GS_lazy_montgomery(&t1, A);
  bls12_fp12_pow_X_compress_montrick(&t2, &t1);
  bls12_fp12_pow_X2_compress_montrick(&t3, &t2);  //t3:=t2^(u2);
  fp12_frobenius_map_p6_montgomery(&t2, &t2);     //t4:=f^(-1);
  fp12_mul_lazy_montgomery(&t2, &t2, &t3);        //t2:=t4*t2;
  fp12_mul_lazy_montgomery(&t2, &t2, A);          //t2:=t4*t2;
  bls12_fp12_pow_X_compress_montrick(&t3, &t2);
  fp12_frobenius_map_p1_lazy_montgomery(&t2, &t2);
  fp12_mul_lazy_montgomery(&t2, &t2, &t3);  //t2:=t4*t2;
  bls12_fp12_pow_X_compress_montrick(&t3, &t2);
  bls12_fp12_pow_X_compress_montrick(&t3, &t3);
  fp12_frobenius_map_p2_montgomery(&t4, &t2);
  fp12_frobenius_map_p6_montgomery(&t2, &t2);  //t4:=f^(-1);
  fp12_mul_lazy_montgomery(&t2, &t2, &t3);     //t2:=t4*t2;
  fp12_mul_lazy_montgomery(&t2, &t2, &t4);     //t2:=t4*t2;
  fp12_mul_lazy_montgomery(&t1, &t1, A);       //t2:=t4*t2;
  fp12_mul_lazy_montgomery(ANS, &t2, &t1);     //t2:=t4*t2;
}
void g1g2_to_g3_final_exp(g3_t *ANS, g3_t *A) {
  fp12_t tmp, t0, t1, t2, t3, t4, t5, test, At;
  fp12_init(&tmp);
  fp12_init(&t0);
  fp12_init(&t1);
  fp12_init(&t2);
  fp12_init(&t3);
  fp12_init(&t5);
  fp12_init(&t4);
  fp12_init(&test);
  fp12_init(&At);

  fp12_set(&At, A);

  //EASY PART
  //f←f^(p^6)*f^-1
  fp12_frobenius_map_p6_montgomery(&t0, &At);  //f^(p^6)
  fp12_inv_lazy_montgomery(&t1, &At);          //f^-1
  fp12_mul_lazy_montgomery(&tmp, &t0, &t1);    //f^(p^6)*f^-1

  //f←f^(p^2)*f
  fp12_frobenius_map_p2_montgomery(&t0, &tmp);  //f^(p^2)
  fp12_mul_lazy_montgomery(&tmp, &t0, &tmp);    //f^(p^2)*f

  g1g2_to_g3_hardpart_compress_montrick(ANS, &tmp);
}
void g1g2_to_g3_pairing(g3_t *ANS, g1_t *P, g2_t *Q) {
#ifdef DEBUG_COST_A
  cost tmp;
#endif

  //Miller's Algo.
  gettimeofday(&tv_start, NULL);
  g1g2_to_g3_miller_algo(ANS, P, Q);
  gettimeofday(&tv_end, NULL);
  MILLER_OPT_PROJECTIVE += timedifference_msec(tv_start, tv_end);

#ifdef DEBUG_COST_A
  cost_check(&tmp);
  cost_addition(&MILLER_OPT_PROJECTIVE_COST, &tmp);
#endif

  //Final Exp.
  gettimeofday(&tv_start, NULL);
  g1g2_to_g3_final_exp(ANS, ANS);
  gettimeofday(&tv_end, NULL);
  FINALEXP_OPT_PROJECTIVE += timedifference_msec(tv_start, tv_end);
}
void g1g2_to_g3_pairing_affine(g3_t *ANS, g1_t *P, g2_t *Q) {
#ifdef DEBUG_COST_A
  cost tmp;
#endif

  //Miller's Algo.
  gettimeofday(&tv_start, NULL);
  g1g2_to_g3_miller_algo_affine(ANS, P, Q);
  gettimeofday(&tv_end, NULL);
  MILLER_OPT_AFFINE += timedifference_msec(tv_start, tv_end);

#ifdef DEBUG_COST_A
  cost_check(&tmp);
  cost_addition(&MILLER_OPT_AFFINE_COST, &tmp);
#endif

  //Final Exp.
  gettimeofday(&tv_start, NULL);
  g1g2_to_g3_final_exp(ANS, ANS);
  gettimeofday(&tv_end, NULL);
  FINALEXP_OPT_AFFINE += timedifference_msec(tv_start, tv_end);
}
void billinear_test() {
  fr_t a, b, c;
  g1_t P, aP;
  g2_t Q, bQ;
  g3_t test1, test2;

  fr_set_random(&a, state);
  fr_set_random(&b, state);
  fr_mul(&c, &a, &b);

  g1_set_random(&P, state);
  g2_set_random(&Q, state);
  g1_scm(&aP, &P, &a);
  g2_scm(&bQ, &Q, &b);

  //test1={Pairng(P,Q)^c
  g1g2_to_g3_pairing(&test1, &P, &Q);
  g3_exp(&test1, &test1, &c);

  //test2=Pairng(aP,bQ)
  g1g2_to_g3_pairing(&test2, &aP, &bQ);

  if (fp12_cmp(&test1, &test2) == 0)
    printf("ok!\n");
  else
    printf("ng\n");
}
int debug_pairing(int pairing) {
  int i, n = 0;
  float opt_time = 0, opt_affine_time = 0;
  cost tmp, opt_cost, opt_affine_cost;
  struct timeval tv_A, tv_B;
  printf("====================================================================================\n");
  printf("bls12_Opt-ate pairing\n\n");
  g1_t P, s1P, s2P;
  g2_t Q, s1Q, s2Q;
  g1_init(&P);
  g1_init(&s1P);
  g1_init(&s2P);
  g2_init(&Q);
  g2_init(&s1Q);
  g2_init(&s2Q);

  g3_t Z, testA, testB, testC, test1, test2, test3;
  g3_init(&Z);
  g3_init(&testA);
  g3_init(&testB);
  g3_init(&testC);
  g3_init(&test1);
  g3_init(&test2);
  g3_init(&test3);

  cost_init(&tmp);
  cost_init(&opt_cost);
  cost_init(&opt_affine_cost);

  fr_t s12, s1, s2;
  fr_init(&s12);
  fr_init(&s1);
  fr_init(&s2);

  fr_set_random(&s1, state);
  fr_set_random(&s2, state);
  fr_mul(&s12, &s1, &s2);

  g1_set_random_schoolbook(&P, state);
  g2_set_random(&Q, state);

  g1_scm(&s1P, &P, &s1);
  g1_scm(&s2P, &P, &s2);
  g2_scm(&s1Q, &Q, &s1);
  g2_scm(&s2Q, &Q, &s2);

  g1g2_to_g3_pairing(&Z, &P, &Q);
  g3_exp(&testA, &Z, &s12);
  g1g2_to_g3_pairing(&testB, &s1P, &s2Q);
  g1g2_to_g3_pairing(&testC, &s2P, &s1Q);
  fp12_printf("testA=", &testA);

  printf("bilinear test\n");
  if (fp12_cmp(&testA, &testB) != 0 || fp12_cmp(&testA, &testC) != 0) {
    printf("bilinear failed!!\n\n");
    return 1;
  }
  printf("bilinear succeced!!\n\n");

  MILLER_OPT_PROJECTIVE = 0;
  FINALEXP_OPT_PROJECTIVE = 0;
  MILLER_OPT_AFFINE = 0;
  FINALEXP_OPT_AFFINE = 0;
  opt_time = 0;
  opt_affine_time = 0;
  cost_init(&MILLER_OPT_PROJECTIVE_COST);
  cost_init(&FINALEXP_OPT_PROJECTIVE_COST);
  cost_init(&MILLER_OPT_AFFINE_COST);
  cost_init(&FINALEXP_OPT_AFFINE_COST);

  g2_set_random(&Q, state);
  for (i = 0; i < pairing; i++) {
    g1_set_random_with_basepoint(&P, &P, state);
    g2_set_random_with_basepoint(&Q, &Q, state);
    //g1g2_to_g3_pairing(&test1,&P,&Q);
    cost_zero();
    gettimeofday(&tv_A, NULL);
    g1g2_to_g3_pairing(&test2, &P, &Q);
    gettimeofday(&tv_B, NULL);
    opt_time += timedifference_msec(tv_A, tv_B);
    cost_check(&tmp);
    cost_addition(&opt_cost, &tmp);

    // cost_zero();
    // gettimeofday(&tv_A,NULL);
    // g1g2_to_g3_pairing_affine(&test2,&P,&Q);
    // gettimeofday(&tv_B,NULL);
    // opt_affine_time+=timedifference_msec(tv_A,tv_B);
    // cost_check(&tmp);
    // cost_addition(&opt_affine_cost,&tmp);

    //if(fp12_cmp(&test1,&test2) != 0){
    //printf("pairing projective failed!\n\n");
    //printf("\n\n");
    //return 1;
    //}
  }
  cost_substruction(&FINALEXP_OPT_PROJECTIVE_COST, &opt_cost, &MILLER_OPT_PROJECTIVE_COST);
  //cost_substruction(&FINALEXP_OPT_AFFINE_COST, &opt_affine_cost, &MILLER_OPT_AFFINE_COST);

  printf("bls12 opt ate           : %.4f[ms]\n", opt_time / pairing);
  printf("bls12 opt ate (MILLER).          : %.4f[ms]\n", MILLER_OPT_PROJECTIVE / pairing);
  printf("bls12 opt ate (FINALEXP).        : %.4f[ms]\n", FINALEXP_OPT_PROJECTIVE / pairing);

  // printf("bls12 opt ate affine.                : %.4f[ms]\n",opt_affine_time/pairing);
  // printf("bls12 opt ate (MILLER).          : %.4f[ms]\n",MILLER_OPT_AFFINE/pairing);
  // printf("bls12 opt ate (FINALEXP).        : %.4f[ms]\n",FINALEXP_OPT_AFFINE/pairing);

#ifdef DEBUG_COST_A
  printf("*********bls12 opt ate montrick fp COST.********         \n");
  cost_printf("bls12 opt ate ", &opt_cost, pairing);
  cost_printf("bls12 opt ate (MILLER)", &MILLER_OPT_PROJECTIVE_COST, pairing);
  cost_printf("bls12 opt ate (FINALEXP)", &FINALEXP_OPT_PROJECTIVE_COST, pairing);
  printf("***************************************         \n");
// printf("*********bls12 opt ate GS fp COST.********         \n");
// cost_printf("bls12 opt ate ", &opt_affine_cost, pairing);
// cost_printf("bls12 opt ate (MILLER)", &MILLER_OPT_AFFINE_COST, pairing);
// cost_printf("bls12 opt ate (FINALEXP)", &FINALEXP_OPT_AFFINE_COST, pairing);
// printf("***************************************         \n");
#endif

  return 0;
}
void g3_test(int exp) {
  g1_t A_g1;
  g2_t A_g2;
  g3_t A_g3;
  fp12_t A_fp12;
  fr_t sca_fr;
  mpz_t sca;
  //initialize
  mpz_init(sca);
  fr_init(&sca_fr);
  g1_init(&A_g1);
  g2_init(&A_g2);
  g3_init(&A_g3);
  fp12_init(&A_fp12);
  //set
  g1_set_random(&A_g1, state);
  g2_set_random(&A_g2, state);
  g1g2_to_g3_pairing(&A_g3, &A_g1, &A_g2);

  //bench
  int i = 0;
  float scm_time = 0;
  cost tmp, scm_cost;
  struct timeval tv_A, tv_B;
  for (i = 0; i < exp; i++) {
    fr_set_random(&sca_fr, state);
    //faster type
    cost_zero();
    gettimeofday(&tv_A, NULL);
    g3_exp(&A_g3, &A_g3, &sca_fr);
    gettimeofday(&tv_B, NULL);
    scm_time += timedifference_msec(tv_A, tv_B);
    cost_check(&tmp);
    cost_addition(&scm_cost, &tmp);
  }
  printf("bls12 g3 exp.     : %.4f[ms]\n", scm_time / exp);

#ifdef DEBUG_COST_A
  printf("*********bls12 g2 scm fp COST.********         \n");
  cost_printf("bls12 g3 exp", &scm_cost, exp);
  printf("***************************************         \n");
#endif
  mpz_clear(sca);
}
