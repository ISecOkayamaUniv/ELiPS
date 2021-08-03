#include <ELiPS/mpn.h>

void mpn_init(mp_limb_t *a, mp_size_t size) {
  mpn_zero(a, size);
}

void mpn_set_char(mp_limb_t *ans, mp_size_t mp_size, char *str) {
  unsigned long int i, sizeL;
  char *str_buf;
  mp_size_t size;

  sizeL = strlen(str);

  str_buf = (char *)malloc(sizeL * sizeof(char));

  for (i = 0; i < sizeL; i++) {
    str_buf[i] = str[i] - 48;
  }

  size = mpn_set_str(ans, (unsigned const char *)str_buf, sizeL, 10);
  for (i = size; i < mp_size; i++) {
    ans[i] = 0;
  }

  free(str_buf);
}
void mpn_set_ui(mp_limb_t *ans, mp_size_t size, unsigned long int ui) {
  unsigned long int i;

  ans[0] = ui;

  for (i = 1; i < size; i++) {
    ans[i] = 0;
  }
}

void mpn_set_mpz(mp_limb_t *ans, mpz_t a) {
  char *str;

  str = (char *)malloc(mpz_sizeinbase(a, 10) + 2);

  str = mpz_get_str(str, 10, a);
  mpn_set_char(ans, FPLIMB, str);

  free(str);
}

void mpn_mod(mp_limb_t *ans, mp_limb_t *a, mp_size_t size_a) {
  mp_limb_t dumy[size_a];
  mpn_tdiv_qr(dumy, ans, 0, a, size_a, prime, FPLIMB);
}

int mpn_cmp_ui(mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  mpn_set_ui(buf, size, ui);
  if (mpn_cmp(a, buf, size) == 0) {
    return 0;
  } else {
    return 1;
  }
}
void mpn_lshift_ext(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, long int L) {
  mp_limb_t tmp[size];

  mpn_copyd(tmp, a, size);
  while (L > 63) {
    mpn_lshift(tmp, tmp, size, 63);
    L = L - 63;
  }
  mpn_lshift(ans, tmp, size, L);
}
void mpn_add_ui(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  mp_limb_t buf[size];

  mpn_set_ui(buf, size, ui);

  mpn_add_n(ans, a, buf, size);
}
void mpn_sub_ui(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  mp_limb_t buf[size];

  mpn_set_ui(buf, size, ui);

  mpn_sub_n(ans, a, buf, size);
}
void mpn_mul_ui(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  mp_limb_t buf[size];

  mpn_set_ui(buf, size, ui);

  mpn_mul_n(ans, a, buf, size);
}
void mpn_invert(mp_limb_t *ANS, mp_limb_t *A, mp_limb_t *p) {
  mp_limb_t prime_tmp[FPLIMB], gp[FPLIMB], sp[FPLIMB], tmp[FPLIMB];
  mp_size_t buf_size;

  mpn_init(gp, FPLIMB);
  mpn_init(sp, FPLIMB);
  mpn_init(tmp, FPLIMB);
  mpn_init(prime_tmp, FPLIMB);

  mpn_copyd(prime_tmp, p, FPLIMB);

  mpn_add_n(buf, A, p, FPLIMB);
  mpn_gcdext(gp, sp, &buf_size, buf, FPLIMB, prime_tmp, FPLIMB);

  if (buf_size < 0) {
    mpn_sub_n(tmp, p, sp, FPLIMB);
  } else {
    mpn_copyd(tmp, sp, FPLIMB);
  }

  mpn_mod(ANS, tmp, FPLIMB);
}

void mpn_mulmod_montgomery(mp_limb_t *ANS, mp_size_t ANS_size, mp_limb_t *A, mp_size_t A_size, mp_limb_t *B, mp_size_t B_size) {
#ifdef DEBUG_COST_A
  cost_mul++;
  cost_mod++;
#endif

  static mp_limb_t T[FPLIMB2];
  mpn_zero(T, FPLIMB2);

  mpn_mul(T, A, A_size, B, B_size);
  for (int i = 0; i < FPLIMB; i++)
    T[i] = mpn_addmul_1(&T[i], prime, FPLIMB, T[i] * Ni_neg);

  mpn_add_n(ANS, T + FPLIMB, T, FPLIMB);
  if (mpn_cmp(ANS, prime, FPLIMB) != -1) mpn_sub_n(ANS, ANS, prime, FPLIMB);
}
//TODO: unuse
// void mpn_sqrmod_montgomery(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size){
//     #ifdef DEBUG_COST_A
//     cost_sqr++;
//     cost_mod++;
//     #endif

//     unsigned long int carry;
//     mp_limb_t r;
//     static mp_limb_t T[FPLIMB2];
//     static unsigned int index=0;
//     static unsigned long int c;
//     static int i,j;
//     mpn_zero(T,FPLIMB2);

//     mpn_sqr(T,A,A_size);
//     index=0;
//     for (i = 0; i < FPLIMB; i++,index++) {
//         r = (mp_limb_t)(T[index] * Ni_neg);
//         T[index] = mpn_addmul_1(T+index,prime,FPLIMB,r);
//     }
//     carry = mpn_add_n(ANS, T+FPLIMB, T, FPLIMB);
//     if (carry || (mpn_cmp(ANS, prime, FPLIMB) != -1)) {
//         carry = mpn_sub_n(ANS,ANS,prime,FPLIMB);
//     }
// }

void mpn_mod_montgomery(mp_limb_t *ANS, mp_size_t ANS_size, mp_limb_t *A, mp_size_t A_size) {
#ifdef DEBUG_COST_A
  cost_mod++;
#endif
  static mp_limb_t T[FPLIMB2];
  mpn_zero(T, FPLIMB2);

  mpn_copyd(T, A, A_size);
  for (int i = 0; i < FPLIMB; i++)
    T[i] = mpn_addmul_1(&T[i], prime, FPLIMB, T[i] * Ni_neg);

  mpn_add_n(ANS, T + FPLIMB, T, FPLIMB);
  while (mpn_cmp(ANS, prime, FPLIMB) != -1) mpn_sub_n(ANS, ANS, prime, FPLIMB);
}

void mpn_to_montgomery(mp_limb_t *ANS, mp_limb_t *A) {
#ifdef DEBUG_COST_A
  //cost_mod++;
  cost_mod_nomal++;
#endif
  static int i;
  static mp_limb_t tmp[FPLIMB2];
  mpn_zero(tmp, FPLIMB2);
  for (i = FPLIMB; i < FPLIMB2; i++) tmp[i] = A[i - FPLIMB];
  mpn_mod(ANS, tmp, FPLIMB2);
}
