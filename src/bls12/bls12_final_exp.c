#include <ELiPS/bls12_final_exp.h>
void bls12_fp12_pow_X_compress_montrick(fp12_t *ANS, fp12_t *A) {
  int i = 0, p = 0;
  fp12_t tmp;
  fp12_t recover_list[4];

  fp12_init(&tmp);
  fp12_set(&tmp, A);
  for (i = 0; i < bls12_X_length; i++) {
    switch (bls12_X_binary[i]) {
      case 0:
        fp12_sqr_compressed_lazy_montgomery(&tmp, &tmp);
        break;
      case 1:
        fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
        fp12_sqr_compressed_lazy_montgomery(&tmp, &tmp);
        p++;
        break;
      case -1:
        fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
        fp12_frobenius_map_p6_montgomery(&recover_list[p], &recover_list[p]);
        fp12_sqr_compressed_lazy_montgomery(&tmp, &tmp);
        p++;
        break;
      default:
        break;
    }
  }
  if (bls12_X_binary[bls12_X_length] == 1) {
    fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
    p++;
  } else {
    fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
    fp12_frobenius_map_p6_montgomery(&recover_list[p], &recover_list[p]);
    p++;
  }

  if (p == 3) {
    fp12_sqr_recover_g1_montrick_montgomery_ham3(&recover_list[0], &recover_list[1], &recover_list[2]);
  } else if (p == 4) {
    fp12_sqr_recover_g1_montrick_montgomery(&recover_list[0], &recover_list[1], &recover_list[2], &recover_list[3]);
  } else {
    printf("not support: hamming 5 or more\n");
    return;
  }
  fp12_sqr_recover_g0_lazy_montgomery(&tmp, &recover_list[0]);
  for (i = 1; i < p; i++) {
    fp12_sqr_recover_g0_lazy_montgomery(&recover_list[i], &recover_list[i]);
    fp12_mul_lazy_montgomery(&tmp, &tmp, &recover_list[i]);
  }
  fp12_set(ANS, &tmp);
}
void bls12_fp12_pow_X2_compress_montrick(fp12_t *ANS, fp12_t *A) {
  int i = 0, p = 0;
  fp12_t tmp;
  fp12_t recover_list[10];

  fp12_init(&tmp);
  fp12_set(&tmp, A);
  for (i = 0; i < bls12_X2_length; i++) {
    switch (bls12_X2_binary[i]) {
      case 0:
        fp12_sqr_compressed_lazy_montgomery(&tmp, &tmp);
        break;
      case 1:
        fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
        fp12_sqr_compressed_lazy_montgomery(&tmp, &tmp);
        p++;
        break;
      case -1:
        fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
        fp12_frobenius_map_p6_montgomery(&recover_list[p], &recover_list[p]);
        fp12_sqr_compressed_lazy_montgomery(&tmp, &tmp);
        p++;
        break;
      default:
        break;
    }
  }
  if (bls12_X2_binary[bls12_X2_length] == 1) {
    fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
    p++;
  } else {
    fp12_sqr_recover_g1_noninv(&recover_list[p], &tmp);
    fp12_frobenius_map_p6_montgomery(&recover_list[p], &recover_list[p]);
    p++;
  }
  if (p == 3) {
    fp12_sqr_recover_g1_montrick_montgomery_ham3(&recover_list[0], &recover_list[1], &recover_list[2]);
  } else if (p == 4) {
    fp12_sqr_recover_g1_montrick_montgomery(&recover_list[0], &recover_list[1], &recover_list[2], &recover_list[3]);
  } else {
    printf("not support: hamming 5 or more\n");
    return;
  }
  fp12_sqr_recover_g0_lazy_montgomery(&tmp, &recover_list[0]);
  for (i = 1; i < p; i++) {
    fp12_sqr_recover_g0_lazy_montgomery(&recover_list[i], &recover_list[i]);
    fp12_mul_lazy_montgomery(&tmp, &tmp, &recover_list[i]);
  }
  fp12_set(ANS, &tmp);
}
void bls12_fp12_pow_X(fp12_t *ANS, fp12_t *A) {
  int i;
  fp12_t tmp, A_inv;
  fp12_init(&tmp);
  fp12_init(&A_inv);
  fp12_frobenius_map_p6(&A_inv, A);

  fp12_set(&tmp, A);
  for (i = bls12_X_length - 1; i >= 0; i--) {
    switch (bls12_X_binary[i]) {
      case 0:
        fp12_sqr_cyclotomic(&tmp, &tmp);
        break;
      case 1:
        fp12_sqr_cyclotomic(&tmp, &tmp);
        fp12_mul(&tmp, &tmp, A);
        break;
      case -1:
        fp12_sqr_cyclotomic(&tmp, &tmp);
        fp12_mul(&tmp, &tmp, &A_inv);
        break;
      default:
        break;
    }
  }
  fp12_set(ANS, &tmp);
}
void bls12_fp12_pow_X2(fp12_t *ANS, fp12_t *A) {
  int i;
  fp12_t tmp, A_inv;
  fp12_init(&tmp);
  fp12_init(&A_inv);
  fp12_frobenius_map_p6(&A_inv, A);

  fp12_set(&tmp, A);
  for (i = bls12_X2_length - 1; i >= 0; i--) {
    switch (bls12_X2_binary[i]) {
      case 0:
        fp12_sqr_cyclotomic(&tmp, &tmp);
        break;
      case 1:
        fp12_sqr_cyclotomic(&tmp, &tmp);
        fp12_mul(&tmp, &tmp, A);
        break;
      case -1:
        fp12_sqr_cyclotomic(&tmp, &tmp);
        fp12_mul(&tmp, &tmp, &A_inv);
        break;
      default:
        break;
    }
  }
  fp12_set(ANS, &tmp);
}
void bls12_final_exp_optimal(fp12_t *ANS, fp12_t *A) {
  fp12_t tmp, t0, t1, t2, t3, t4, t5, test;
  fp12_init(&tmp);
  fp12_init(&t0);
  fp12_init(&t1);
  fp12_init(&t2);
  fp12_init(&t3);
  fp12_init(&t5);
  fp12_init(&t4);
  fp12_init(&test);
  mpz_t positive_X, positive_X2;
  mpz_init(positive_X);
  mpz_init(positive_X2);

  //EASY PART
  mpz_set(positive_X, X_z);
  mpz_tdiv_q_ui(positive_X2, positive_X, 2);

  //f←f^(p^6)*f^-1
  fp12_frobenius_map_p6(&t0, A);  //f^(p^6)
  fp12_inv(&t1, A);               //f^-1
  fp12_mul(&tmp, &t0, &t1);       //f^(p^6)*f^-1

  //f←f^(p^2)*f
  fp12_frobenius_map_p2(&t0, &tmp);  //f^(p^2)
  fp12_mul(&tmp, &t0, &tmp);         //f^(p^2)*f

  //HARD PART
  fp12_sqr_cyclotomic(&t0, &tmp);
  bls12_fp12_pow_X(&t1, &t0);

  bls12_fp12_pow_X2(&t2, &t1);       //t2:=t1^(u2);
  fp12_frobenius_map_p6(&t3, &tmp);  //t3:=f^(-1);
  fp12_mul(&t1, &t3, &t1);           //t1:=t3*t1;
  fp12_frobenius_map_p6(&t1, &t1);   //t1:=t1^(-1);
  fp12_mul(&t1, &t1, &t2);           //t1:=t1*t2;

  bls12_fp12_pow_X(&t2, &t1);       //t2:=t1^(u);
  bls12_fp12_pow_X(&t3, &t2);       //t3:=t2^(u);
  fp12_frobenius_map_p6(&t1, &t1);  //t1:=t1^(-1);

  fp12_mul(&t3, &t1, &t3);          //t3:=t1*t3;
  fp12_frobenius_map_p6(&t1, &t1);  //t1:=t1^(-1);
  fp12_frobenius_map_p3(&t1, &t1);  //t1:=t1^(p^3);
  fp12_frobenius_map_p2(&t2, &t2);  //t2:=t2^(p^2);

  fp12_mul(&t1, &t1, &t2);     //t1:=t1*t2;
  bls12_fp12_pow_X(&t2, &t3);  //t2:=t3^(u);
  fp12_mul(&t2, &t2, &t0);     //t2:=t2*t0;
  fp12_mul(&t2, &t2, &tmp);    //t2:=t2*f;
  fp12_mul(&t1, &t1, &t2);     //t1:=t1*t2;

  fp12_frobenius_map_p1(&t2, &t3);  //t2:=t3^p;
  fp12_mul(ANS, &t1, &t2);          //t1:=t1*t2;

  mpz_clear(positive_X);
  mpz_clear(positive_X2);
}
void bls12_pow_hardpart_schoolbook_montgomery(fp12_t *ANS, fp12_t *A) {
  mpz_t power, tmp;
  mpz_init(power);
  mpz_init(tmp);
  fp12_t B;
  fp12_init(&B);

  mpz_pow_ui(power, prime_z, 4);
  mpz_pow_ui(tmp, prime_z, 2);
  mpz_sub(power, power, tmp);
  mpz_add_ui(power, power, 1);
  mpz_div(power, power, order_z);
  fp12_mod_montgomery(&B, A);
  fp12_pow(ANS, &B, power);

  mpz_clear(power);
  mpz_clear(tmp);
}
void bls12_pow_hardpart_compress_montgomery(fp12_t *ANS, fp12_t *A) {
  fp12_t t0, t1, t2, t3;
  //fp12_init(&tmp);
  fp12_init(&t0);
  fp12_init(&t1);
  fp12_init(&t2);
  fp12_init(&t3);
  // fp12_init(&t5);
  // fp12_init(&t4);
  // fp12_init(&test);
  // fp12_init(&At);
  //HARDPART
  fp12_sqr_GS_lazy_montgomery(&t0, A);
  bls12_fp12_pow_X_compress_montrick(&t1, &t0);

  //lamda3
  bls12_fp12_pow_X2_compress_montrick(&t2, &t1);  //t2:=t1^(u2);
  fp12_frobenius_map_p6_montgomery(&t3, A);       //t3:=f^(-1);
  fp12_mul_lazy_montgomery(&t1, &t3, &t1);        //t1:=t3*t1;
  fp12_frobenius_map_p6_montgomery(&t1, &t1);     //t1:=t1^(-1);
  fp12_mul_lazy_montgomery(&t1, &t1, &t2);        //t1:=t1*t2;

  bls12_fp12_pow_X_compress_montrick(&t2, &t1);  //t2:=t1^(u);
  bls12_fp12_pow_X_compress_montrick(&t3, &t2);  //t3:=t2^(u);

  //(lamda3)^(-1)
  fp12_frobenius_map_p6_montgomery(&t1, &t1);  //t1:=t1^(-1);

  fp12_mul_lazy_montgomery(&t3, &t1, &t3);  //t3:=t1*t3;

  fp12_frobenius_map_p6_montgomery(&t1, &t1);  //t1:=t1^(-1); //not need ??

  fp12_frobenius_map_p3_lazy_montgomery(&t1, &t1);  //t1:=t1^(p^3);
  fp12_frobenius_map_p2_montgomery(&t2, &t2);       //t2:=t2^(p^2);

  fp12_mul_lazy_montgomery(&t1, &t1, &t2);       //t1:=t1*t2;
  bls12_fp12_pow_X_compress_montrick(&t2, &t3);  //t2:=t3^(u);
  fp12_mul_lazy_montgomery(&t2, &t2, &t0);       //t2:=t2*t0;
  fp12_mul_lazy_montgomery(&t2, &t2, A);         //t2:=t2*f;
  fp12_mul_lazy_montgomery(&t1, &t1, &t2);       //t1:=t1*t2;

  fp12_frobenius_map_p1_lazy_montgomery(&t2, &t3);  //t2:=t3^p;
  fp12_mul_lazy_montgomery(ANS, &t1, &t2);          //t1:=t1*t2;
  fp12_mod_montgomery(ANS, ANS);
}
void bls12_final_exp_optimal_compress(fp12_t *ANS, fp12_t *A) {
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

  fp12_to_montgomery(&At, A);

  //EASY PART
  //mpz_set(positive_X,X_z);
  //mpz_tdiv_q_ui(positive_X2,positive_X,2);

  //f←f^(p^6)*f^-1
  fp12_frobenius_map_p6_montgomery(&t0, &At);  //f^(p^6)
  fp12_inv_lazy_montgomery(&t1, &At);          //f^-1
  fp12_mul_lazy_montgomery(&tmp, &t0, &t1);    //f^(p^6)*f^-1

  //f←f^(p^2)*f
  fp12_frobenius_map_p2_montgomery(&t0, &tmp);  //f^(p^2)
  fp12_mul_lazy_montgomery(&tmp, &t0, &tmp);    //f^(p^2)*f

  bls12_pow_hardpart_compress_montgomery(ANS, &tmp);
  //bls12_pow_hardpart_schoolbook_montgomery(ANS,&tmp);
}
/*---------------------------------------using cvma---------------------------------------*/
void bls12_final_exp_optimal_cvma(fp12_t *ANS, fp12_t *A) {
  fpm_t a1;
  fp12cv_t t1, t2;
  //!!!preparete

  //transformation fp12 -> fp12cv
  //fp12_to_fpm(&a1,A);
  //fpm_to_fp12cv(&t1,&a1);
  fp12_to_fp12cv(&t1, A);
  bls12_pow_easypart_cvma(&t2, &t1);
  bls12_pow_hardpart_cvma(&t1, &t2);
  //transformation fp12cv -> fp12
  fp12cv_to_fp12(ANS, &t1);
  //fp12cv_to_fpm(&a1,&t1);
  //fpm_to_fp12(ANS,&a1);
}
void bls12_pow_easypart_cvma(fp12cv_t *ANS, fp12cv_t *A) {
  fp12cv_t t1, t2, t3;
  fp12cv_frobenius_times(&t1, A, 6);
  fp12cv_inv(&t2, A);
  fp12cv_mul(&t3, &t1, &t2);
  fp12cv_frobenius_times(&t1, &t3, 2);
  fp12cv_mul(ANS, &t1, &t3);
}
void bls12_pow_hardpart_cvma(fp12cv_t *ANS, fp12cv_t *A) {
  fp12cv_t t0, t1, t2, t3, t4;
  fp12cv_sqr_GS(&t0, A);
  //fp12cv_sqr(&t0,A);
  bls12_fp12cv_pow_X(&t1, &t0);
  bls12_fp12cv_pow_X2(&t2, &t1);
  fp12cv_frobenius_times(&t3, A, 6);
  fp12cv_mul(&t1, &t1, &t3);
  fp12cv_frobenius_times(&t1, &t1, 6);
  fp12cv_mul(&t1, &t1, &t2);     //t1 = x^2-2x+1 = ramda3
  bls12_fp12cv_pow_X(&t2, &t1);  //t2 = ramda3^x = ramda2
  bls12_fp12cv_pow_X(&t3, &t2);
  fp12cv_frobenius_times(&t4, &t1, 6);  //t4 = -ramda3
  fp12cv_mul(&t3, &t3, &t4);            //t3 = ramda1
  bls12_fp12cv_pow_X(&t4, &t3);
  fp12cv_mul(&t4, &t4, &t0);
  fp12cv_mul(&t4, &t4, A);  //t4 = ramda0

  fp12cv_frobenius_times(&t1, &t1, 3);
  fp12cv_frobenius_times(&t2, &t2, 2);
  fp12cv_frobenius_times(&t3, &t3, 1);

  fp12cv_mul(&t1, &t1, &t2);
  fp12cv_mul(&t1, &t1, &t3);
  fp12cv_mul(ANS, &t1, &t4);
}
void bls12_fp12cv_pow_X(fp12cv_t *ANS, fp12cv_t *A) {
  int i;
  fp12cv_t tmp, A_inv;
  fp12cv_init(&tmp);
  fp12cv_init(&A_inv);
  fp12cv_frobenius_times(&A_inv, A, 6);

  fp12cv_set(&tmp, A);
  for (i = bls12_X_length - 1; i >= 0; i--) {
    switch (bls12_X_binary[i]) {
      case 0:
        //fp12cv_sqr(&tmp,&tmp);
        fp12cv_sqr_GS(&tmp, &tmp);
        break;
      case 1:
        //fp12cv_sqr(&tmp,&tmp);
        fp12cv_sqr_GS(&tmp, &tmp);
        fp12cv_mul(&tmp, &tmp, A);
        break;
      case -1:
        //fp12cv_sqr(&tmp,&tmp);
        fp12cv_sqr_GS(&tmp, &tmp);
        fp12cv_mul(&tmp, &tmp, &A_inv);
        break;
      default:
        break;
    }
  }
  fp12cv_set(ANS, &tmp);
}
void bls12_fp12cv_pow_X2(fp12cv_t *ANS, fp12cv_t *A) {
  int i;
  fp12cv_t tmp, A_inv;
  fp12cv_init(&tmp);
  fp12cv_init(&A_inv);
  fp12cv_frobenius_times(&A_inv, A, 6);

  fp12cv_set(&tmp, A);
  for (i = bls12_X2_length - 1; i >= 0; i--) {
    switch (bls12_X2_binary[i]) {
      case 0:
        //fp12cv_sqr(&tmp,&tmp);
        fp12cv_sqr_GS(&tmp, &tmp);
        break;
      case 1:
        //fp12cv_sqr(&tmp,&tmp);
        fp12cv_sqr_GS(&tmp, &tmp);
        fp12cv_mul(&tmp, &tmp, A);
        break;
      case -1:
        //fp12cv_sqr(&tmp,&tmp);
        fp12cv_sqr_GS(&tmp, &tmp);
        fp12cv_mul(&tmp, &tmp, &A_inv);
        break;
      default:
        break;
    }
  }
  fp12cv_set(ANS, &tmp);
}