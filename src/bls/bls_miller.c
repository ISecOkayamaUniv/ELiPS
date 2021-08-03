#include <ELiPS/bls_miller.h>

void bls_miller_for_optate_basic(fpm2_t *ANS, efpm2_t *p, efpm2_t *q) {
  int i, length;
  fpm2_t f, tmp1, tmp2;
  efpm2_t p1, q1, q1_neg, t;
  //init
  efpm2_init(&p1);
  efpm2_init(&q1);
  efpm2_init(&q1_neg);
  efpm2_init(&t);
  //set
  //length=(int)mpz_sizeinbase(X_z,2);
  //printf("log2(X_z) = %d\n",length);
  efpm2_set(&p1, p);
  efpm2_set(&q1, q);
  fpm2_set(&q1_neg.x, &q1.x);
  fpm2_set_neg(&q1_neg.y, &q1.y);

  //step1
  fpm2_set_ui(&f, 1);  //f = 1
  efpm2_set(&t, &q1);  //t = q

  for (i = bls_X_length - 1; i >= 0; i--) {
    fpm2_mul(&tmp1, &f, &f);
    bls_f_ltt(&tmp2, &p1, &t);
    fpm2_mul(&f, &tmp1, &tmp2);
    efpm2_ecd(&t, &t);
    if (bls_X_binary[i] == 1) {
      bls_f_ltp(&tmp1, &q1, &p1, &t);
      fpm2_mul(&f, &f, &tmp1);
      efpm2_eca(&t, &t, &q1);
    } else if (bls_X_binary[i] == -1) {
      bls_f_ltp(&tmp1, &q1_neg, &p1, &t);
      fpm2_mul(&f, &f, &tmp1);
      efpm2_eca(&t, &t, &q1_neg);
    }
  }

  fpm2_set(ANS, &f);
}
void bls_2i3_miller_for_optate(fpm2_t *ANS, efpm2_t *P, efpm2_t *Q) {
  int i, length;
  fpm2_t f;
  efp2cv_t mapped_P, mapped_Q, mapped_Q_neg, mapped_T;
  //init
  efp2cv_init(&mapped_P);
  efp2cv_init(&mapped_Q);
  efp2cv_init(&mapped_Q_neg);
  efp2cv_init(&mapped_T);
  //twist and set
  efpm2_to_efp2cv_for_g1(&mapped_P, P);
  efpm2_to_efp2cv_for_g2(&mapped_Q, Q);
  fp2cv_set(&mapped_Q_neg.x, &mapped_Q.x);
  fp2cv_set_neg(&mapped_Q_neg.y, &mapped_Q.y);
  efp2cv_set(&mapped_T, &mapped_Q);

  fpm2_set_ui(&f, 1);
  //miller
  for (i = bls_X_length - 1; i >= 0; i--) {
    switch (bls_X_binary[i]) {
      case 0:
        bls_2i3_ff_ltt(&f, &mapped_P, &mapped_T);
        break;
      case 1:
        bls_2i3_ff_ltt(&f, &mapped_P, &mapped_T);
        bls_2i3_ff_ltq(&f, &mapped_P, &mapped_Q, &mapped_T);
        break;
      case -1:
        bls_2i3_ff_ltt(&f, &mapped_P, &mapped_T);
        bls_2i3_ff_ltq(&f, &mapped_P, &mapped_Q_neg, &mapped_T);
      default:
        break;
    }
  }
  fpm2_set(ANS, &f);
}