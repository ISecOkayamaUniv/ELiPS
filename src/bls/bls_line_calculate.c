#include <ELiPS/bls_line_calculate.h>

void bls_f_ltp(fpm2_t *f, efpm2_t *p, efpm2_t *q, efpm2_t *t) {
  fpm2_t tmp1, tmp2, tmp3;

  if (t->infinity) {
    fpm2_set_ui(f, 1);
  } else if (fpm2_cmp(&q->x, &t->x) == 0) {
    fpm2_sub(f, &q->x, &t->x);
  } else {
    fpm2_sub(&tmp1, &p->y, &t->y);
    fpm2_sub(&tmp2, &p->x, &t->x);
    fpm2_inv(&tmp3, &tmp2);
    fpm2_mul(&tmp2, &tmp1, &tmp3);
    fpm2_sub(&tmp3, &q->x, &p->x);
    fpm2_mul(&tmp2, &tmp2, &tmp3);
    fpm2_sub(&tmp1, &q->y, &p->y);
    fpm2_sub(f, &tmp1, &tmp2);
  }
}

void bls_f_ltt(fpm2_t *ANS, efpm2_t *q, efpm2_t *t) {
  fpm2_t tmp1, tmp2, tmp3;

  if (t->infinity) {
    fpm2_set_ui(ANS, 1);
  } else if (fpm2_cmp_ui(&t->y, 0) == 0) {
    fpm2_sub(ANS, &q->x, &t->x);
  } else {
    fpm2_mul(&tmp1, &t->x, &t->x);
    fpm2_add(&tmp2, &tmp1, &tmp1);
    fpm2_add(&tmp1, &tmp1, &tmp2);

    fpm2_add(&tmp2, &t->y, &t->y);
    fpm2_inv(&tmp3, &tmp2);
    fpm2_mul(&tmp2, &tmp1, &tmp3);

    fpm2_sub(&tmp1, &q->x, &t->x);
    fpm2_mul(&tmp3, &tmp1, &tmp2);
    fpm2_sub(&tmp1, &q->y, &t->y);
    fpm2_sub(ANS, &tmp1, &tmp3);
  }
}

void bls_2i3_ff_ltt(fpm2_t *f, efp2cv_t *P, efp2cv_t *T) {
  efp2cv_t Tmp_T;
  fp2cv_t t1, t2, t3;
  fp2cv_t tmp1, tmp2, tmp3;
  fpm2_t x, y, s1, s2, s3;
  efp2cv_init(&Tmp_T);
  efp2cv_set(&Tmp_T, T);

  fp2cv_add(&tmp1, &Tmp_T.y, &Tmp_T.y);
  fp2cv_inv(&tmp1, &tmp1);
  fp2cv_mul(&tmp2, &Tmp_T.x, &Tmp_T.x);
  fp2cv_add(&tmp3, &tmp2, &tmp2);
  fp2cv_add(&tmp2, &tmp3, &tmp2);
  fp2cv_mul(&t1, &tmp1, &tmp2);        //t1 = 3xt^2/2yt
  fp2cv_add(&t2, &Tmp_T.x, &Tmp_T.x);  //t2 = 2xt
  fp2cv_mul(&tmp1, &t1, &Tmp_T.x);
  fp2cv_sub(&t3, &tmp1, &Tmp_T.y);  //t3 = t1*xt - yt
  //T <- 2T
  fp2cv_mul(&tmp1, &t1, &t1);
  fp2cv_sub(&T->x, &tmp1, &t2);
  fp2cv_mul(&tmp1, &t1, &T->x);
  fp2cv_sub(&T->y, &t3, &tmp1);
  //ltt
  fpm2_set_fp2cv(&x, &P->x);
  fpm2_set_fp2cv(&y, &P->y);
  fpm2_set_fp2cv(&s1, &t1);
  fpm2_mul(&s2, &s1, &x);
  fpm2_mul(&s1, &s2, &twist_g_sxrt_inv);
  fpm2_set_fp2cv(&s3, &t3);
  fpm2_mul(&s2, &s3, &twist_g_sqrt_inv);
  fpm2_sub(&s3, &y, &s1);
  fpm2_add(&s3, &s3, &s2);
  //f <- f^2*ltt
  fpm2_mul(&s1, f, f);
  fpm2_mul(f, &s1, &s3);
}
void bls_2i3_ff_ltq(fpm2_t *f, efp2cv_t *P, efp2cv_t *Q, efp2cv_t *T) {
  fp2cv_t t1, t2, t3;
  fp2cv_t tmp1, tmp2, tmp3;
  efp2cv_t Tmp_T;
  fpm2_t x, y, s1, s2, s3;
  efp2cv_init(&Tmp_T);
  efp2cv_set(&Tmp_T, T);

  fp2cv_sub(&tmp1, &Q->x, &Tmp_T.x);
  fp2cv_inv(&tmp1, &tmp1);
  fp2cv_sub(&tmp2, &Q->y, &Tmp_T.y);
  fp2cv_mul(&t1, &tmp1, &tmp2);     //t1 = (yq-yt)/(xq-xt)
  fp2cv_add(&t2, &Q->x, &Tmp_T.x);  //t2 = xq + xt
  fp2cv_mul(&tmp1, &t1, &Tmp_T.x);
  fp2cv_sub(&t3, &tmp1, &Tmp_T.y);  //t3 = t1*xt - yt
  //T <- T + Q
  fp2cv_mul(&tmp1, &t1, &t1);
  fp2cv_sub(&T->x, &tmp1, &t2);
  fp2cv_mul(&tmp1, &t1, &T->x);
  fp2cv_sub(&T->y, &t3, &tmp1);
  //ltq
  fpm2_set_fp2cv(&x, &P->x);
  fpm2_set_fp2cv(&y, &P->y);
  fpm2_set_fp2cv(&s1, &t1);
  fpm2_mul(&s2, &s1, &x);
  fpm2_mul(&s1, &s2, &twist_g_sxrt_inv);
  fpm2_set_fp2cv(&s3, &t3);
  fpm2_mul(&s2, &s3, &twist_g_sqrt_inv);
  fpm2_sub(&s3, &y, &s1);
  fpm2_add(&s3, &s3, &s2);
  //f <- f*ltq
  fpm2_mul(f, f, &s3);
}
