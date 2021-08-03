#include <ELiPS/efp2cv.h>
//efp2cv
void efp2cv_init(efp2cv_t *p) {
  fp2cv_init(&p->x);
  fp2cv_init(&p->y);
  p->infinity = 0;
}
void efp2cv_set(efp2cv_t *q, efp2cv_t *p) {
  fp2cv_set(&q->x, &p->x);
  fp2cv_set(&q->y, &p->y);
  q->infinity = p->infinity;
}
void efp2cv_printf(char *str, efp2cv_t *p) {
  printf("%s", str);
  if (p->infinity == 0) {
    printf("x = \n");
    fp2cv_println("", &p->x);
    printf("y = \n");
    fp2cv_println("", &p->y);
  } else {
    printf("infinity\n");
  }
}
void efp2cv_println(char *str, efp2cv_t *p) {
  printf("%s", str);
  if (p->infinity == 0) {
    printf("x = \n");
    fp2cv_println("", &p->x);
    printf("y = \n");
    fp2cv_println("", &p->y);
  } else {
    printf("infinity\n");
  }
}
void efp2cv_rational_point(efp2cv_t *p) {
  fp2cv_t b;
  fp2cv_t tmp1, tmp2;

  fp2cv_set_mpn(&b, curve_b);
  while (1) {
    fp2cv_set_random(&p->x, state);
    fp2cv_mul(&tmp1, &p->x, &p->x);
    fp2cv_mul(&tmp2, &tmp1, &p->x);
    fp2cv_add(&tmp2, &tmp2, &b);
    if (fp2cv_legendre(&tmp2) == 1) {
      fp2cv_sqrt(&p->y, &tmp2);
      break;
    }
  }
}
void efp2cv_eca(efp2cv_t *p3, efp2cv_t *p1, efp2cv_t *p2) {
  fp2cv_t x, y, ramda, tmp1, tmp2;

  //p1 = p2 ならefp2_DBLを実行
  if (fp2cv_cmp(&p1->x, &p2->x) == 0 && fp2cv_cmp(&p1->y, &p2->y) == 0) {
    efp2cv_ecd(p3, p1);
  }
  //無限遠点 + 無限遠点 = 無限遠点
  if (p1->infinity && p2->infinity) {
    p3->infinity = 1;
  }
  //無限遠点 + p2 = p2
  else if (p1->infinity) {
    efp2cv_set(p3, p2);
  }
  //無限遠点 + p1 = p1
  else if (p2->infinity) {
    efp2cv_set(p3, p1);
  }
  //p1x = p2x なら無限遠点
  else if (fp2cv_cmp(&p1->x, &p2->x) == 0) {
    p3->infinity = 1;
  } else {
    fp2cv_sub(&tmp1, &p2->y, &p1->y);
    fp2cv_sub(&tmp2, &p2->x, &p1->x);
    fp2cv_inv(&tmp2, &tmp2);
    fp2cv_mul(&ramda, &tmp1, &tmp2);
    fp2cv_sqr(&tmp1, &ramda);
    fp2cv_sub(&x, &tmp1, &p1->x);
    fp2cv_sub(&x, &x, &p2->x);

    fp2cv_sub(&tmp1, &p1->x, &x);
    fp2cv_mul(&tmp1, &tmp1, &ramda);
    fp2cv_sub(&y, &tmp1, &p1->y);

    fp2cv_set(&p3->x, &x);
    fp2cv_set(&p3->y, &y);
  }
}
void efp2cv_ecd(efp2cv_t *p3, efp2cv_t *p1) {
  fp2cv_t x, y, ramda, tmp1, tmp2;
  //p1が無限遠点ならば無限遠点
  if (p1->infinity) {
    p3->infinity = 1;
  }
  //p1yが0ならば無限遠点
  else if (fp2cv_cmp_ui(&p1->y, 0) == 0) {
    p3->infinity = 1;
  } else {
    //ramda = 3x1^2/2y1
    fp2cv_sqr(&ramda, &p1->x);
    fp2cv_set_ui(&tmp1, 3);
    fp2cv_mul(&ramda, &ramda, &tmp1);
    fp2cv_set_ui(&tmp1, 2);
    fp2cv_mul(&tmp2, &tmp1, &p1->y);
    fp2cv_inv(&tmp2, &tmp2);
    fp2cv_mul(&ramda, &ramda, &tmp2);
    //x = ramda^2 - x1 - x1
    fp2cv_sqr(&x, &ramda);
    fp2cv_sub(&x, &x, &p1->x);
    fp2cv_sub(&x, &x, &p1->x);
    //y = (x1 - x3)*ramda - y1
    fp2cv_sub(&tmp1, &p1->x, &x);
    fp2cv_mul(&tmp1, &tmp1, &ramda);
    fp2cv_sub(&y, &tmp1, &p1->y);
    //セット　efp_DBL(&p,&p)などの入力にも対応できるように
    fp2cv_set(&p3->x, &x);
    fp2cv_set(&p3->y, &y);
  }
}
void efp2cv_scm(efp2cv_t *q, efp2cv_t *p, mpz_t scalar) {
  if (mpz_cmp_ui(scalar, 0) == 0) {
    q->infinity = 1;
  } else if (mpz_cmp_ui(scalar, 1) == 0) {
    efp2cv_set(q, p);
  } else {
    efp2cv_t tmp_p, next_p;
    int i, length;
    length = (int)mpz_sizeinbase(scalar, 2);
    char binary[length + 1];
    mpz_get_str(binary, 2, scalar);

    efp2cv_set(&tmp_p, p);
    efp2cv_set(&next_p, &tmp_p);
    for (i = 1; i < length; i++) {
      efp2cv_ecd(&next_p, &next_p);
      if (binary[i] == '1') {
        efp2cv_eca(&next_p, &next_p, &tmp_p);
      }
    }
    efp2cv_set(q, &next_p);
  }
}
