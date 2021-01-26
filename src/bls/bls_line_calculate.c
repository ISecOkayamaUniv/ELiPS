#include <ELiPS/bls_line_calculate.h>

void bls_f_ltp(fpm2_t *f, efpm2_t *p, efpm2_t *q, efpm2_t *t){
    fpm2_t tmp1,tmp2,tmp3;

    if(t->infinity){
        fpm2_set_ui(f,1);
    }
    else if(fpm2_cmp(&q->x,&t->x)==0){
        fpm2_sub(f,&q->x,&t->x);
    }
    else{
    fpm2_sub(&tmp1,&p->y,&t->y);
    fpm2_sub(&tmp2,&p->x,&t->x);
    fpm2_inv(&tmp3,&tmp2);
    fpm2_mul(&tmp2,&tmp1,&tmp3);
    fpm2_sub(&tmp3,&q->x,&p->x);
    fpm2_mul(&tmp2,&tmp2,&tmp3);
    fpm2_sub(&tmp1,&q->y,&p->y);
    fpm2_sub(f,&tmp1,&tmp2);
  }
}

void bls_f_ltt(fpm2_t *ANS, efpm2_t *q, efpm2_t *t){
  fpm2_t tmp1,tmp2,tmp3;

  if(t->infinity){
    fpm2_set_ui(ANS,1);
  }
  else if(fpm2_cmp_ui(&t->y,0)==0){
    fpm2_sub(ANS,&q->x,&t->x);
  }
  else{
    fpm2_mul(&tmp1,&t->x,&t->x);
    fpm2_add(&tmp2,&tmp1,&tmp1);
    fpm2_add(&tmp1,&tmp1,&tmp2);

    fpm2_add(&tmp2,&t->y,&t->y);
    fpm2_inv(&tmp3,&tmp2);
    fpm2_mul(&tmp2,&tmp1,&tmp3);

    fpm2_sub(&tmp1,&q->x,&t->x);
    fpm2_mul(&tmp3,&tmp1,&tmp2);
    fpm2_sub(&tmp1,&q->y,&t->y);
    fpm2_sub(ANS,&tmp1,&tmp3);
  }
}