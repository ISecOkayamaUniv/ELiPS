#include <ELiPS/efpm.h>
//efpm
void efpm_init(efpm_t *p){
  fpm_init(&p->x);
  fpm_init(&p->y);
  p->infinity = 0;
}
void efpm_set(efpm_t *q, efpm_t *p){
  fpm_set(&q->x,&p->x);
  fpm_set(&q->y,&p->y);
  q->infinity = p->infinity;
}
void efpm_printf(char *str, efpm_t *p){
  printf("%s",str);
  if(p->infinity==0){
    printf("x = ");
    fpm_printf("",&p->x);
    printf("y = ");
    fpm_printf("",&p->y);
  }
  else{
    printf("infinity\n");
  }
}
void efpm_rational_point(efpm_t *p){
  fpm_t b;
  fpm_t tmp1,tmp2;

  fpm_set_mpn(&b,curve_b);
  while(1){
    fpm_set_random(&p->x,state);
    fpm_mul(&tmp1,&p->x,&p->x);
    fpm_mul(&tmp2,&tmp1,&p->x);
    fpm_add(&tmp2,&tmp2,&b);
    if(fpm_legendre(&tmp2)==1){
      fpm_sqrt(&p->y,&tmp2);
      break;
    }
  }
}
void efpm_eca(efpm_t *p3, efpm_t *p1, efpm_t *p2){
  fpm_t x,y,ramda,tmp1,tmp2;

  //p1 = p2 ならefp2_DBLを実行
  if(fpm_cmp(&p1->x,&p2->x) && fpm_cmp(&p1->y,&p2->y)){
    efpm_ecd(p3,p1);
  }
  //無限遠点 + 無限遠点 = 無限遠点
  if(p1->infinity && p2->infinity){
    p3->infinity = 1;
  }
  //無限遠点 + p2 = p2
  else if(p1->infinity){
    efpm_set(p3,p2);
  }
  //無限遠点 + p1 = p1
  else if(p2->infinity){
    efpm_set(p3,p1);
  }
  //p1x = p2x なら無限遠点
  else if(fpm_cmp(&p1->x,&p2->x)){
    p3->infinity = 1;
  }
  else{
    fpm_sub(&tmp1,&p2->y,&p1->y);
    fpm_sub(&tmp2,&p2->x,&p1->x);
    fpm_inv(&tmp2,&tmp2);
    fpm_mul(&ramda,&tmp1,&tmp2);
    fpm_sqr(&tmp1,&ramda);
    fpm_sub(&x,&tmp1,&p1->x);
    fpm_sub(&x,&x,&p2->x);

    fpm_sub(&tmp1,&p1->x,&x);
    fpm_mul(&tmp1,&tmp1,&ramda);
    fpm_sub(&y,&tmp1,&p1->y);

    fpm_set(&p3->x,&x);
    fpm_set(&p3->y,&y);
  }
}
void efpm_ecd(efpm_t *p3, efpm_t *p1){
  fpm_t x,y,ramda,tmp1,tmp2;
  //p1が無限遠点ならば無限遠点
  if(p1->infinity){
    p3->infinity = 1;
  }
  //p1yが0ならば無限遠点
  else if(fpm_cmp_ui(&p1->y,0)==0){
    p3->infinity = 1;
  }
  else{
    //ramda = 3x1^2/2y1
    fpm_sqr(&ramda,&p1->x);
    fpm_set_ui(&tmp1,3);
    fpm_mul(&ramda,&ramda,&tmp1);
    fpm_set_ui(&tmp1,2);
    fpm_mul(&tmp2,&tmp1,&p1->y);
    fpm_inv(&tmp2,&tmp2);
    fpm_mul(&ramda,&ramda,&tmp2);
    //x = ramda^2 - x1 - x1
    fpm_sqr(&x,&ramda);
    fpm_sub(&x,&x,&p1->x);
    fpm_sub(&x,&x,&p1->x);
    //y = (x1 - x3)*ramda - y1
    fpm_sub(&tmp1,&p1->x,&x);
    fpm_mul(&tmp1,&tmp1,&ramda);
    fpm_sub(&y,&tmp1,&p1->y);
    //セット　efp_DBL(&p,&p)などの入力にも対応できるように
    fpm_set(&p3->x,&x);
    fpm_set(&p3->y,&y);
  }
}
void efpm_scm(efpm_t *q, efpm_t *p, mpz_t scalar){
  if(mpz_cmp_ui(scalar,0)==0){
      q->infinity = 1;
  }
  else if(mpz_cmp_ui(scalar,1)==0){
      efpm_set(q,p);
  }
  else{
    efpm_t tmp_p,next_p;
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);

    efpm_set(&tmp_p,p);
    efpm_set(&next_p,&tmp_p);
    for(i=1;i<length;i++){
      efpm_ecd(&next_p,&next_p);
      if(binary[i]=='1'){
        efpm_eca(&next_p,&next_p,&tmp_p);
      }
    }
    efpm_set(q,&next_p);
  }
}
