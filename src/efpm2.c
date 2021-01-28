#include <ELiPS/efpm2.h>
//efpm2
void efpm2_init(efpm2_t *p){
  fpm2_init(&p->x);
  fpm2_init(&p->y);
  p->infinity = 0;
}
void efpm2_set(efpm2_t *q, efpm2_t *p){
  fpm2_set(&q->x,&p->x);
  fpm2_set(&q->y,&p->y);
  q->infinity = p->infinity;
}
void efpm2_printf(char *str, efpm2_t *p){
  printf("%s",str);
  if(p->infinity==0){
    printf("x = \n");
    fpm2_println("",&p->x);
    printf("y = \n");
    fpm2_println("",&p->y);
  }
  else{
    printf("infinity\n");
  }
}
void efpm2_println(char *str, efpm2_t *p){
  printf("%s",str);
  if(p->infinity==0){
    printf("x = \n");
    fpm2_println("",&p->x);
    printf("y = \n");
    fpm2_println("",&p->y);
  }
  else{
    printf("infinity\n");
  }
}
void efpm2_rational_point(efpm2_t *p){
  
  fpm2_t b;
  fpm2_t tmp1,tmp2;

  fpm2_set_mpn(&b,curve_b);
  while(1){
    fpm2_set_random(&p->x,state);
    fpm2_mul(&tmp1,&p->x,&p->x);
    fpm2_mul(&tmp2,&tmp1,&p->x);
    fpm2_add(&tmp2,&tmp2,&b);
    if(fpm2_legendre(&tmp2)==1){
      fpm2_sqrt(&p->y,&tmp2);
      break;
    }
  }
}
void efpm2_eca(efpm2_t *p3, efpm2_t *p1, efpm2_t *p2){
  fpm2_t x,y,ramda,tmp1,tmp2;

  //p1 = p2 ならefp2_DBLを実行
  if(fpm2_cmp(&p1->x,&p2->x)==0 && fpm2_cmp(&p1->y,&p2->y)==0){
    efpm2_ecd(p3,p1);
  }
  //無限遠点 + 無限遠点 = 無限遠点
  if(p1->infinity && p2->infinity){
    p3->infinity = 1;
  }
  //無限遠点 + p2 = p2
  else if(p1->infinity){
    efpm2_set(p3,p2);
  }
  //無限遠点 + p1 = p1
  else if(p2->infinity){
    efpm2_set(p3,p1);
  }
  //p1x = p2x なら無限遠点
  else if(fpm2_cmp(&p1->x,&p2->x)==0){
    p3->infinity = 1;
  }
  else{
    fpm2_sub(&tmp1,&p2->y,&p1->y);
    fpm2_sub(&tmp2,&p2->x,&p1->x);
    fpm2_inv(&tmp2,&tmp2);
    fpm2_mul(&ramda,&tmp1,&tmp2);
    fpm2_sqr(&tmp1,&ramda);
    fpm2_sub(&x,&tmp1,&p1->x);
    fpm2_sub(&x,&x,&p2->x);

    fpm2_sub(&tmp1,&p1->x,&x);
    fpm2_mul(&tmp1,&tmp1,&ramda);
    fpm2_sub(&y,&tmp1,&p1->y);

    fpm2_set(&p3->x,&x);
    fpm2_set(&p3->y,&y);
  }
}
void efpm2_ecd(efpm2_t *p3, efpm2_t *p1){
  fpm2_t x,y,ramda,tmp1,tmp2;
  //p1が無限遠点ならば無限遠点
  if(p1->infinity){
    p3->infinity = 1;
  }
  //p1yが0ならば無限遠点
  else if(fpm2_cmp_ui(&p1->y,0)==0){
    p3->infinity = 1;
  }
  else{
    //ramda = 3x1^2/2y1
    fpm2_sqr(&ramda,&p1->x);
    fpm2_set_ui(&tmp1,3);
    fpm2_mul(&ramda,&ramda,&tmp1);
    fpm2_set_ui(&tmp1,2);
    fpm2_mul(&tmp2,&tmp1,&p1->y);
    fpm2_inv(&tmp2,&tmp2);
    fpm2_mul(&ramda,&ramda,&tmp2);
    //x = ramda^2 - x1 - x1
    fpm2_sqr(&x,&ramda);
    fpm2_sub(&x,&x,&p1->x);
    fpm2_sub(&x,&x,&p1->x);
    //y = (x1 - x3)*ramda - y1
    fpm2_sub(&tmp1,&p1->x,&x);
    fpm2_mul(&tmp1,&tmp1,&ramda);
    fpm2_sub(&y,&tmp1,&p1->y);
    //セット　efp_DBL(&p,&p)などの入力にも対応できるように
    fpm2_set(&p3->x,&x);
    fpm2_set(&p3->y,&y);
  }
}
void efpm2_scm(efpm2_t *q, efpm2_t *p, mpz_t scalar){
  if(mpz_cmp_ui(scalar,0)==0){
      q->infinity = 1;
  }
  else if(mpz_cmp_ui(scalar,1)==0){
      efpm2_set(q,p);
  }
  else{
    efpm2_t tmp_p,next_p;
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);

    efpm2_set(&tmp_p,p);
    efpm2_set(&next_p,&tmp_p);
    for(i=1;i<length;i++){
      efpm2_ecd(&next_p,&next_p);
      if(binary[i]=='1'){
        efpm2_eca(&next_p,&next_p,&tmp_p);
      }
    }
    efpm2_set(q,&next_p);
  }
}
