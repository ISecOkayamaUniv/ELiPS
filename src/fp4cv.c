#include <ELiPS/fp4cv.h>
//fp4_cvma
//can use p=2,3(mod 5)
void fp4cv_init(fp4cv_t *a){
  fp_init(&a->x0);
  fp_init(&a->x1);
  fp_init(&a->x2);
  fp_init(&a->x3);
}

void fp4cv_printf(char *str,fp4cv_t *a){
  gmp_printf("%s(",str);
  fp_printf("",&a->x0);
  gmp_printf(",");
  fp_printf("",&a->x1);
  gmp_printf(",");
  fp_printf("",&a->x2);
  gmp_printf(",");
  fp_printf("",&a->x3);
  gmp_printf(")");
}

void fp4cv_println(char *str,fp4cv_t *a){
  gmp_printf("%s(",str);
  fp_printf("",&a->x0);
  gmp_printf(",");
  fp_printf("",&a->x1);
  gmp_printf(",");
  fp_printf("",&a->x2);
  gmp_printf(",");
  fp_printf("",&a->x3);
  gmp_printf(")\n");
}

void fpd4cv_println(char *str,fpd4cv_t *a){
  gmp_printf("%s(",str);
  fpd_printf("",&a->x0);
  gmp_printf(",");
  fpd_printf("",&a->x1);
  gmp_printf(",");
  fpd_printf("",&a->x2);
  gmp_printf(",");
  fpd_printf("",&a->x3);
  gmp_printf(")\n");
}

void fp4cv_set(fp4cv_t *b, fp4cv_t *a){
  fp_set(&b->x0,&a->x0);
  fp_set(&b->x1,&a->x1);
  fp_set(&b->x2,&a->x2);
  fp_set(&b->x3,&a->x3);
}

void fpd4cv_set(fpd4cv_t *b, fpd4cv_t *a){
  fpd_set(&b->x0,&a->x0);
  fpd_set(&b->x1,&a->x1);
  fpd_set(&b->x2,&a->x2);
  fpd_set(&b->x3,&a->x3);
}

void fp4cv_set_ui(fp4cv_t *a, unsigned long int b){
  fp_t tmp1;
  fp_set_ui(&tmp1,b);
  fp_set_neg(&tmp1,&tmp1);
  fp_set(&a->x0,&tmp1);
  fp_set(&a->x1,&tmp1);
  fp_set(&a->x2,&tmp1);
  fp_set(&a->x3,&tmp1);
}

void fp4cv_set_mpn(fp4cv_t *b,mp_limb_t *a){
  fp_t tmp1;
  fp_set_mpn(&tmp1,a);
  fp_set_neg(&tmp1,&tmp1);
  fp_set(&b->x0,&tmp1);
  fp_set(&b->x1,&tmp1);
  fp_set(&b->x2,&tmp1);
  fp_set(&b->x3,&tmp1);
}

void fp4cv_set_neg(fp4cv_t *b, fp4cv_t *a){
  fp_set_neg(&b->x0,&a->x0);
  fp_set_neg(&b->x1,&a->x1);
  fp_set_neg(&b->x2,&a->x2);
  fp_set_neg(&b->x3,&a->x3);
}

void fp4cv_set_random(fp4cv_t *a,gmp_randstate_t state){
    fp_set_random(&a->x0,state);
    fp_set_random(&a->x1,state);
    fp_set_random(&a->x2,state);
    fp_set_random(&a->x3,state);
}

void fp4cv_add(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b){
  fp_add(&c->x0,&a->x0,&b->x0);
  fp_add(&c->x1,&a->x1,&b->x1);
  fp_add(&c->x2,&a->x2,&b->x2);
  fp_add(&c->x3,&a->x3,&b->x3);
}

void fp4cv_add_nonmod_single(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b){
  fp_add_nonmod_single(&c->x0,&a->x0,&b->x0);
  fp_add_nonmod_single(&c->x1,&a->x1,&b->x1);
  fp_add_nonmod_single(&c->x2,&a->x2,&b->x2);
  fp_add_nonmod_single(&c->x3,&a->x3,&b->x3);
}

void fp4cv_add_nonmod_double(fpd4cv_t *c, fpd4cv_t *a, fpd4cv_t *b){
  fp_add_nonmod_double(&c->x0,&a->x0,&b->x0);
  fp_add_nonmod_double(&c->x1,&a->x1,&b->x1);
  fp_add_nonmod_double(&c->x2,&a->x2,&b->x2);
  fp_add_nonmod_double(&c->x3,&a->x3,&b->x3);
}

void fp4cv_add_ui(fp4cv_t *c, fp4cv_t *a, unsigned long int b){
  fp_t tmp1;
  fp_set_ui(&tmp1,b);
  fp_set_neg(&tmp1,&tmp1);
  fp_add(&c->x0,&a->x0,&tmp1);
  fp_add(&c->x1,&a->x1,&tmp1);
  fp_add(&c->x2,&a->x2,&tmp1);
  fp_add(&c->x3,&a->x3,&tmp1);
}

void fp4cv_sub(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b){
  fp_sub(&c->x0,&a->x0,&b->x0);
  fp_sub(&c->x1,&a->x1,&b->x1);
  fp_sub(&c->x2,&a->x2,&b->x2);
  fp_sub(&c->x3,&a->x3,&b->x3);
}

void fp4cv_sub_nonmod_single(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b){
  fp_sub_nonmod_single(&c->x0,&a->x0,&b->x0);
  fp_sub_nonmod_single(&c->x1,&a->x1,&b->x1);
  fp_sub_nonmod_single(&c->x2,&a->x2,&b->x2);
  fp_sub_nonmod_single(&c->x3,&a->x3,&b->x3);
}

void fp4cv_sub_nonmod_double(fpd4cv_t *c, fpd4cv_t *a, fpd4cv_t *b){
  fp_sub_nonmod_double(&c->x0,&a->x0,&b->x0);
  fp_sub_nonmod_double(&c->x1,&a->x1,&b->x1);
  fp_sub_nonmod_double(&c->x2,&a->x2,&b->x2);
  fp_sub_nonmod_double(&c->x3,&a->x3,&b->x3);
}

void fp4cv_sub_ui(fp4cv_t *c, fp4cv_t *a, unsigned long int b){
  fp_t tmp1;
  fp_set_ui(&tmp1,b);
  fp_set_neg(&tmp1,&tmp1);
  fp_sub(&c->x0,&a->x0,&tmp1);
  fp_sub(&c->x1,&a->x1,&tmp1);
  fp_sub(&c->x2,&a->x2,&tmp1);
  fp_sub(&c->x3,&a->x3,&tmp1);
}

void fp4cv_mul(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b){
  fp_t t1,t2,t3,t4,t5,t6,t7,t8,u1,tmp1,tmp2;
  //prepare
  fp_sub(&t1,&a->x0,&a->x2);
  fp_sub(&t2,&a->x1,&a->x3);
  fp_sub(&t3,&b->x0,&b->x2);
  fp_sub(&t4,&b->x1,&b->x3);
  fp_mul(&t5,&t2,&t4);
  fp_mul(&t6,&t1,&t3);
  fp_sub(&tmp1,&a->x0,&a->x1);
  fp_sub(&tmp2,&b->x0,&b->x1);
  fp_mul(&t7,&tmp1,&tmp2);
  fp_sub(&tmp1,&a->x2,&a->x3);
  fp_sub(&tmp2,&b->x2,&b->x3);
  fp_mul(&t8,&tmp1,&tmp2);
  fp_add(&tmp1,&t1,&t2);
  fp_add(&tmp2,&t3,&t4);
  fp_mul(&u1,&tmp1,&tmp2);
  fp_sub(&u1,&u1,&t5);
  fp_sub(&u1,&u1,&t6);
  fp_add(&u1,&u1,&t7);
  fp_add(&u1,&u1,&t8);
  //x0
  fp_mul(&tmp1,&a->x0,&b->x0);
  fp_sub(&tmp2,&u1,&t5);
  fp_sub(&c->x0,&tmp2,&tmp1);
  //x1
  fp_mul(&tmp1,&a->x1,&b->x1);
  fp_sub(&tmp2,&u1,&t8);
  fp_sub(&c->x1,&tmp2,&tmp1);
  //x2
  fp_mul(&tmp1,&a->x2,&b->x2);
  fp_sub(&tmp2,&u1,&t7);
  fp_sub(&c->x2,&tmp2,&tmp1);
  //x3
  fp_mul(&tmp1,&a->x3,&b->x3);
  fp_sub(&tmp2,&u1,&t6);
  fp_sub(&c->x3,&tmp2,&tmp1);
}

void fp4cv_mul_lazy(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b){
  fp_t t1,t2,t3,t4,tmp1,tmp2;
  fpd_t t5,t6,t7,t8,u1,tmp3,tmp4;
  //prepare
  fp_sub(&t1,&a->x0,&a->x2);
  fp_sub(&t2,&a->x1,&a->x3);
  fp_sub(&t3,&b->x0,&b->x2);
  fp_sub(&t4,&b->x1,&b->x3);

  fp_mul_nonmod(&t5,&t2,&t4);
  fp_mul_nonmod(&t6,&t1,&t3);

  fp_sub(&tmp1,&a->x0,&a->x1);
  fp_sub(&tmp2,&b->x0,&b->x1);
  fp_mul_nonmod(&t7,&tmp1,&tmp2);

  fp_sub(&tmp1,&a->x2,&a->x3);
  fp_sub(&tmp2,&b->x2,&b->x3);
  fp_mul_nonmod(&t8,&tmp1,&tmp2);

  fp_add_nonmod_single(&tmp1,&t1,&t2);
  fp_add_nonmod_single(&tmp2,&t3,&t4);
  fp_mul_nonmod(&u1,&tmp1,&tmp2);

  fp_sub_nonmod_double(&u1,&u1,&t5);
  fp_sub_nonmod_double(&u1,&u1,&t6);
  fp_add_nonmod_double(&u1,&u1,&t7);
  fp_add_nonmod_double(&u1,&u1,&t8);
  //x0
  fp_mul_nonmod(&tmp3,&a->x0,&b->x0);
  fp_sub_nonmod_double(&tmp4,&u1,&t5);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  fp_mod(&c->x0,tmp4.x0,FPLIMB2);
  //x1
  fp_mul_nonmod(&tmp3,&a->x1,&b->x1);
  fp_sub_nonmod_double(&tmp4,&u1,&t8);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  fp_mod(&c->x1,tmp4.x0,FPLIMB2);
  //x2
  fp_mul_nonmod(&tmp3,&a->x2,&b->x2);
  fp_sub_nonmod_double(&tmp4,&u1,&t7);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  fp_mod(&c->x2,tmp4.x0,FPLIMB2);
  //x3
  fp_mul_nonmod(&tmp3,&a->x3,&b->x3);
  fp_sub_nonmod_double(&tmp4,&u1,&t6);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  fp_mod(&c->x3,tmp4.x0,FPLIMB2);
}

/*
void fp4cv_mul_lazy_montgomery(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b){
  fp_t t1,t2,t3,t4,tmp1,tmp2;
  fpd_t t5,t6,t7,t8,u1,tmp3,tmp4;
  //prepare
  fp_sub(&t1,&a->x0,&a->x2);
  fp_sub(&t2,&a->x1,&a->x3);
  fp_sub(&t3,&b->x0,&b->x2);
  fp_sub(&t4,&b->x1,&b->x3);

  fp_mul_nonmod(&t5,&t2,&t4);
  fp_mul_nonmod(&t6,&t1,&t3);

  fp_sub(&tmp1,&a->x0,&a->x1);
  fp_sub(&tmp2,&b->x0,&b->x1);
  fp_mul_nonmod(&t7,&tmp1,&tmp2);

  fp_sub(&tmp1,&a->x2,&a->x3);
  fp_sub(&tmp2,&b->x2,&b->x3);
  fp_mul_nonmod(&t8,&tmp1,&tmp2);

  fp_add_nonmod_single(&tmp1,&t1,&t2);
  fp_add_nonmod_single(&tmp2,&t3,&t4);
  fp_mul_nonmod(&u1,&tmp1,&tmp2);

  fp_sub_nonmod_double(&u1,&u1,&t5);
  fp_sub_nonmod_double(&u1,&u1,&t6);
  fp_add_nonmod_double(&u1,&u1,&t7);
  fp_add_nonmod_double(&u1,&u1,&t8);
  //x0
  fp_mul_nonmod(&tmp3,&a->x0,&b->x0);
  fp_sub_nonmod_double(&tmp4,&u1,&t5);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  mpn_mod_montgomery(c->x0.x0,FPLIMB,tmp4.x0,FPLIMB2);
  //x1
  fp_mul_nonmod(&tmp3,&a->x1,&b->x1);
  fp_sub_nonmod_double(&tmp4,&u1,&t8);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  mpn_mod_montgomery(c->x1.x0,FPLIMB,tmp4.x0,FPLIMB2);
  //x2
  fp_mul_nonmod(&tmp3,&a->x2,&b->x2);
  fp_sub_nonmod_double(&tmp4,&u1,&t7);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  mpn_mod_montgomery(c->x2.x0,FPLIMB,tmp4.x0,FPLIMB2);
  //x3
  fp_mul_nonmod(&tmp3,&a->x3,&b->x3);
  fp_sub_nonmod_double(&tmp4,&u1,&t6);
  fp_sub_nonmod_double(&tmp4,&tmp4,&tmp3);
  mpn_mod_montgomery(c->x3.x0,FPLIMB,tmp4.x0,FPLIMB2);
}
*/
void fp4cv_sqr(fp4cv_t *b, fp4cv_t *a){
  fp_t t1,t2,t3,t4,t5,t6,tmp1,tmp2;

  fp_sub(&t1,&a->x0,&a->x2);
  fp_sub(&t2,&a->x0,&a->x1);
  fp_sub(&t3,&a->x2,&a->x3);
  fp_sub(&t4,&a->x1,&a->x3);
  fp_mul(&tmp1,&t2,&t3);
  fp_add(&t5,&tmp1,&tmp1);
  fp_mul(&tmp1,&t1,&t4);
  fp_add(&t6,&tmp1,&tmp1);
  //x0
  fp_add(&tmp1,&a->x0,&t1);
  fp_mul(&tmp2,&a->x2,&tmp1);
  fp_sub(&b->x0,&t5,&tmp2);
  //x1
  fp_sub(&tmp1,&a->x1,&t2);
  fp_mul(&tmp2,&a->x0,&tmp1);
  fp_sub(&b->x1,&t6,&tmp2);
  //x2
  fp_add(&tmp1,&a->x2,&t3);
  fp_mul(&tmp2,&a->x3,&tmp1);
  fp_sub(&b->x2,&t6,&tmp2);
  //x3
  fp_sub(&tmp1,&a->x3,&t4);
  fp_mul(&tmp2,&a->x1,&tmp1);
  fp_sub(&b->x3,&t5,&tmp2);
}

void fp4cv_sqr_lazy(fp4cv_t *b, fp4cv_t *a){
  fp_t t1,t2,t3,t4,tmp1,tmp2;
  fpd_t t5,t6,tmp3,tmp4;

  fp_sub(&t1,&a->x0,&a->x2);
  fp_sub(&t2,&a->x0,&a->x1);
  fp_sub(&t3,&a->x2,&a->x3);
  fp_sub(&t4,&a->x1,&a->x3);

  fp_mul_nonmod(&tmp3,&t2,&t3);
  fp_add_nonmod_double(&t5,&tmp3,&tmp3);
  fp_mul_nonmod(&tmp4,&t1,&t4);
  fp_add_nonmod_double(&t6,&tmp4,&tmp4);
  //x0
  fp_add_nonmod_single(&tmp1,&a->x0,&t1);
  fp_mul_nonmod(&tmp3,&a->x2,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t5,&tmp3);
  fp_mod(&b->x0,tmp4.x0,FPLIMB2);
  //x1
  fp_sub_nonmod_single(&tmp1,&a->x1,&t2);
  fp_mul_nonmod(&tmp3,&a->x0,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t6,&tmp3);
  fp_mod(&b->x1,tmp4.x0,FPLIMB2);
  //x2
  fp_add_nonmod_single(&tmp1,&a->x2,&t3);
  fp_mul_nonmod(&tmp3,&a->x3,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t6,&tmp3);
  fp_mod(&b->x2,tmp4.x0,FPLIMB2);
  //x3
  fp_sub_nonmod_single(&tmp1,&a->x3,&t4);
  fp_mul_nonmod(&tmp3,&a->x1,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t5,&tmp3);
  fp_mod(&b->x3,tmp4.x0,FPLIMB2);
}
/*
void fp4cv_sqr_lazy_montgomery(fp4cv_t *b, fp4cv_t *a){
  fp_t t1,t2,t3,t4,tmp1,tmp2;
  fpd_t t5,t6,tmp3,tmp4;

  fp_sub(&t1,&a->x0,&a->x2);
  fp_sub(&t2,&a->x0,&a->x1);
  fp_sub(&t3,&a->x2,&a->x3);
  fp_sub(&t4,&a->x1,&a->x3);

  fp_mul_nonmod(&tmp3,&t2,&t3);
  fp_add_nonmod_double(&t5,&tmp3,&tmp3);
  fp_mul_nonmod(&tmp4,&t1,&t4);
  fp_add_nonmod_double(&t6,&tmp4,&tmp4);
  
  //x0
  fp_add_nonmod_single(&tmp1,&a->x0,&t1);
  fp_mul_nonmod(&tmp3,&a->x2,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t5,&tmp3);
  mpn_mod_montgomery(b->x0.x0,FPLIMB,tmp4.x0,FPLIMB2);
  //x1
  fp_sub_nonmod_single(&tmp1,&a->x1,&t2);
  fp_mul_nonmod(&tmp3,&a->x0,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t6,&tmp3);
  mpn_mod_montgomery(b->x1.x0,FPLIMB,tmp4.x0,FPLIMB2);
  //x2
  fp_add_nonmod_single(&tmp1,&a->x2,&t3);
  fp_mul_nonmod(&tmp3,&a->x3,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t6,&tmp3);
  mpn_mod_montgomery(b->x2.x0,FPLIMB,tmp4.x0,FPLIMB2);
  //x3
  fp_sub_nonmod_single(&tmp1,&a->x3,&t4);
  fp_mul_nonmod(&tmp3,&a->x1,&tmp1);
  fp_sub_nonmod_double(&tmp4,&t5,&tmp3);
  mpn_mod_montgomery(b->x3.x0,FPLIMB,tmp4.x0,FPLIMB2);
}
*/
void fp4cv_frobenius(fp4cv_t *b, fp4cv_t *a){
  fp4cv_t tmp1;
  //!!!don't change fp_set(&b->x0,&a->x2);!!!
  //use tmp for input fp4cv_frobenius(x,x);
  fp_set(&tmp1.x0,&a->x2);
  fp_set(&tmp1.x1,&a->x0);
  fp_set(&tmp1.x2,&a->x3);
  fp_set(&tmp1.x3,&a->x1);

  fp4cv_set(b,&tmp1);
}

void fp4cv_inv(fp4cv_t *b, fp4cv_t *a){
  fp4cv_t ans,tmp1;
  fp_t t1,t2,t3,t4,t5,t6,t7,t8;
  fp_t b0,b1,i0,i1;

  fp_sub(&t1,&a->x0,&a->x2);//t1 = a0-a2
  fp_sub(&t2,&a->x1,&a->x3);//t2 = a1-a3
  fp_sub(&t3,&a->x0,&a->x1);//t3 = a0-a1
  fp_sub(&t4,&a->x2,&a->x3);//t4 = a2-a3
  fp_mul(&t5,&t1,&t2);//t5 = (a0-a2)*(a1-a3)
  fp_mul(&t6,&t3,&t4);//t6 = (a0-a1)*(a2-a3)
  fp_add(&t7,&t1,&t2);//t7 = a0+a1-a2-a3
  fp_sqr(&t8,&t7);//t8 = (a0+a1-a2-a3)^2

  //b1
  fp_sub(&t1,&t5,&t8);//t1 = -(a0+a1-a2-a3)^2 + (a0-a2)*(a1-a3)
  fp_add(&t2,&t1,&t5);//t2 = -(a0+a1-a2-a3)^2 + 2(a0-a2)*(a1-a3)
  fp_sub(&t3,&t2,&t6);//t3 = -(a0+a1-a2-a3)^2 + 2(a0-a2)*(a1-a3) - (a0-a1)*(a2-a3)
  fp_mul(&t4,&a->x1,&a->x2);//t4 = a1*a2
  fp_sub(&b1,&t3,&t4);//b1 = -(a0+a1-a2-a3)^2 + 2(a0-a2)*(a1-a3) - (a0-a1)*(a2-a3) - a1*a2
  //b0
  fp_add(&t1,&t3,&t5);//t1 = -(a0+a1-a2-a3)^2 + 3(a0-a2)*(a1-a3) - (a0-a1)*(a2-a3)
  fp_sub(&t2,&t1,&t6);//t2 = -(a0+a1-a2-a3)^2 + 3(a0-a2)*(a1-a3) - 2(a0-a1)*(a2-a3)
  fp_mul(&t3,&a->x0,&a->x3);//t3 = a0*a3
  fp_sub(&b0,&t2,&t3);//b0 = -(a0+a1-a2-a3)^2 + 3(a0-a2)*(a1-a3) - 2(a0-a1)*(a2-a3) - a0*a3

  //s
  fp_sub(&t1,&b0,&b1);//t1 = b0-b1
  fp_sqr(&t2,&t1);//t2 = (b0-b1)^2
  fp_mul(&t3,&b0,&b1);//t3 = b0*b1
  fp_sub(&t4,&t2,&t3);//t4 = (b0-b1)^2 - b0*b1 = s
  //(s,s,s,s) = p-s
  fp_set_neg(&t5,&t4);
  fp_inv(&t6,&t5);//t6 = (p-s)^-1
  

  //B   = (b0,b1,b1,b0)
  //B^p = (b1,b0,b0,b1);
  //B^-1 = s^-1*B^p = (p-s)^-1 * (b1,b0,b0,b1) = (i0,i1,i1,i0)
  fp_mul(&i0,&b1,&t6);
  fp_mul(&i1,&b0,&t6);
  //A^-1 = B^-1*A^p2
  //A^-1 = (i0,i1,i1,i0)*(a3,a2,a1,a0)

  fp_sub(&t1,&i0,&i1);//t1 = i0-i1
  fp_sub(&t2,&a->x1,&a->x0);//t2 = a1-a0
  fp_sub(&t3,&a->x0,&a->x3);//t3 = a0-a3
  fp_sub(&t4,&a->x2,&a->x1);//t4 = a2-a1
  fp_mul(&t5,&t1,&t2);//t5 = (i0-i1)*(a1-a0)
  fp_mul(&t6,&t1,&t3);//t6 = (i0-i1)*(a0-a3)
  fp_mul(&t7,&t1,&t4);//t7 = (i0-i1)*(a2-a1)
  //x0
  fp_add(&t8,&t5,&t7);//t8 = (i0-i1)*(a1-a0) + (i0-i1)*(a2-a1)
  fp_mul(&t1,&i0,&a->x3);
  fp_sub(&ans.x0,&t8,&t1);
  //x1
  fp_mul(&t1,&i1,&a->x2);
  fp_sub(&ans.x1,&t5,&t1);
  //x2
  fp_add(&t2,&t8,&t6);
  fp_mul(&t1,&i1,&a->x1);
  fp_sub(&ans.x2,&t2,&t1);
  //x3
  fp_add(&t2,&t5,&t6);
  fp_mul(&t1,&i0,&a->x0);
  fp_sub(&ans.x3,&t2,&t1);
  //set
  fp_set(&b->x0,&ans.x0);
  fp_set(&b->x1,&ans.x1);
  fp_set(&b->x2,&ans.x2);
  fp_set(&b->x3,&ans.x3);
}

void fp4cv_pow(fp4cv_t *b,fp4cv_t *a,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp4cv_t tmp;

    fp4cv_set(&tmp,a);

    for(i=1;i<length; i++){
        fp4cv_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            fp4cv_mul(&tmp,a,&tmp);
        }
    }
    fp4cv_set(b,&tmp);
}

void fp4cv_sqrt(fp4cv_t *b,fp4cv_t *a){
    fp4cv_t x,y,t,k,n,tmp;
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);

    fp4cv_set_random(&n,state);
    while(fp4cv_legendre(&n)!=-1){
        fp4cv_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,4);//p^4-1
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    fp4cv_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp4cv_pow(&x,a,exp);
    fp4cv_mul(&tmp,&x,&x);
    fp4cv_mul(&k,&tmp,a);
    fp4cv_mul(&x,&x,a);
    while(fp4cv_cmp_ui(&k,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fp4cv_pow(&tmp,&k,exp);
        while(fp4cv_cmp_ui(&tmp,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fp4cv_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fp4cv_pow(&t,&y,result);
        fp4cv_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fp4cv_mul(&x,&x,&t);
        fp4cv_mul(&k,&k,&y);
    }
    fp4cv_set(b,&x);

    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}

int fp4cv_legendre(fp4cv_t *a){
    fp4cv_t tmp;
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime_z,4);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp4cv_pow(&tmp,a,exp);

    if(fp4cv_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }
    else{
        mpz_clear(exp);
        return -1;
    }
}

int fp4cv_cmp(fp4cv_t *a, fp4cv_t *b){
  if(fp_cmp(&a->x0,&b->x0)==0 && fp_cmp(&a->x1,&b->x1)==0 && fp_cmp(&a->x2,&b->x2)==0 && fp_cmp(&a->x3,&b->x3)==0){
    return 0;
  }
  else{
    return 1;
  }
}

int fp4cv_cmp_ui(fp4cv_t *a, unsigned long int b){
  fp4cv_t tmp1;
  fp4cv_set_ui(&tmp1,b);
  return fp4cv_cmp(a,&tmp1);
}

int fp4cv_cmp_zero(fp4cv_t *a){
  if(fp_cmp_zero(&a->x0)==0 && fp_cmp_zero(&a->x1)==0 && fp_cmp_zero(&a->x2)==0 && fp_cmp_zero(&a->x3)==0){
      return 0;
  }
  return 1;
}

int fp4cv_cmp_one(fp4cv_t *a){
  fp_t tmp1;
  fp_set_ui(&tmp1,1);
  fp_set_neg(&tmp1,&tmp1);
  if(fp_cmp(&a->x0,&tmp1)==0 && fp_cmp(&a->x1,&tmp1)==0 && fp_cmp(&a->x2,&tmp1)==0 && fp_cmp(&a->x3,&tmp1)==0){
    return 0;
  }
  return 1;
}
