#include <ELiPS/fp12cv.h>
//fp12_cvma
//can use p = 2,3(mod 5) && p = 5 (mod 7)
void fp12cv_init(fp12cv_t *a){
  fp4cv_init(&a->x0);
  fp4cv_init(&a->x1);
  fp4cv_init(&a->x2);
}

void fp12cv_printf(char *str, fp12cv_t *a){
  gmp_printf("%s(",str);
  fp4cv_printf("",&a->x0);
  gmp_printf(",");
  fp4cv_printf("",&a->x1);
  gmp_printf(",");
  fp4cv_printf("",&a->x2);
  gmp_printf(")");
}

void fp12cv_println(char *str, fp12cv_t *a){
  gmp_printf("%s(",str);
  fp4cv_printf("",&a->x0);
  gmp_printf(",");
  fp4cv_printf("",&a->x1);
  gmp_printf(",");
  fp4cv_printf("",&a->x2);
  gmp_printf(")\n");
}

void fp12cv_set(fp12cv_t *b, fp12cv_t *a){
  fp4cv_set(&b->x0,&a->x0);
  fp4cv_set(&b->x1,&a->x1);
  fp4cv_set(&b->x2,&a->x2);
}

void fp12cv_set_ui(fp12cv_t *a, unsigned long int b){
  fp4cv_t tmp1;
  if(b==0){
    fp4cv_init(&tmp1);
    fp4cv_set(&a->x0,&tmp1);
    fp4cv_set(&a->x1,&tmp1);
    fp4cv_set(&a->x2,&tmp1);
  }
  else{
    fp4cv_set_ui(&tmp1,b);
    fp4cv_set_neg(&tmp1,&tmp1);
    fp4cv_set(&a->x0,&tmp1);
    fp4cv_set(&a->x1,&tmp1);
    fp4cv_set(&a->x2,&tmp1);
  }
}

void fp12cv_set_random(fp12cv_t *a, gmp_randstate_t state){
  fp4cv_set_random(&a->x0,state);
  fp4cv_set_random(&a->x1,state);
  fp4cv_set_random(&a->x2,state);
}

void fp12cv_add(fp12cv_t *c, fp12cv_t *a, fp12cv_t *b){
  fp4cv_add(&c->x0,&a->x0,&b->x0);
  fp4cv_add(&c->x1,&a->x1,&b->x1);
  fp4cv_add(&c->x2,&a->x2,&b->x2);
}

void fp12cv_sub(fp12cv_t *c, fp12cv_t *a, fp12cv_t *b){
  fp4cv_sub(&c->x0,&a->x0,&b->x0);
  fp4cv_sub(&c->x1,&a->x1,&b->x1);
  fp4cv_sub(&c->x2,&a->x2,&b->x2);
}

void fp12cv_mul(fp12cv_t *c, fp12cv_t *a, fp12cv_t *b){
  fp4cv_t s1,s2,s3,tmp1,tmp2;
 
  fp4cv_sub(&tmp1,&a->x0,&a->x1);
  fp4cv_sub(&tmp2,&b->x1,&b->x0);
  fp4cv_mul(&s1,&tmp1,&tmp2);//s1 = (a0-a1)*(b1-b0)

  fp4cv_sub(&tmp1,&a->x1,&a->x2);
  fp4cv_sub(&tmp2,&b->x2,&b->x1);
  fp4cv_mul(&s2,&tmp1,&tmp2);//s2 = (a1-a2)*(b2-b1)

  fp4cv_sub(&tmp1,&a->x0,&a->x2);
  fp4cv_sub(&tmp2,&b->x2,&b->x0);
  fp4cv_mul(&s3,&tmp1,&tmp2);//s3 = (a0-a2)*(b2-b0)

  //x0
  fp4cv_add(&tmp1,&s1,&s2);
  fp4cv_mul(&tmp2,&a->x0,&b->x0);
  fp4cv_sub(&c->x0,&tmp1,&tmp2);//x0 = s1 + s2 - a0*b0
  //x1
  fp4cv_add(&tmp1,&s2,&s3);
  fp4cv_mul(&tmp2,&a->x1,&b->x1);
  fp4cv_sub(&c->x1,&tmp1,&tmp2);//x1 = s2 + s3 - a1*b1
  //x2
  fp4cv_add(&tmp1,&s3,&s1);
  fp4cv_mul(&tmp2,&a->x2,&b->x2);
  fp4cv_sub(&c->x2,&tmp1,&tmp2);//x2 = s3 + s1 - a2*b2
}

void fp12cv_sqr(fp12cv_t *b, fp12cv_t *a){
  fp4cv_t t1,t2,t3,t4,t5,tmp1,tmp2;

  //t1
  fp4cv_sub(&t1,&a->x1,&a->x0);//t1 = a1-a0
  //t2
  fp4cv_add(&tmp1,&t1,&a->x2);
  fp4cv_sqr(&t2,&tmp1);//t2 = (t1-a2)^2
  //t3
  fp4cv_sub(&tmp1,&a->x2,&a->x0);
  fp4cv_sub(&tmp2,&a->x1,&a->x2);
  fp4cv_mul(&t3,&tmp1,&tmp2);//t3 = (a2-a0)*(a1-a2)
  //t4
  fp4cv_sqr(&tmp1,&t1);
  fp4cv_sqr(&tmp2,&a->x1);
  fp4cv_add(&t4,&tmp1,&tmp2);//t4 = t1^2 + a1^2
  //t5
  fp4cv_add(&tmp1,&a->x0,&a->x2);
  fp4cv_mul(&tmp2,&t1,&tmp1);
  fp4cv_add(&t5,&tmp2,&t3);//t5 = t1*(a0+a2) + t3

  //x0
  fp4cv_sub(&b->x0,&t5,&t4);//x0 = t5 - t4
  //x1
  fp4cv_add(&tmp1,&t3,&t3);
  fp4cv_sub(&b->x1,&tmp1,&t4);//x1 = 2t3 - t4
  //x2
  fp4cv_sub(&b->x2,&t5,&t2);//x2 = t5 - t2
}
/*
void fp12cv_sqr_GS(fp12cv_t *b, fp12cv_t *a){
  fp4cv_t t0,t1,t2,t3,t4,t5,tmp1,tmp2;
  fp4cv_frobenius_times(&t3,&a->x0,2);
  fp4cv_frobenius_times(&t4,&a->x1,2);
  fp4cv_frobenius_times(&t5,&a->x2,2);

  fp4cv_mul(&t0,&a->x0,&a->x2);
  fp4cv_mul(&t1,&a->x0,&a->x1);
  fp4cv_mul(&t2,&a->x1,&a->x2);

  fp4cv_add(&t0,&t0,&t4);//t0 = (a0-a1)^2
  fp4cv_add(&t1,&t1,&t5);//t1 = (a1-a2)^2
  fp4cv_add(&t2,&t2,&t3);//t2 = (a0-a2)^2

  fp4cv_sqr(&t3,&a->x0);
  fp4cv_sqr(&t4,&a->x1);
  fp4cv_sqr(&t5,&a->x2);

  fp4cv_add(&tmp1,&t0,&t1);
  fp4cv_add(&tmp2,&tmp1,&t3);
  fp4cv_set_neg(&b->x0,&tmp2);
  fp4cv_add(&tmp1,&t2,&t1);
  fp4cv_add(&tmp2,&tmp1,&t4);
  fp4cv_set_neg(&b->x1,&tmp2);
  fp4cv_add(&tmp1,&t2,&t0);
  fp4cv_add(&tmp2,&tmp1,&t5);
  fp4cv_set_neg(&b->x2,&tmp2);
}
*/
void fp12cv_sqr_GS(fp12cv_t *b, fp12cv_t *a){
  fp4cv_t t0,t1,t2,t3,t4,t5,t6;
  fp4cv_add(&t0,&a->x0,&a->x1);
  fp4cv_add(&t0,&t0,&a->x2);
  fp4cv_mul(&t1,&t0,&a->x0);
  fp4cv_mul(&t2,&t0,&a->x1);
  fp4cv_mul(&t3,&t0,&a->x2);
  fp4cv_frobenius_times(&t4,&a->x0,2);
  fp4cv_frobenius_times(&t5,&a->x1,2);
  fp4cv_frobenius_times(&t6,&a->x2,2);
  fp4cv_add(&b->x0,&t1,&t5);
  fp4cv_add(&b->x0,&t6,&b->x0);
  fp4cv_set_neg(&b->x0,&b->x0);
  fp4cv_add(&b->x1,&t2,&t4);
  fp4cv_add(&b->x1,&t6,&b->x1);
  fp4cv_set_neg(&b->x1,&b->x1);
  fp4cv_add(&b->x2,&t3,&t4);
  fp4cv_add(&b->x2,&t5,&b->x2);
  fp4cv_set_neg(&b->x2,&b->x2);

}

void fp12cv_frobenius(fp12cv_t *b, fp12cv_t *a){
  fp4cv_t tmp1;
  //!!!Caution!!!
  //use tmp for input fp12cv_frobenius(x,x);

  //p = 5 (mod 7)
  //b0
  fp4cv_set(&tmp1,&a->x0);
  fp4cv_frobenius(&b->x0,&a->x2);//b0 = a2
  //b2
  fp4cv_frobenius(&b->x2,&a->x1);//b2 = a1
  //b1
  fp4cv_frobenius(&b->x1,&tmp1);//b1 = a0
}

void fp12cv_frobenius_times(fp12cv_t *b, fp12cv_t *a, unsigned long int number){
  int i;
  fp12cv_t tmp;
  fp12cv_set(&tmp,a);
  for(i=0;i<number;i++){
    fp12cv_frobenius(&tmp,&tmp);
  }
  fp12cv_set(b,&tmp);
}

void fp12cv_inv(fp12cv_t *b, fp12cv_t *a){
  fp4cv_t t1,t2,t3,u1,u2,u3,s1,s2,tmp1,tmp2;
  //T = A^p^4 * A^p^8 = (t1,t2,t3)
  //t1
  fp4cv_sub(&u1,&a->x2,&a->x0);
  fp4cv_sqr(&tmp1,&u1);
  fp4cv_mul(&tmp2,&a->x1,&a->x2);
  fp4cv_sub(&t1,&tmp1,&tmp2);//t1 = (a2-a0)^2 - a1*a2
  //t2
  fp4cv_sub(&u2,&a->x0,&a->x1);
  fp4cv_sqr(&tmp1,&u2);
  fp4cv_mul(&tmp2,&a->x0,&a->x2);
  fp4cv_sub(&t2,&tmp1,&tmp2);//t2 = (a0-a1)^2 - a0*a2
  //t3
  fp4cv_sub(&u3,&a->x2,&a->x1);
  fp4cv_sqr(&tmp1,&u3);
  fp4cv_mul(&tmp2,&a->x0,&a->x1);
  fp4cv_sub(&t3,&tmp1,&tmp2);//t3 = (a2-a1)^2 - a0*a1

  //S = A*T (=fp4)
  //S = (s,s,s)
  //s = (a0-a1)*(t2-t1) + (a1-a2)*(t3-t2) - a0*t1
  //s1
  fp4cv_sub(&tmp1,&t2,&t1);
  fp4cv_mul(&s1,&u2,&tmp1);//s1 = (a0-a1)*(t2-t1)
  //s2
  fp4cv_sub(&tmp1,&t3,&t2);
  fp4cv_mul(&s2,&u3,&tmp1);//s2 = -(a1-a2)*(t3-t2)
  //s
  fp4cv_sub(&tmp1,&s1,&s2);
  fp4cv_mul(&tmp2,&a->x0,&t1);
  fp4cv_sub(&tmp1,&tmp1,&tmp2);

  //inv
  fp4cv_inv(&tmp2,&tmp1);

  fp4cv_mul(&b->x0,&t1,&tmp2);
  fp4cv_mul(&b->x1,&t2,&tmp2);
  fp4cv_mul(&b->x2,&t3,&tmp2);
  fp4cv_set_neg(&b->x0,&b->x0);
  fp4cv_set_neg(&b->x1,&b->x1);
  fp4cv_set_neg(&b->x2,&b->x2);
  
}

void fp12cv_pow(fp12cv_t *b, fp12cv_t *a, mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp12cv_t tmp;

    fp12cv_set(&tmp,a);

    for(i=1;i<length; i++){
        fp12cv_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            fp12cv_mul(&tmp,a,&tmp);
        }
    }
    fp12cv_set(b,&tmp);
}

void fp12cv_sqrt(fp12cv_t *b, fp12cv_t *a){
    fp12cv_t x,y,t,k,n,tmp;
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);

    fp12cv_set_random(&n,state);
    while(fp12cv_legendre(&n)!=-1){
        fp12cv_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,12);//p^12-1
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    fp12cv_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp12cv_pow(&x,a,exp);
    fp12cv_mul(&tmp,&x,&x);
    fp12cv_mul(&k,&tmp,a);
    fp12cv_mul(&x,&x,a);
    while(fp12cv_cmp_ui(&k,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fp12cv_pow(&tmp,&k,exp);
        while(fp12cv_cmp_ui(&tmp,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fp12cv_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fp12cv_pow(&t,&y,result);
        fp12cv_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fp12cv_mul(&x,&x,&t);
        fp12cv_mul(&k,&k,&y);
    }
    fp12cv_set(b,&x);

    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}

int fp12cv_legendre(fp12cv_t *a){
    fp12cv_t tmp;
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime_z,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp12cv_pow(&tmp,a,exp);

    if(fp12cv_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }
    else{
        mpz_clear(exp);
        return -1;
    }
}

int fp12cv_cmp(fp12cv_t *a, fp12cv_t *b){
  if(fp4cv_cmp(&a->x0,&b->x0)==0 && fp4cv_cmp(&a->x1,&b->x1)==0 && fp4cv_cmp(&a->x2,&b->x2)==0){
    return 0;
  }
  else{
    return 1;
  }
}

int fp12cv_cmp_ui(fp12cv_t *a, unsigned long int b){
  fp12cv_t tmp1;
  fp12cv_set_ui(&tmp1,b);
  return fp12cv_cmp(a,&tmp1);
}

int fp12cv_cmp_zero(fp12cv_t *a){
  if(fp4cv_cmp_zero(&a->x0)==0 && fp4cv_cmp_zero(&a->x1)==0 && fp4cv_cmp_zero(&a->x2)==0){
    return 0;
  }
  return 1;
}

int fp12cv_cmp_one(fp12cv_t *a){
  fp4cv_t tmp1;
  fp4cv_set_ui(&tmp1,1);
  fp4cv_set_neg(&tmp1,&tmp1);
  if(fp4cv_cmp(&a->x0,&tmp1)==0 && fp4cv_cmp(&a->x1,&tmp1)==0 && fp4cv_cmp(&a->x2,&tmp1)==0){
    return 0;
  }
  return 1;
}

