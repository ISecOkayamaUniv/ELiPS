#include <ELiPS/fpm.h>

void fpm_make_cvma(){
  int flag = 0;
  int flag_random_h = 1;
  int r;
  int h;
  int e;
  mpz_t tmp1;
  mpz_init(tmp1);
  //
  //
  for(h=0;;h++){
    r = h*DEGREE_EXTENTION_FIELD_M + 1;
    if(mpz_cmp_ui(prime_z,r)!=0){
      mpz_set_ui(tmp1,r);
      if(mpz_probab_prime_p(tmp1,50)){
        //flag
        flag = 1;
        //check gcd
        for(e=1;;e++){
          mpz_pow_ui(tmp1,prime_z,e);
          mpz_tdiv_r_ui(tmp1,tmp1,r);
          if(mpz_cmp_ui(tmp1,1)==0){
            mpz_set_ui(tmp1,h*DEGREE_EXTENTION_FIELD_M/e);
            mpz_gcd_ui(tmp1,tmp1,DEGREE_EXTENTION_FIELD_M);
            if(mpz_cmp_ui(tmp1,1)!=0)
              flag = 0;
            break;
          }
        }
        if(flag==1 && flag_random_h==1)
          break;
      }
    }
  }
  //set
  r_degree_m = r;
  h_degree_m = h;
  mpz_clear(tmp1);
}

void fpm_preparate_mul(){
  int i,j,k,l;
  int x1,x2,x3,x4;//代入変数
  int flag = 0;
  int d;
  int r = r_degree_m;
  int h = h_degree_m;
  int m = DEGREE_EXTENTION_FIELD_M;
  int e[SIZE_ALL];//サイズの修正
  int count;
  mpz_t tmp1;
  //init
  for(i=0;i<SIZE_ALL;i++){
    e[i] = 0;
    n_ijk_degree_m[i] = 0;
  }
  mpz_init(tmp1);
  //prime_zreparation step
  //step1
  for(d=1;d<r;d++){
    x1 = 1;
    for(i=0;i<h;i++){
      x1 = x1*d;
      x1 = x1%r;
      if(x1==1){
        if(i==h-1)
          flag = 1;
        else
          flag = 0;
        break;
      }
    }
    if(flag)
      break;
  }
  //step2
  e[0] = m;
  //step3,4

  for(i=0;i<=m-1;i++){
    for(k=0;k<=h-1;k++){
      mpz_pow_ui(tmp1,prime_z,i);
      mpz_tdiv_r_ui(tmp1,tmp1,r);
      x1 = mpz_get_ui(tmp1);

      x2 = 1;
      for(l=0;l<k;l++){
        x2 = x2*d;
        x2 = x2%r;
      }
      x3 = x1*x2;
      x4 = x3%r;
      e[x4] = i;
    }
  }
  //step5,6,7
  count = 0;
  for(i=0;i<=m-2;i++){
    mpz_pow_ui(tmp1,prime_z,i);
    mpz_tdiv_r_ui(tmp1,tmp1,r);
    x1 = mpz_get_ui(tmp1);//x = p^i (mod r)
    for(j=i+1;j<=m-1;j++){
      mpz_pow_ui(tmp1,prime_z,j);
      mpz_tdiv_r_ui(tmp1,tmp1,r);
      x2 = mpz_get_ui(tmp1);//y = p^j (mod r)
      for(k=0;k<=h-1;k++){
        //z = d^k
        x3 = 1;
        for(l=0;l<k;l++){
          x3 = x3*d;
          x3 = x3%r;
        }
        x4 = x2*x3;
        x4 = x1 + x4;
        x4 = x4%r;
        n_ijk_degree_m[count] = e[x4];
        count++;
      }
    }
  }
  //clear
  mpz_clear(tmp1);
}

void fpm_print_paramerter(){
  printf("\n-----  1->%d use cvma type<h,m>  -----\n",DEGREE_EXTENTION_FIELD_M);
  printf("h = %d r = %d\n",h_degree_m,r_degree_m);
}

void fpm_build(){
  fpm_make_cvma();
  fpm_preparate_mul();
  fpm_print_paramerter();
}

void fpm_init(fpm_t *a){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_init(&a->x[i]);
  }
}

void fpm_printf(char *str,fpm_t *a){
  int i;
  gmp_printf("%s(",str);
  fp_printf("",&a->x[0]);
  for(i=1;i<DEGREE_EXTENTION_FIELD_M;i++){
    gmp_printf(",");
    fp_printf("",&a->x[i]);
  }
  gmp_printf(")");
}

void fpm_println(char *str,fpm_t *a){
  int i;
  gmp_printf("%s",str);
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_println("",&a->x[i]);
  }
}

void fpm_set(fpm_t *b, fpm_t *a){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_set(&b->x[i],&a->x[i]);
  }
}

void fpm_set_ui(fpm_t *a, unsigned long int b){
  int i;
  fp_t tmp1;

  if(b==0){
    fpm_init(a);
  }
  else{
    for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
      fp_set_ui(&tmp1,b);
      fp_set_neg(&a->x[i],&tmp1);
    }
  }
}

void fpm_set_mpn(fpm_t *a, mp_limb_t *b){
  fp_t tmp1;
  fp_set_mpn(&tmp1,b);
  fpm_set_fp(a,&tmp1);
}

void fpm_set_fp(fpm_t *a, fp_t *b){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_set_neg(&a->x[i],b);
  }
}

void fpm_set_neg(fpm_t *b, fpm_t *a){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_set_neg(&b->x[i],&a->x[i]);
  }
}

void fpm_set_random(fpm_t *a,gmp_randstate_t state){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_set_random(&a->x[i],state);
  }
}

int fpm_cmp(fpm_t *a, fpm_t *b){
  int flag = 0;
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    if(fp_cmp(&a->x[i],&b->x[i])!=0){
      flag = 1;
    }
  }
  return flag;
}
int fpm_cmp_ui(fpm_t *a, unsigned long int b){
  int flag = 0;
  fpm_t tmp1;

  fpm_set_ui(&tmp1,b);
  flag = fpm_cmp(a,&tmp1);
  return flag;
}

void fpm_add(fpm_t *c, fpm_t *a, fpm_t *b){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_add(&c->x[i],&a->x[i],&b->x[i]);
  }
}
void fpm_sub(fpm_t *c, fpm_t *a, fpm_t *b){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_sub(&c->x[i],&a->x[i],&b->x[i]);
  }
}
void fpm_mul(fpm_t *c, fpm_t *a, fpm_t *b){
  int i,j,k,l,count;
  int x;
  int h = h_degree_m;
  int m = DEGREE_EXTENTION_FIELD_M;
  fp_t v[DEGREE_EXTENTION_FIELD_M+1];
  fp_t tmp1,tmp2,tmp3;

  //Evaluation step
  //step8

  for(l=0;l<=m-1;l++){
    fp_mul(&v[l],&a->x[l],&b->x[l]);
  }
  //step9
  fp_set_ui(&v[m],0);

  //step10~13
  count = 0;
  for(i=0;i<=m-2;i++){
    for(j=i+1;j<=m-1;j++){
      fp_sub(&tmp1,&a->x[i],&a->x[j]);//tmp1 = a[i]-a[j]
      fp_sub(&tmp2,&b->x[i],&b->x[j]);//tmp2 = b[i]-b[j]
      fp_mul(&tmp3,&tmp1,&tmp2);//tmp3 =  (a[i]-a[j])(b[i]-b[j])
      for(k=0;k<=h-1;k++){
        x = n_ijk_degree_m[count];
        fp_add(&v[x],&v[x],&tmp3);
        count++;
      }
    }
  }
  //step14
  x = h%2;
  if(x==1){
    fp_set_ui(&tmp1,h);
    fp_mul(&tmp2,&tmp1,&v[m]);
    for(l=0;l<=m-1;l++){
      fp_sub(&c->x[l],&tmp2,&v[l]);

    }
  }
  else{
    for(l=0;l<=m-1;l++){
      fp_set_neg(&c->x[l],&v[l]);
    }
  }
}
void fpm_sqr(fpm_t *b, fpm_t *a){
  fpm_mul(b,a,a);
}
void fpm_frobenius(fpm_t *b, fpm_t *a){
  int i;
  fpm_t answer;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M-1;i++){
    fp_set(&answer.x[i+1],&a->x[i]);
  }
  fp_set(&answer.x[0],&a->x[DEGREE_EXTENTION_FIELD_M-1]);
  fpm_set(b,&answer);
}
void fpm_inv(fpm_t *b, fpm_t *a){
  //Itoh-Tsujii inversion algorithm
  int i;
  fpm_t answer;
  fpm_t tmp1;
  fp_t scalar;

  fpm_frobenius(&tmp1,a);
  fpm_set(&answer,&tmp1);
  for(i=2;i<DEGREE_EXTENTION_FIELD_M;i++){
    fpm_frobenius(&tmp1,&tmp1);
    fpm_mul(&answer,&answer,&tmp1);
  }
  fpm_mul(&tmp1,a,&answer);

  fp_inv(&scalar,&tmp1.x[0]);
  fp_set_neg(&scalar,&scalar);
  for(i=0;i<DEGREE_EXTENTION_FIELD_M;i++){
    fp_mul(&answer.x[i],&answer.x[i],&scalar);
  }

  fpm_set(b,&answer);


}
void fpm_pow(fpm_t *b,fpm_t *a,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fpm_t tmp;

    fpm_set(&tmp,a);

    for(i=1;i<length; i++){
        fpm_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            fpm_mul(&tmp,a,&tmp);
        }
    }
    fpm_set(b,&tmp);
}
int fpm_legendre(fpm_t *a){
    fpm_t tmp;
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime_z,DEGREE_EXTENTION_FIELD_M);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fpm_pow(&tmp,a,exp);

    if(fpm_cmp_ui(&tmp,1)==0){
        mpz_clear(exp);
        return 1;
    }
    else{
        mpz_clear(exp);
        return -1;
    }
}
void fpm_sqrt(fpm_t *b,fpm_t *a){
    fpm_t x,y,t,k,n,tmp;
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);

    fpm_set_random(&n,state);
    while(fpm_legendre(&n)!=-1){
        fpm_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,DEGREE_EXTENTION_FIELD_M);//p^m-1
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    fpm_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fpm_pow(&x,a,exp);
    fpm_mul(&tmp,&x,&x);
    fpm_mul(&k,&tmp,a);
    fpm_mul(&x,&x,a);
    while(fpm_cmp_ui(&k,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fpm_pow(&tmp,&k,exp);
        while(fpm_cmp_ui(&tmp,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fpm_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fpm_pow(&t,&y,result);
        fpm_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fpm_mul(&x,&x,&t);
        fpm_mul(&k,&k,&y);
    }
    fpm_set(b,&x);

    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}
