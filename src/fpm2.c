#include <ELiPS/fpm2.h>


void fpm2_make_cvma(){
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
    r = h*DEGREE_EXTENTION_FIELD_M2 + 1;
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
            mpz_set_ui(tmp1,h*DEGREE_EXTENTION_FIELD_M2/e);
            mpz_gcd_ui(tmp1,tmp1,DEGREE_EXTENTION_FIELD_M2);
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
  r_degree_m2 = r;
  h_degree_m2 = h;
  mpz_clear(tmp1);
}
void fpm2_preparate_mul(){
  int i,j,k,l;
  int x1,x2,x3,x4;//代入変数
  int flag = 0;
  int d;
  int r = r_degree_m2;
  int h = h_degree_m2;
  int m = DEGREE_EXTENTION_FIELD_M2;
  int e[SIZE_ALL];//サイズの修正
  int count;
  mpz_t tmp1;
  //init
  for(i=0;i<SIZE_ALL;i++){
    e[i] = 0;
    n_ijk_degree_m2[i] = 0;
  }
  //init
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
    //printf("p^i = %d (mod r)\n",x);
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
        //printf("p^j = %d (mod r)\n",y);
        x4 = x1 + x4;
        //printf("p^i + p^j = %d\n",y);
        x4 = x4%r;
        //
        //printf("x = %d\n",x);
        n_ijk_degree_m2[count] = e[x4];
        count++;
      }
    }
  }

  //clear
  mpz_clear(tmp1);
}
void fpm2_print_paramerter(){
  printf("\n-----  %d->%d use cvma type<h,m>  -----\n",DEGREE_EXTENTION_FIELD_M,DEGREE_EXTENTION_FIELD_M*DEGREE_EXTENTION_FIELD_M2);
  printf("h = %d r = %d\n",h_degree_m2,r_degree_m2);
}
void fpm2_build(){
  fpm2_make_cvma();
  fpm2_preparate_mul();
  fpm2_print_paramerter();
}
void fpm2_init(fpm2_t *a){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_init(&a->x[i]);
  }
}
void fpm2_printf(char *str,fpm2_t *a){
  int i;
  gmp_printf("%s(",str);
  fpm_printf("",&a->x[0]);
  for(i=1;i<DEGREE_EXTENTION_FIELD_M2;i++){
    gmp_printf(",");
    fpm_printf("",&a->x[i]);
  }
  gmp_printf(")");
}
void fpm2_println(char *str,fpm2_t *a){
  int i;
  gmp_printf("%s",str);
  
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_println("",&a->x[i]);
  }
}
void fpm2_set_random(fpm2_t *a,gmp_randstate_t state){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_set_random(&a->x[i],state);
  }
}
void fpm2_set(fpm2_t *b, fpm2_t *a){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_set(&b->x[i],&a->x[i]);
  }
}
void fpm2_set_ui(fpm2_t *a, unsigned long int b){
  int i;
  fpm2_t tmp1;

  if(b==0){
    fpm2_init(&tmp1);
    fpm2_set(a,&tmp1);
  }
  else{
    for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
      fpm_set_ui(&tmp1.x[i],b);
      fpm_set_neg(&a->x[i],&tmp1.x[i]);
    }
  }
}
void fpm2_set_fp(fpm2_t *a, fp_t *b){
  int i;
  fpm2_t answer;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_set_fp(&answer.x[i],b);
  }
  fpm2_set_neg(a,&answer);
}
void fpm2_set_mpn(fpm2_t *a, mp_limb_t *b){
  fp_t tmp1;
  fp_set_mpn(&tmp1,b);
  fpm2_set_fp(a,&tmp1);
}
void fpm2_set_neg(fpm2_t *b, fpm2_t *a){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_set_neg(&b->x[i],&a->x[i]);
  }
}
int fpm2_cmp(fpm2_t *a, fpm2_t *b){
  int flag = 0;
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    if(fpm_cmp(&a->x[i],&b->x[i])!=0){
      flag = 1;
    }
  }
  return flag;
}
int fpm2_cmp_ui(fpm2_t *a, unsigned long int b){
  int flag = 0;
  fpm2_t tmp1;

  fpm2_set_ui(&tmp1,b);
  flag = fpm2_cmp(a,&tmp1);
  return flag;
}
void fpm2_add(fpm2_t *c, fpm2_t *a, fpm2_t *b){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_add(&c->x[i],&a->x[i],&b->x[i]);
  }
}
void fpm2_sub(fpm2_t *c, fpm2_t *a, fpm2_t *b){
  int i;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_sub(&c->x[i],&a->x[i],&b->x[i]);
  }
}
void fpm2_mul(fpm2_t *c, fpm2_t *a, fpm2_t *b){
  int i,j,k,l,count;
  int x;
  int h = h_degree_m2;
  int m = DEGREE_EXTENTION_FIELD_M2;
  fpm_t v[DEGREE_EXTENTION_FIELD_M2+1];
  fpm_t tmp1,tmp2,tmp3;
  fpm2_t answer;

  //Evaluation step
  //step8

  for(l=0;l<=m-1;l++){
    fpm_mul(&v[l],&a->x[l],&b->x[l]);
  }
  //step9
  fpm_set_ui(&v[m],0);

  //step10~13
  count = 0;
  for(i=0;i<=m-2;i++){
    for(j=i+1;j<=m-1;j++){
      fpm_sub(&tmp1,&a->x[i],&a->x[j]);//tmp1 = a[i]-a[j]
      fpm_sub(&tmp2,&b->x[i],&b->x[j]);//tmp2 = b[i]-b[j]
      fpm_mul(&tmp3,&tmp1,&tmp2);//tmp3 =  (a[i]-a[j])(b[i]-b[j])
      for(k=0;k<=h-1;k++){
        x = n_ijk_degree_m2[count];
        fpm_add(&v[x],&v[x],&tmp3);
        count++;
      }
    }
  }
  //step14
  x = h%2;
  if(x==1){
    fpm_set_ui(&tmp1,h);
    fpm_mul(&tmp2,&tmp1,&v[m]);
    for(l=0;l<=m-1;l++){
      fpm_sub(&answer.x[l],&tmp2,&v[l]);
    }
  }
  else{
    for(l=0;l<=m-1;l++){
      fpm_set_neg(&answer.x[l],&v[l]);
    }
  }
  //mod
  fpm2_set(c,&answer);
}
void fpm2_sqr(fpm2_t *b, fpm2_t *a){
  fpm2_mul(b,a,a);
}
void fpm2_frobenius(fpm2_t *b, fpm2_t *a){
  int i;
  fpm2_t answer;
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2-1;i++){
    fpm_frobenius(&answer.x[i+1],&a->x[i]);
  }
  fpm_frobenius(&answer.x[0],&a->x[DEGREE_EXTENTION_FIELD_M2-1]);
  fpm2_set(b,&answer);
}
void fpm2_inv(fpm2_t *b, fpm2_t *a){
  //Itoh-Tsujii inversion algorithm
  int i,j;
  fpm2_t answer;
  fpm2_t tmp1,tmp2;
  fpm_t scalar;

  fpm2_set(&tmp1,a);
  fpm2_set_ui(&tmp2,1);
  fpm2_set(&answer,a);
  for(i=1;i<DEGREE_EXTENTION_FIELD_M2;i++){
    for(j=0;j<DEGREE_EXTENTION_FIELD_M;j++){
      fpm2_frobenius(&tmp1,&tmp1);
    }
    fpm2_mul(&tmp2,&tmp2,&tmp1);
  }
  fpm2_mul(&tmp1,a,&tmp2);
  fpm_inv(&scalar,&tmp1.x[0]);
  for(i=0;i<DEGREE_EXTENTION_FIELD_M2;i++){
    fpm_set(&tmp1.x[i],&scalar);
  }
  fpm2_mul(b,&tmp1,&tmp2);
}
void fpm2_pow(fpm2_t *b,fpm2_t *a,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fpm2_t tmp;

    fpm2_set(&tmp,a);

    for(i=1;i<length; i++){
        fpm2_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            fpm2_mul(&tmp,a,&tmp);
        }
    }
    fpm2_set(b,&tmp);
}
int fpm2_legendre(fpm2_t *a){
    fpm2_t tmp;
    mpz_t exp;
    mpz_init(exp);
    //exp = (p^m*m2-1)/2
    mpz_pow_ui(exp,prime_z,DEGREE_EXTENTION_FIELD_M*DEGREE_EXTENTION_FIELD_M2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fpm2_pow(&tmp,a,exp);

    if(fpm2_cmp_ui(&tmp,1)==0){
        mpz_clear(exp);
        return 1;
    }
    else{
        mpz_clear(exp);
        return -1;
    }
}
void fpm2_sqrt(fpm2_t *b,fpm2_t *a){
    fpm2_t x,y,t,k,n,tmp;
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);

    fpm2_set_random(&n,state);
    while(fpm2_legendre(&n)!=-1){
        fpm2_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,DEGREE_EXTENTION_FIELD_M*DEGREE_EXTENTION_FIELD_M2);//p^m*m2-1
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    fpm2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fpm2_pow(&x,a,exp);
    fpm2_mul(&tmp,&x,&x);
    fpm2_mul(&k,&tmp,a);
    fpm2_mul(&x,&x,a);
    while(fpm2_cmp_ui(&k,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fpm2_pow(&tmp,&k,exp);
        while(fpm2_cmp_ui(&tmp,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fpm2_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fpm2_pow(&t,&y,result);
        fpm2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fpm2_mul(&x,&x,&t);
        fpm2_mul(&k,&k,&y);
    }
    fpm2_set(b,&x);

    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}
