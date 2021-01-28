#include <ELiPS/bls_generate.h>

void bls_init(){
    mpz_init(X_z);
    mpz_init(prime_z);
    mpz_init(order_z);
    mpz_init(trace_z);
    mpn_zero(curve_b,FPLIMB);

}

void bls_clear(){
    mpz_clear(X_z);
    mpz_clear(prime_z);
    mpz_clear(order_z);
    mpz_clear(trace_z);
}

void bls_build(){
  int reps = 50;

  bls_init();

  bls_generate_X();
  bls_generate_parameter_of_bls12();

  //bls6
  //mpz_setbit(X_z,115);
  //mpz_setbit(X_z,21);
  //bls_search_parameter_all();


  if(mpz_probab_prime_p(prime_z,reps) && mpz_probab_prime_p(order_z,reps)){
    printf("Successful generate parameter\n");
    bls_weil();
    //bls_weil_of_bls12();
    bls_search_curve();
    bls_print_parameter();
  }
  else{
    printf("Error\n");
    printf("you cannot use this paramerter\n");
    bls_clear();
  }
}

void bls_print_parameter(){
  gmp_printf("\n--------------------BLS--------------------\n");
  gmp_printf("parameter\n");
  gmp_printf("X = %Zd\n",X_z);
  gmp_printf("prime = %Zd\n",prime_z);
  gmp_printf("trace = %Zd\n",trace_z);
  gmp_printf("order = %Zd\n",order_z);
  gmp_printf("\nelliptic curve\n");
  gmp_printf("BLS12 : y^2 = x^3 + %Nu\n", curve_b,FPLIMB);
}

void bls_generate_X(){
  int i;
  mpz_t tmp1;
  mpz_init(tmp1);

  //set bit

  bls_X_length = 77;
  bls_X_binary[77] = 1;
  bls_X_binary[11] = 1;
  bls_X_binary[9] = -1;
  bls_X_binary[6] = -1;

  bls_X2_binary[76] = 1;
  bls_X2_binary[10] = 1;
  bls_X2_binary[8] = -1;
  bls_X2_binary[5] =- 1;


/*
  bls_X_length = 76;
  bls_X_binary[76] = 1;
  bls_X_binary[59] = -1;
  bls_X_binary[11] = -1;
  bls_X_binary[5] = -1;
*/

  //set X_z
  mpz_set_ui(X_z,0);
  for(i=0;i<SIZE_BLS;i++){
    if(bls_X_binary[i]==1){
      mpz_setbit(tmp1,i);
      mpz_add(X_z,X_z,tmp1);
      mpz_clrbit(tmp1,i);
    }
    else if(bls_X_binary[i]==-1){
      mpz_setbit(tmp1,i);
      mpz_sub(X_z,X_z,tmp1);
      mpz_clrbit(tmp1,i);
    }
  }
  mpz_clear(tmp1);
}

void bls_generate_parameter_of_bls12(){
  mpz_t tmp1,tmp2;
  mpz_init(tmp1);
  mpz_init(tmp2);

  //t = X + 1
  mpz_add_ui(trace_z,X_z,1);
  //r = X^4 - X^2 + 1
  mpz_pow_ui(tmp1,X_z,4);
  mpz_pow_ui(tmp2,X_z,2);
  mpz_sub(order_z,tmp1,tmp2);
  mpz_add_ui(order_z,order_z,1);
  //p = (X-1)^2(X^4-X^2+1)/3 + X
  mpz_sub_ui(tmp1,X_z,1);
  mpz_pow_ui(tmp1,tmp1,2);
  mpz_mul(tmp2,tmp1,order_z);
  mpz_tdiv_qr_ui(tmp2,tmp1,tmp2,3);
  mpz_add(prime_z,tmp2,X_z);

  mpn_set_mpz(prime,prime_z);

  mpz_clear(tmp1);
  mpz_clear(tmp2);
}

void bls_search_parameter_all(){
  int length = 0;
  int i,j,r,reps = 50;
  mpz_t psy[DEGREE_MN + 1];
  mpz_t z;
  mpz_t tmp1,tmp2;
  //init
  mpz_init(z);
  mpz_init(tmp1);
  mpz_init(tmp2);
  for(i=0;i<=DEGREE_MN;i++){
    mpz_init(psy[i]);
  }

  while(1){
    //z = (X-1)/6
    mpz_sub_ui(tmp1,X_z,1);
    mpz_tdiv_qr_ui(z,tmp2,tmp1,6);
    if(mpz_cmp_ui(tmp2,0)==0){
      //psy
      mpz_sub_ui(psy[1],X_z,1);//Φ1 = X - 1
      mpz_add_ui(psy[2],X_z,1);//Φ2 = X + 1
      for(i=3;i<=DEGREE_MN;i++){
        mpz_pow_ui(tmp1,X_z,i);
        mpz_sub_ui(psy[i],tmp1,1);
        for(j=1;j<i;j++){
          r = i%j;
          if(r==0){
            mpz_tdiv_q(psy[i],psy[i],psy[j]);
          }
        }
      }

      mpz_add_ui(trace_z,X_z,1);
      mpz_set(order_z,psy[DEGREE_MN]);
      i = DEGREE_MN/6;
      mpz_pow_ui(tmp1,X_z,i);
      mpz_mul_ui(tmp1,tmp1,2);
      mpz_sub_ui(tmp1,tmp1,1);
      mpz_mul_ui(tmp2,z,2);
      mpz_mul(tmp1,tmp1,tmp2);
      mpz_pow_ui(tmp1,tmp1,2);
      mpz_mul_ui(tmp1,tmp1,3);
      mpz_mul_ui(tmp2,z,6);
      mpz_pow_ui(tmp2,tmp2,2);
      mpz_add(tmp1,tmp1,tmp2);
      mpz_tdiv_q_ui(tmp1,tmp1,4);
      mpz_add(prime_z,tmp1,X_z);
      if(mpz_probab_prime_p(prime_z,reps) && mpz_probab_prime_p(order_z,reps)){
        break;
      }
    }
    mpz_add_ui(X_z,X_z,1);
  }

  length = (int)mpz_sizeinbase(prime_z,2);
  printf("log2(prime) = %d\n",length);
  bls_X_length = (int)mpz_sizeinbase(X_z,2);
  bls_X_length = bls_X_length - 1;
  printf("log2(X) = %d\n",bls_X_length);
  //set
  mpz_t s;
  mpz_init(s);
  mpz_set(s,X_z);
  for(i=0;mpz_popcount(s);i++){
    if(mpz_scan1(s,0)==i){
      bls_X_binary[i] = 1;
      mpz_clrbit(s,i);
    }
    printf(" %dbit  %d\n",i,bls_X_binary[i]);
  }
  mpz_clear(s);
  mpn_set_mpz(prime,prime_z);
  mpn_mul_ui(prime_carry,prime,FPLIMB,2);

  //clear
  mpz_clear(z);
  mpz_clear(tmp1);
  mpz_clear(tmp2);
}

void bls_weil_of_bls12(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);

    //efp_total
    mpz_add_ui(buf,prime_z,1);
    mpz_sub(efp_total,buf,trace_z);

    //t2←α^2+β^2
    mpz_pow_ui(t2,trace_z,2);
    mpz_mul_ui(buf,prime_z,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,prime_z,2);

    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);

    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);

    //efp12_232_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(efpm2_total,buf,1);

    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}
void bls_weil_of_bls48(){
  mpz_t t2,t6,t12,t24,t48,p2,p6,p12,p24,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(t24);
    mpz_init(t48);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(p12);
    mpz_init(p24);
    mpz_init(buf);

    //efp_total
    mpz_add_ui(buf,prime_z,1);
    mpz_sub(efp_total,buf,trace_z);

    //t2←α^2+β^2
    mpz_pow_ui(t2,trace_z,2);
    mpz_mul_ui(buf,prime_z,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,prime_z,2);

    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);

    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    mpz_pow_ui(p12,p6,2);

    //α^24+β^24
    mpz_pow_ui(t24,t12,2);
    mpz_mul_ui(buf,p12,2);
    mpz_sub(t24,t24,buf);
    mpz_pow_ui(p24,p12,2);

    //α^48+β^48
    mpz_pow_ui(t48,t24,2);
    mpz_mul_ui(buf,p24,2);
    mpz_sub(t48,t48,buf);

    //efpm2_total
    mpz_pow_ui(buf,p24,2);
    mpz_sub(buf,buf,t48);
    mpz_add_ui(efpm2_total,buf,1);

    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(t24);
    mpz_clear(t48);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(p12);
    mpz_clear(p24);
    mpz_clear(buf);
}

void bls_weil(){
  int i,j;
  int r,count,cmb;
  mpz_t t[DEGREE_MN + 1];
  mpz_t tmp1,tmp2;
  //init
  for(i=0;i<=DEGREE_MN;i++){
    mpz_init(t[i]);
  }
  mpz_init(tmp1);
  mpz_init(tmp2);
  //t0
  mpz_set_ui(t[0],1);
  //t1
  mpz_set(t[1],trace_z);
  //t2
  mpz_pow_ui(tmp1,t[1],2);
  mpz_mul_ui(tmp2,prime_z,2);
  mpz_sub(t[2],tmp1,tmp2);

  //binomial theorem
  for(i=3;i<=DEGREE_MN;i++){
    r = i%2;
    if(r==0){
      mpz_pow_ui(t[i],t[1],i);
      count = 1;
      for(j=i-2;j>=0;j=j-2){
        mpz_pow_ui(tmp1,prime_z,count);
        cmb = combination(i,count);
        mpz_mul_ui(tmp1,tmp1,cmb);
        mpz_mul(tmp1,tmp1,t[j]);
        mpz_sub(t[i],t[i],tmp1);
        count++;
      }
    }
    else{
      mpz_pow_ui(t[i],t[1],i);
      count = 1;
      for(j=i-2;j>=1;j=j-2){
        mpz_pow_ui(tmp1,prime_z,count);
        cmb = combination(i,count);
        mpz_mul_ui(tmp1,tmp1,cmb);
        mpz_mul(tmp1,tmp1,t[j]);
        mpz_sub(t[i],t[i],tmp1);
        count++;
      }
    }
  }
  //
  //efp_total
  mpz_add_ui(tmp1,prime_z,1);
  mpz_sub(efp_total,tmp1,trace_z);
  //efpm2_total
  mpz_pow_ui(tmp1,prime_z,DEGREE_MN);
  mpz_add_ui(tmp2,tmp1,1);
  mpz_sub(efpm2_total,tmp2,t[DEGREE_MN]);

  //clear
  for(i=0;i<=DEGREE_MN;i++){
    mpz_clear(t[i]);
  }
}

int combination(int n, int r){
    if (r == 0 || r == n)
        return 1;
    else if (r == 1)
        return n;
    return (combination(n - 1, r - 1) + combination(n - 1, r));
}

void bls_search_curve(){
  efp_t p;
  efp_init(&p);

  while(1){
    mpn_add_ui(curve_b,curve_b,FPLIMB,1);
    efp_set_random(&p);
    efp_scm(&p,&p,efp_total);
    if(p.infinity == 1){
      break;
    }
  }
}
