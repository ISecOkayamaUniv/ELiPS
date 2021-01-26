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

  if(mpz_probab_prime_p(prime_z,reps) && mpz_probab_prime_p(order_z,reps)){
    printf("Successful generate parameter\n");
    bls_weil_of_bls12();
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
  bls_X_binary[77] = 1;
  bls_X_binary[11] = 1;
  bls_X_binary[9] = -1;
  bls_X_binary[6] = -1;

/*
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
  //
  mpn_set_mpz(prime,prime_z);

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

void bls_search_curve(){
  efp_t p;
  efp_init(&p);

  while(1){
    mpn_add_ui(curve_b,curve_b,FPLIMB,1);
    efp_rational_point(&p);
    efp_scm(&p,&p,efp_total);
    if(p.infinity == 1){
      break;
    }
  }
}