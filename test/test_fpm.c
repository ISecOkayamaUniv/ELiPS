#include <gmp.h>
#include "ELiPS/fpm.h"
#include "ELiPS/fpm2.h"
#include "ELiPS/bls.h"

void test_fpm();
void test_fpm2();
void test_pairing();

int main(void){
  gmp_randinit_default (state);
  gmp_randseed_ui(state,(unsigned long)time(NULL));


  bls_build();
  fpm_build();
  fpm2_build();
  test_pairing();
  //test_fpm();
  //test_fpm2();
  return 0;
}
void test_fpm(){
  fpm_t tmp1,tmp2,tmp3;
  fpm_build();
  fpm_set_random(&tmp1,state);
  fpm_set_random(&tmp2,state);
  fpm_println("tmp1 = ",&tmp1);
  fpm_println("tmp2 = ",&tmp2);
  //add
  fpm_add(&tmp3,&tmp1,&tmp2);
  fpm_println("tmp1 + tmp2 = ",&tmp3);
  //sub
  fpm_sub(&tmp3,&tmp1,&tmp2);
  fpm_println("tmp1 - tmp2 = ",&tmp3);
  //mul
  fpm_mul(&tmp3,&tmp1,&tmp2);
  fpm_println("tmp1 * tmp2 = ",&tmp3);
  /*
  //set
  fpm_set_ui(&tmp1,1);
  fpm_println("tmp1 = ",&tmp1);
  fpm_set_neg(&tmp1,&tmp2);
  fpm_add(&tmp3,&tmp1,&tmp2);
  fpm_println("tmp1 + tmp2 = ",&tmp3);
  //frobenius
  fpm_set_random(&tmp1,state);
  fpm_println("tmp1 = ",&tmp1);
  fpm_frobenius(&tmp3,&tmp1);
  fpm_println("tmp1^f = ",&tmp3);

  fpm_mul(&tmp3,&tmp1,&tmp1);
  for(int i=2;i<17;i++){
    fpm_mul(&tmp3,&tmp3,&tmp1);
  }
  fpm_println("tmp1^p = ",&tmp3);

  fpm_frobenius(&tmp3,&tmp1);
  fpm_frobenius(&tmp3,&tmp3);
  fpm_println("tmp1^2f = ",&tmp3);

  fpm_mul(&tmp3,&tmp1,&tmp1);
  for(int i=2;i<289;i++){
    fpm_mul(&tmp3,&tmp3,&tmp1);
  }
  fpm_println("tmp1^p2 = ",&tmp3);

  mpz_t exp;
  mpz_init(exp);
  mpz_set_ui(exp,289);
  fpm_pow(&tmp3,&tmp1,exp);
  fpm_println("tmp1^289 = ",&tmp3);
  //inv
  fpm_inv(&tmp2,&tmp1);
  fpm_mul(&tmp3,&tmp1,&tmp2);
  fpm_println("tmp1 * tmp1^-1 = ",&tmp3);
  //sqrt
  while(1){
    fpm_set_random(&tmp1,state);
    if(fpm_legendre(&tmp1)==1){
      break;
    }
  }
  fpm_println("tmp1 = ",&tmp1);
  fpm_sqrt(&tmp2,&tmp1);
  fpm_sqr(&tmp3,&tmp2);
  fpm_println("tmp3 = ",&tmp3);
  */
}
void test_fpm2(){
  fpm2_t tmp1,tmp2,tmp3;
  fpm2_build();
  fpm2_set_random(&tmp1,state);
  fpm2_set_random(&tmp2,state);
  fpm2_println("tmp1 = ",&tmp1);
  fpm2_println("tmp2 = ",&tmp2);
  //add
  fpm2_add(&tmp3,&tmp1,&tmp2);
  fpm2_println("tmp1 + tmp2 = ",&tmp3);
  //sub
  fpm2_sub(&tmp3,&tmp1,&tmp2);
  fpm2_println("tmp1 - tmp2 = ",&tmp3);
  //set
  fpm2_set_ui(&tmp1,1);
  fpm2_println("tmp1 = ",&tmp1);
  fpm2_set_neg(&tmp1,&tmp2);
  fpm2_add(&tmp3,&tmp1,&tmp2);
  fpm2_println("tmp1 + tmp2 = ",&tmp3);
  //frobenius
  fpm2_set_random(&tmp1,state);
  fpm2_println("tmp1 = ",&tmp1);
  fpm2_frobenius(&tmp3,&tmp1);
  fpm2_println("tmp1^f = ",&tmp3);

  fpm2_mul(&tmp3,&tmp1,&tmp1);
  for(int i=2;i<17;i++){
    fpm2_mul(&tmp3,&tmp3,&tmp1);
  }
  fpm2_println("tmp1^p = ",&tmp3);

  fpm2_frobenius(&tmp3,&tmp1);
  fpm2_frobenius(&tmp3,&tmp3);
  fpm2_println("tmp1^2f = ",&tmp3);

  fpm2_mul(&tmp3,&tmp1,&tmp1);
  for(int i=2;i<289;i++){
    fpm2_mul(&tmp3,&tmp3,&tmp1);
  }
  fpm2_println("tmp1^p2 = ",&tmp3);

  mpz_t exp;
  mpz_init(exp);
  mpz_set_ui(exp,289);
  fpm2_pow(&tmp3,&tmp1,exp);
  fpm2_println("tmp1^289 = ",&tmp3);
  //inv
  fpm2_inv(&tmp2,&tmp1);
  fpm2_mul(&tmp3,&tmp1,&tmp2);
  fpm2_println("tmp1 * tmp1^-1 = ",&tmp3);
  //check
  fpm2_set_random(&tmp1,state);
  fpm2_frobenius(&tmp2,&tmp1);
  fpm2_frobenius(&tmp2,&tmp2);
  fpm2_frobenius(&tmp2,&tmp2);
  fpm2_frobenius(&tmp2,&tmp2);
  fpm2_mul(&tmp3,&tmp1,&tmp2);
  fpm2_frobenius(&tmp2,&tmp2);
  fpm2_frobenius(&tmp2,&tmp2);
  fpm2_frobenius(&tmp2,&tmp2);
  fpm2_frobenius(&tmp2,&tmp2);
  fpm2_mul(&tmp3,&tmp3,&tmp2);
  fpm2_println("tmp1 * tmp1^p4 * tmp1^p8 = ",&tmp3);
  //sqrt
  while(1){
    fpm2_set_random(&tmp1,state);
    if(fpm2_legendre(&tmp1)==1){
      break;
    }
  }
  fpm2_println("tmp1 = ",&tmp1);
  fpm2_sqrt(&tmp2,&tmp1);
  fpm2_sqr(&tmp3,&tmp2);
  fpm2_println("tmp3 = ",&tmp3);
}
void test_pairing(){
  fpm2_t a;
  efpm2_t g1,g2;
  efpm2_t p1,p2;

  efpm2_init(&g1);
  efpm2_init(&g2);
  efpm2_init(&p1);
  efpm2_init(&p2);

  bls_generate_g1_for_efpm2(&g1);
  efpm2_printf("g1 = ",&g1);
  bls_generate_g2_for_efpm2(&g2);
  efpm2_printf("g2 = ",&g2);

  bls_miller_for_optate_basic(&a,&g1,&g2);



 /*
  efpm2_init(&p2);
  efpm2_init(&g2);
  bls_make_g2(&g2);
  efpm2_printf("g2 = ",&g2);
  efpm2_scm(&p2,&g2,order_z);
  efpm2_printf("[r]g2 = ",&p2);
  efpm2_scm(&p2,&g2,efp_total);
  efpm2_printf("[E]g2 = ",&p2);
  */
}