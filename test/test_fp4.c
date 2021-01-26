#include <gmp.h>
#include "ELiPS/fp12cv.h"
#include "ELiPS/matrix.h"
#include "ELiPS/bls12.h"

void test_fp4cv();
void test_fp12cv();
void test_matrix();

int main(void){
  gmp_randinit_default (state);
  gmp_randseed_ui(state,(unsigned long)time(NULL));
  mpn_set_ui(prime,FPLIMB,467);
  mpz_set_ui(prime_z,467);

  mpz_t p,r1,r2,r3;
  mpz_init(p);
  mpz_init(r1);
  mpz_init(r2);
  mpz_init(r3);
  mpz_set_ui(p,1);
  for(int i=1;i<1000;i++){
    mpz_set_ui(p,i);
    if(mpz_probab_prime_p(p,50)>=1){
      mpz_tdiv_r_ui(r1,p,8);
      mpz_tdiv_r_ui(r2,p,5);
      mpz_tdiv_r_ui(r3,p,7);
      if(mpz_cmp_ui(r1,3)==0 && mpz_cmp_ui(r3,5)==0){
        if(mpz_cmp_ui(r2,2)==0 || mpz_cmp_ui(r2,3)==0){
          gmp_printf("%Zd\n",p);
        }
        
      }
      
    }
  }
  
  //bls12_init();
  //test_fp4cv();
  test_fp12cv();
  test_matrix();
   return 0;
}
void test_fp4cv(){
  fp4cv_t tmp1,tmp2,tmp3;

  fp4cv_set_random(&tmp1,state);
  fp4cv_println("tmp1 = ",&tmp1);
  //inv
  fp4cv_inv(&tmp2,&tmp1);
  fp4cv_mul(&tmp3,&tmp1,&tmp2);
  fp4cv_println("tmp1 * tmp1^-1 = ",&tmp3);

}
void test_fp12cv(){
  fp12cv_t tmp1,tmp2,tmp3;
  fp12cv_init(&tmp1);

  fp12cv_set_random(&tmp1,state);
  fp12cv_set_random(&tmp2,state);

  fp12cv_println("tmp1 = ",&tmp1);
  
  fp12cv_frobenius(&tmp3,&tmp1);
  fp12cv_frobenius(&tmp3,&tmp3);
  fp12cv_frobenius(&tmp3,&tmp3);
  fp12cv_frobenius(&tmp3,&tmp3);
  fp12cv_println("tmp1^p4 = ",&tmp3);

  //inv
  fp12cv_inv(&tmp2,&tmp1);
  fp12cv_mul(&tmp3,&tmp1,&tmp2);
  fp12cv_println("tmp1*tmp1^-1 = ",&tmp3);

  fp12cv_frobenius(&tmp3,&tmp1);
  fp12cv_println("tmp1^f = ",&tmp3);
  fp12cv_pow(&tmp2,&tmp1,prime_z);
  fp12cv_println("tmp1^p = ",&tmp2);

  while(1){
    fp12cv_set_random(&tmp1,state);
    if(fp12cv_legendre(&tmp1)==1){
      break;
    }
  }
  fp12cv_println("tmp1 = ",&tmp1);
  fp12cv_sqrt(&tmp2,&tmp1);
  fp12cv_println("tmp1^1/2 = ",&tmp2);

  fp12cv_sqr(&tmp3,&tmp2);
  fp12cv_println("(tmp1^1/2)2 = ",&tmp3);
  //
  fp12cv_set_ui(&tmp2,1);
  fp12cv_println("1 = ",&tmp2);
  printf("%d\n",fp12cv_cmp_one(&tmp2));

}
void test_matrix(){
  matrix_build_fp12_and_fp12cv();
}