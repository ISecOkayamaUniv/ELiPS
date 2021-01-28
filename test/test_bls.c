#include <gmp.h>
#include "ELiPS/bls.h"
#include "ELiPS/bls12.h"
#include <ELiPS/fpm2.h>
#include "ELiPS/time.h"



void test_pairing();
void test_twist();
void test_X();

int main(void){
  gmp_randinit_default (state);
  gmp_randseed_ui(state,(unsigned long)time(NULL));

  
  bls_build();
  fpm_build();
  fpm2_build();
  fp2cv_build();
  set_twist_g();
  test_twist();
  test_pairing();
  

  return 0;
}

void test_pairing(){
  mpz_t a,b;
  efpm2_t g1,g2;
  efpm2_t p1,p2,q1,q2;
  fpm2_t f1,f2;
  struct timeval tv_A,tv_B;
  float miller_time=0,final_time=0,pairing_time=0;

  mpz_init(a);
  mpz_init(b);
  efpm2_init(&g1);
  efpm2_init(&g2);
  efpm2_init(&p1);
  efpm2_init(&p2);
  efpm2_init(&q1);
  efpm2_init(&q2);
  //set
  mpz_set_ui(a,123);
  mpz_set_ui(b,45);
  bls_generate_g1_for_efpm2(&g1);
  printf("g1\n");
  bls_generate_g2_for_efpm2(&g2);
  printf("g2\n");
  efpm2_scm(&p1,&g1,a);//a[P]
  efpm2_scm(&q1,&g2,b);//b[Q]
  efpm2_scm(&p2,&g1,b);//b[P]
  efpm2_scm(&q2,&g2,a);//a[Q]

  
  bls_optate_pairing_basic(&f1,&p1,&q1);
  bls_optate_pairing_basic(&f2,&p2,&q2);
  

  if(fpm2_cmp(&f1,&f2)==0){
    printf("success optimal pairing\n");
  }
  else{

    printf("failed optimal pairing\n");
  }
  //fpm2_println("f1 = \n",&f1);
  //fpm2_println("f2 = \n",&f2);

  
  gettimeofday(&tv_A,NULL);
  //bls_miller_for_optate_basic(&f1,&p1,&q1);
  bls_2i3_miller_for_optate(&f1,&p1,&q1);
  gettimeofday(&tv_B,NULL);
  miller_time=timedifference_msec(tv_A,tv_B);
  printf("miller time = %f\n",miller_time);

  gettimeofday(&tv_A,NULL);
  bls_final_exp(&f1,&f1);
  gettimeofday(&tv_B,NULL);
  final_time=timedifference_msec(tv_A,tv_B);
  printf("final time = %f\n",final_time);
  
}

void test_twist(){
  mpz_t a;
  efpm2_t g1,g2;
  efpm2_t p1,p2,q1,q2;
  efp2cv_t cv1;

  mpz_init(a);
  efpm2_init(&g1);
  efpm2_init(&g2);
  efpm2_init(&p1);
  efpm2_init(&p2);
  efpm2_init(&q1);
  efpm2_init(&q2);
  efp2cv_init(&cv1);
  //set
  mpz_set_ui(a,12344);
  
  bls_generate_g2_for_efpm2(&g2);
  efpm2_scm(&q1,&g2,a);//a[Q]
  efpm2_println("[a]g2 = \n",&q1);
  efpm2_to_efp2cv_for_g2(&cv1,&q1);
  efp2cv_to_efpm2_for_g2(&q2,&cv1);
  efpm2_println("twisted [a]g2 = \n",&q2);
  fpm2_t tmp1,tmp2;
  fpm2_mul(&tmp1,&g2.x,&twist_g_cbrt);
  fpm2_mul(&tmp2,&g2.y,&twist_g_sqrt);
  fpm2_println("x = \n",&tmp1);
  fpm2_println("y = \n",&tmp2);

}
void test_X(){
  fpm2_t tmp1,tmp2,tmp3;
  fpm2_set_random(&tmp1,state);
  fpm2_pow(&tmp2,&tmp1,X_z);
  bls_fpm2_pow_X(&tmp3,&tmp1);
  printf("cmp = %d\n",fpm2_cmp(&tmp2,&tmp3));
}
