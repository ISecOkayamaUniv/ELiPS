#include <gmp.h>
#include "ELiPS/bls12.h"
#include <ELiPS/matrix.h>
#include <ELiPS/fp12cv.h>

void test_13();
void test_operation();
void test_matrix_cv();
void test_matrix_fpm();
void test_matrix_fpm2();
void test_matrix_fpm_fp12cv();

int main(void){
    bls12_init();
    bls12_print_parameters();
    fpm_build();
    
    //fpm2_build();
    //test_matrix_fpm2();
    //test_matrix_fpm();
    //test_matrix_fpm_fp12cv();
    test_matrix_cv();

    return 0;
}
void test_13(){
    int i;
    mpz_t exp,r,q;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(r);
    mpz_set(exp,prime_z);
    for(i=1;i<13;i++){
        mpz_pow_ui(exp,prime_z,i);
        mpz_sub_ui(exp,exp,1);
        mpz_tdiv_qr_ui(q,r,exp,13);//exp = (p^12-1)/13
        gmp_printf("i = %d r = %Zd\n",i,r);
    }
  
}
void test_matrix_cv(){
    printf("\n------test fp12,fpm,fp12cv------\n");
    fp12_t a1,a2,a3,a4;
    fpm_t b1,b2,b3;
    fp12cv_t c1,c2,c3;

    matrix_build();
    //first
    fp12_set_random(&a1,state);
    fp12_to_fpm(&b1,&a1);
    fpm_to_fp12cv(&c1,&b1);
    fp12cv_to_fpm(&b3,&c1);
    fpm_to_fp12(&a3,&b3);
    printf("first cmp = %d\n",fp12_cmp(&a1,&a3));
    //first 2
    fp12_to_fp12cv(&c1,&a1);
    fp12cv_to_fp12(&a3,&c1);
    printf("first2 cmp = %d\n",fp12_cmp(&a1,&a3));

    //2 pattern
    fp12_set_random(&a1,state);
    fp12_to_fpm(&b1,&a1);
    fpm_to_fp12cv(&c1,&b1);
    fp12_to_fp12cv(&c2,&a1);
    printf("2 cmp = %d\n",fp12cv_cmp(&c1,&c2));

    //mul
    fp12_set_random(&a1,state);
    fp12_set_random(&a2,state);
    fp12_mul(&a4,&a1,&a2);
    fp12_to_fpm(&b1,&a1);
    fp12_to_fpm(&b2,&a2);
    fpm_to_fp12cv(&c1,&b1);
    fpm_to_fp12cv(&c2,&b2);
    fp12cv_mul(&c3,&c1,&c2);
    fp12cv_to_fpm(&b3,&c3);
    fpm_to_fp12(&a3,&b3);
    printf("mul cmp = %d\n",fp12_cmp(&a3,&a4));

}
void test_matrix_fpm2(){
    fp12_t a1,a2,a3,a4;
    fpm2_t b1,b2,b3;
    //init
    fp12_init(&a1);
    fp12_init(&a2);
    fp12_init(&a3);
    fp12_init(&a4);
    fpm2_init(&b1);
    fpm2_init(&b2);
    fpm2_init(&b3);
    //build matrix
    matrix_build_fp12_and_fpm2();
    //set random
    fp12_set_random(&a1,state);
    fp12_set_random(&a2,state);
    //first
    fp12_to_fpm2(&b1,&a1);
    fpm2_to_fp12(&a2,&b1);
    printf("first cmp = %d\n",fp12_cmp(&a1,&a2));
    //add
    fp12_add(&a3,&a1,&a2);
    fp12_to_fpm2(&b1,&a1);
    fp12_to_fpm2(&b2,&a2);
    fpm2_add(&b3,&b1,&b2);
    fpm2_to_fp12(&a4,&b3);
    printf("add cmp = %d\n",fp12_cmp(&a3,&a4));
    //mul
    fp12_mul(&a3,&a1,&a2);
    fp12_to_fpm2(&b1,&a1);
    fp12_to_fpm2(&b2,&a2);
    fpm2_mul(&b3,&b1,&b2);
    fpm2_to_fp12(&a4,&b3);
    printf("mul cmp = %d\n",fp12_cmp(&a3,&a4));
    //frobenius
    fp12_frobenius_map_p1(&a3,&a1);
    fp12_to_fpm2(&b1,&a1);
    fpm2_frobenius(&b2,&b1);
    fpm2_to_fp12(&a4,&b2);
    printf("frobenius cmp = %d\n",fp12_cmp(&a3,&a4));

}
void test_matrix_fpm(){
    printf("\n------test fp12 to fpm(12)------\n");
    fp12_t a1,a2,a3,a4;
    fpm_t b1,b2,b3;
    //init
    fp12_init(&a1);
    fp12_init(&a2);
    fp12_init(&a3);
    fp12_init(&a4);
    fpm_init(&b1);
    fpm_init(&b2);
    fpm_init(&b3);
    //build matrix
    matrix_build_fp12_and_fpm();
    //set random
    fp12_set_random(&a1,state);
    fp12_set_random(&a2,state);
    //check
    fp12_to_fpm(&b1,&a1);
    fpm_to_fp12(&a2,&b1);
    printf("first cmp = %d\n",fp12_cmp(&a1,&a2));
    //add
    fp12_add(&a3,&a1,&a2);
    fp12_to_fpm(&b1,&a1);
    fp12_to_fpm(&b2,&a2);
    fpm_add(&b3,&b1,&b2);
    fpm_to_fp12(&a4,&b3);
    printf("add cmp = %d\n",fp12_cmp(&a3,&a4));
    //mul
    fp12_mul(&a3,&a1,&a2);
    fp12_to_fpm(&b1,&a1);
    fp12_to_fpm(&b2,&a2);
    fpm_mul(&b3,&b1,&b2);
    fpm_to_fp12(&a4,&b3);
    printf("mul cmp = %d\n",fp12_cmp(&a3,&a4));
    //frobenius
    fp12_frobenius_map_p1(&a3,&a1);
    fp12_to_fpm(&b1,&a1);
    fpm_frobenius(&b2,&b1);
    fpm_to_fp12(&a4,&b2);
    printf("frobenius cmp = %d\n",fp12_cmp(&a3,&a4));

}
void test_matrix_fpm_fp12cv(){
    printf("\n------test fp12cv to fpm(12)------\n");
    int mul_cmp=100;
    fp12cv_t a1,a2,a3,a4;
    fpm_t b1,b2,b3;
    //init
    fp12cv_init(&a1);
    fp12cv_init(&a2);
    fp12cv_init(&a3);
    fp12cv_init(&a4);
    fpm_init(&b1);
    fpm_init(&b2);
    fpm_init(&b3);
    //build matrix
    matrix_build_fpm_and_fp12cv();
    //matrix_printf("M = ",&matrix_of_fpm_to_fp12cv);
    //set random
    fp12cv_set_random(&a1,state);
    fp12cv_set_random(&a2,state);
    //check
    fp12cv_to_fpm(&b1,&a1);
    fpm_to_fp12cv(&a2,&b1);
    printf("first cmp = %d\n",fp12cv_cmp(&a1,&a2));
    fp12cv_set_random(&a1,state);
    fp12cv_set_random(&a2,state);
    //add
    fp12cv_add(&a3,&a1,&a2);
    fp12cv_to_fpm(&b1,&a1);
    fp12cv_to_fpm(&b2,&a2);
    fpm_add(&b3,&b1,&b2);
    fpm_to_fp12cv(&a4,&b3);
    mul_cmp = fp12cv_cmp(&a3,&a4);
    printf("add cmp = %d\n",mul_cmp);
    //mul
    fp12cv_mul(&a3,&a1,&a2);
    fp12cv_to_fpm(&b1,&a1);
    fp12cv_to_fpm(&b2,&a2);
    fpm_mul(&b3,&b1,&b2);
    fpm_to_fp12cv(&a4,&b3);
    printf("mul cmp = %d\n",fp12cv_cmp(&a3,&a4));
    //frobenius
    fp12cv_frobenius(&a3,&a1);
    fp12cv_to_fpm(&b1,&a1);
    fpm_frobenius(&b2,&b1);
    fpm_to_fp12cv(&a4,&b2);
    printf("frobenius cmp = %d\n",fp12cv_cmp(&a3,&a4));

}
void test_operation(){
    matrix_t a1,a2,a3,a4,e;
    matrix_init(&a1);
    matrix_init(&a2);
    matrix_init(&a3);
    matrix_init(&a4);
    matrix_set_unit(&e);

    matrix_set_random(&a1,state);
    matrix_inv(&a2,&a1);
    
    matrix_mul(&a3,&a1,&a2);
    matrix_mul(&a4,&a2,&a1);
    printf("matrix_cmp = %d\n",matrix_cmp(&a3,&a4));
}