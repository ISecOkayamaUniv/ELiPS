#include <gmp.h>
#include "ELiPS/bls12.h"
#include <ELiPS/fp12cv.h>
void test_fp12cv();
void test_fp12cv_GS();

int main(void){
    bls12_init();
    bls12_print_parameters();
    test_fp12cv_GS();
    return 0;
}
void test_fp12cv(){
    fp12cv_t a1,a2,a3,a4;
    fp12cv_init(&a1);
    fp12cv_init(&a2);
    fp12cv_init(&a3);
    fp12cv_init(&a4);
    //random
    fp12cv_set_random(&a1,state);
    fp12cv_set_random(&a2,state);
    fp12cv_set_random(&a3,state);
    fp12cv_set_random(&a4,state);
    //mul inv
    fp12cv_inv(&a2,&a1);
    fp12cv_mul(&a3,&a1,&a2);
    fp12cv_println("fp12cv a3 = ",&a3);
    printf("fp12cv cmp 1 = %d\n",fp12cv_cmp_one(&a3));
    fp12cv_mul(&a4,&a1,&a3);
    printf("a1 cmp a1*1  = %d\n",fp12cv_cmp(&a1,&a4));
    fp12cv_mul(&a3,&a1,&a2);
    fp12cv_mul(&a4,&a2,&a1);
    printf("fp12cv cmp a3a4 = %d\n",fp12cv_cmp(&a3,&a4));
    //
    fp12cv_set_ui(&a3,1);
    fp12cv_println("set ui 1 a3 = ",&a3);
    fp12cv_mul(&a4,&a1,&a3);
    printf("set ui a1 cmp a1*1  = %d\n",fp12cv_cmp(&a1,&a4));

    //sqr
    fp12cv_set_random(&a1,state);
    fp12cv_mul(&a3,&a1,&a1);
    fp12cv_sqr(&a2,&a1);
    printf("fp12cv cmp sqr = %d\n",fp12cv_cmp(&a3,&a2));
    //test pow X_z
    fp12cv_set_random(&a1,state);
    bls12_fp12cv_pow_X(&a2,&a1);
    fp12cv_pow(&a3,&a1,X_z);
    printf("fp12cv cmp X_z = %d\n",fp12cv_cmp(&a3,&a2));
}
void test_fp12cv_GS(){
    fp12cv_t t1,t2,t3;
    fp4cv_t a1,a2,a3,a4;
    fp12cv_set_random(&t1,state);
    fp12cv_set_random(&t2,state);
    fp12cv_set_random(&t3,state);

    bls12_pow_easypart_cvma(&t1,&t1);
    fp4cv_sub(&a1,&t1.x0,&t1.x1);
    fp4cv_sqr(&a2,&a1);
    fp4cv_mul(&a3,&t1.x0,&t1.x2);
    fp4cv_sub(&a1,&a2,&a3);
    fp4cv_frobenius_times(&a2,&t1.x1,2);
    printf("ccmp = %d\n",fp4cv_cmp(&a1,&a2));
}