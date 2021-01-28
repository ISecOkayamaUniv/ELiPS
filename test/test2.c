#include "ELiPS/bls12.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    bls12_init();
    bls12_print_parameters();
    //test_field(0,1000,0,0,0);
    //test_fp(1000);
    //test_fp2(10000);
    //test_fp6(100);
    //test_mod(1000000);
    //test_Frobenius_map();
    //test_skew_frobenius_map();
    //test_twist();
    //test_mod(10000000,100000);

    //test_efp(100,100);
    //test_efp2(1000,1000);
    //test_efp12(10,10,10);

    //bls12_test_rational_point();
    fp2_t A,B,C;
    cost A_cost;
    int i;
    fp2_set_random(&A,state);
    fp2_set_random(&B,state);
    fp2_to_montgomery(&A,&A);
    fp2_to_montgomery(&B,&B);

    cost_zero(&A_cost);
    fp2_mul_lazy_montgomery(&C,&A,&B);
    cost_check(&A_cost);
    cost_printf("mul",&A_cost,1);

    cost_zero(&A_cost);
    fp2_sqr_lazy_montgomery(&C,&A);
    cost_check(&A_cost);
    cost_printf("sqr",&A_cost,1);

    cost_zero(&A_cost);
    fp2_inv_lazy_montgomery(&C,&A);
    cost_check(&A_cost);
    cost_printf("inv",&A_cost,1);

    //fp2_sqr_lazy_montgomery(&C,&A);

    //bls12_test_opt_ate_pairing(10000);
    //bls12_test_g1_scm(100);
    //bls12_test_g2_scm(1);
    //bls12_test_g3_exp(100);

    //test_All();
    return 0;
}
