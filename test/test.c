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
    //bls12_test_g1_scm(1);

    //bls12_test_rational_point();
    bls12_test_opt_ate_pairing(10);
    //bls12_test_symmmetric_opt_ate_pairing(10);
    //bls12_test_g1_scm(100);
    //bls12_test_g2_scm(1);
    //bls12_test_g3_exp(100);

    //test_All();
    return 0;
}
