#include "ELiPS/bn12.h"
#include "ELiPS/bls12.h"


/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    
    //BN12_init();
    //BN12_print_parameters();
    BLS12_init();
    BLS12_print_parameters();
    
    //test_Field();
    //test_Frobenius_map();
    //test_skew_frobenius_map();
    //BN12_test_rational_point();
    //BLS12_test_rational_point();
    //test_twist();

    
    //test_Field_Lazy(10000,10000,10000);
    //test_EFp_Lazy();
    //test_EFp2_Lazy();
    //test_EFp12_Lazy();

    //test_BN12_G1_SCM_Lazy();
    //test_BN12_G2_SCM_Lazy();
    //test_BN12_G3_exp_Lazy();
    //test_BN12_opt_ate_pairing_Lazy();

    //test_BLS12_G1_SCM_Lazy();
    //test_BLS12_G2_SCM_Lazy();
    //test_BLS12_G3_exp_Lazy();
    //test_BLS12_opt_ate_pairing_Lazy();

    //test_EFp_Jacobian();
    //test_EFp2_Jacobian();
    //test_BLS12_G1_SCM_Jacobian();
    //test_BLS12_G2_SCM_Jacobian();

    //BN12_test_plain_ate_pairing();
    //BN12_test_opt_ate_pairing();
    //BN12_test_x_ate_pairing();
    //BN12_test_G1_SCM();
    //BN12_test_G2_SCM();
    //BN12_test_G3_EXP();
    //BN12_compare_pairings();

    //BLS12_test_plain_ate_pairing();
    //BLS12_test_opt_ate_pairing();
    //BLS12_test_x_ate_pairing();
    //BLS12_test_G1_SCM();
    //BLS12_test_G2_SCM();
    //BLS12_test_G3_EXP();
    //BLS12_compare_pairings();

    Fast_test_BLS12_G1_SCM(100);
    Fast_test_BLS12_G2_SCM(1);
    Fast_test_BLS12_G3_EXP(1);
    Fast_test_BLS12_pairing(1);

    //test_All();
    return 0;
}

