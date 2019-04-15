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
    //test_mod(10000000,100000);
    
    //test_Field_Lazy(100000,100000,10000);
    //test_EFp_Lazy();
    //test_EFp2_Lazy();
    //test_EFp12_Lazy();

    //test_BN12_G1_SCM_Lazy();
    //test_BN12_G2_SCM_Lazy();
    //test_BN12_G3_exp_Lazy();
    //test_BN12_opt_ate_pairing_Lazy();

    //test_BLS12_G1_SCM_Lazy(1000);
    //test_BLS12_G2_SCM_Lazy(10);
    //test_BLS12_G3_exp_Lazy(10);
    //test_BLS12_opt_ate_pairing_Lazy(10);

    //test_EFp_Jacobian(1,1,1);
    //test_EFp2_Jacobian(1,1,10);
    //test_BLS12_G1_SCM_Jacobian(100);
    //test_BLS12_G2_SCM_Jacobian(10);

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

    //Fast_test_BLS12_G1_SCM(1000);
    //Fast_test_BLS12_G2_SCM(100);
    //Fast_test_BLS12_G3_EXP(100);
    //Fast_test_BLS12_pairing(100);
    
    //test_compressed(10);
    test_BLS12_opt_ate_pairing_compress(1);
    //test_All();
    return 0;
}

