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
    
    //test_Field(10000,10000,10000,100);
    //test_Frobenius_map();
    //test_skew_frobenius_map();
    //test_twist();
    //test_mod(10000000,100000);
    
    //test_EFp(10,10,10);
    //test_EFp2(10,10,10);
    //test_EFp12(10,10,10);

    //BN12_test_rational_point();
    //BN12_test_plain_ate_pairing();
    //BN12_test_opt_ate_pairing();
    //BN12_test_x_ate_pairing();
    //BN12_test_G1_SCM();
    //BN12_test_G2_SCM();
    //BN12_test_G3_EXP();
    //BN12_compare_pairings();

    //BLS12_test_rational_point();
    //BLS12_test_plain_ate_pairing();
    //BLS12_test_opt_ate_pairing(10);
    //BLS12_test_x_ate_pairing();
    //BLS12_test_G1_SCM(100);
    BLS12_test_G2_SCM(10);
    //BLS12_test_G3_EXP(10);

    //Fast_test_BLS12(100,100,100,100);
    //test_All();
    return 0;
}

