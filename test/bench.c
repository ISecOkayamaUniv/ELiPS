#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"
/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    bls12_init();
    test_fp_montgomery(100);
    // g1_test(1000);
    // g2_test(1000);
    // g3_test(1000);
    // g1_set_random_test(1000);
    //debug_pairing(1000);

    //billinear_test();
    return 0;
}
