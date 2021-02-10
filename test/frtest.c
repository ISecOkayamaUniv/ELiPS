#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"
/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    bls12_init();
    //gmp_randinit_default (state);
    // bls12_print_parameters();
    // fr_order_init();
    //cost test;
    //cost_init(&test);
    //cost_mul=0;
    //billlinear_test();
    //g2_test(1000);
    //g3_test(1000);

    //bls12_test_g2_scm(100);
    debug_pairing(1000);
    //cost_check(&test);
    //printf("sub=%d",cost_mul);
    return 0;
}
