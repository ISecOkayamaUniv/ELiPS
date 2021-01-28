#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"
/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    bls12_init();
    fr_order_init();
    //gmp_randinit_default (state);
    // bls12_print_parameters();
    // fr_order_init();
    //cost test;
    //cost_init(&test);
    //cost_mul=0;
    billlinear_test();
    //g1_test(10);
    //debug_pairing(10000);
    //cost_check(&test);
    //printf("sub=%d",cost_mul);
    return 0;
}
