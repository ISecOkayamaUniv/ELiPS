#include "pairing_precompute.h"

int main(){
    //initialize ELiPS
    bls12_init();
    fr_order_init();

    //declear 
    g1_t P,aP;
    g2_t Q,bQ;
    g3_t ANS1,ANS2;
    fr_t a;

    //initialize
    g1_init(&P);
    g1_init(&aP);
    g2_init(&Q);
    g2_init(&bQ);
    g3_init(&ANS1);
    g3_init(&ANS2);
    fr_init(&a);

    //random
    g1_set_random(&P,state);
    g2_set_random(&Q,state);
    fr_set_random(&a,state);

    //compute scm
    g1_scm(&aP,&P,&a);

    //initialize precompute
    g1g2_to_g3_pairing_precompute_init(&Q);

    //pairing using precompute
    g1g2_to_g3_pairing_precompute(&ANS1,&P);

    //check linear
    g3_exp(&ANS1,&ANS1,&a);
    g1g2_to_g3_pairing_precompute(&ANS2,&aP);

    fp12_printf_montgomery("ANS1=",&ANS1);
    printf("\n");
    fp12_printf_montgomery("ANS2=",&ANS2);
    if(fp12_cmp(&ANS1,&ANS2)==0) printf("ok!\n");
    else printf("ng!\n");
    

    return 0;
}