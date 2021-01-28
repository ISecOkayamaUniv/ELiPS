#include <ELiPS/bls_pairing.h>

void bls_optate_pairing_basic(fpm2_t *ANS, efpm2_t *p, efpm2_t *q){
    //Miller
    bls_miller_for_optate_basic(ANS,p,q);
    //bls_2i3_miller_for_optate(ANS,p,q);

    //Final exp
    //bls_final_exp_optimal(ANS,ANS);
    bls_final_exp(ANS,ANS);
}