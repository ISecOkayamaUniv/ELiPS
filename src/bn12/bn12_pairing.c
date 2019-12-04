#include <ELiPS/bn12_pairing.h>

/*----------------------------------------------------------------------------*/
//bn12
void bn12_ate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bn12_miller_algo_for_plain_ate(ANS,P,Q);
    
    //Final Exp.
    bn12_final_exp_optimal(ANS,ANS);
}

void bn12_optate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bn12_miller_algo_for_opt_ate(ANS,P,Q);
    
    //Final Exp.
    bn12_final_exp_optimal(ANS,ANS);
}


void bn12_optate_pairing_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bn12_miller_algo_for_opt_ate_lazy(ANS,P,Q);
    
    //Final Exp.
    bn12_final_exp_optimal_lazy(ANS,ANS);
}