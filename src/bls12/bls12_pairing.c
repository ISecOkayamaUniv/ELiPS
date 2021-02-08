#include <ELiPS/bls12_pairing.h>

/*----------------------------------------------------------------------------*/
//bls12
void bls12_optate_pairing_basic(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bls12_miller_algo_for_optate_basic(ANS,P,Q);
    
    //Final Exp.
    bls12_final_exp_optimal(ANS,ANS);
}
void bls12_optate_pairing_affine(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    #ifdef DEBUG_COST_A
    cost tmp;
    #endif

    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    bls12_miller_algo_for_optate_affine(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT_AFFINE+=timedifference_msec(tv_start,tv_end);
    
    #ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_AFFINE_COST,&tmp);
    #endif
   
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    bls12_final_exp_optimal_compress(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_AFFINE+=timedifference_msec(tv_start,tv_end);
}

void bls12_optate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    #ifdef DEBUG_COST_A
    cost tmp;
    #endif

    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    bls12_miller_algo_for_optate_projective(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT_PROJECTIVE+=timedifference_msec(tv_start,tv_end);
    
    #ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_PROJECTIVE_COST,&tmp);
    #endif
   
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    bls12_final_exp_optimal_compress(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_PROJECTIVE+=timedifference_msec(tv_start,tv_end);
}

void bls12_symmetric_optate_pairing(fp12_t *ANS,sym_t *A,sym_t *B){
    efp12_t P;
    efp12_t Q;
    
    efp12_set(&P,&A->p);
    efp12_set(&Q,&B->q);
    
    //Miller's Algo.
    bls12_miller_algo_for_optate_projective(ANS,&P,&Q);
    
    //Final Exp.
    bls12_final_exp_optimal_compress(ANS,ANS);
}

void bls12_optate_pairing_basic_cvma(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bls12_miller_algo_for_optate_basic(ANS,P,Q);
    
    //Final Exp.
    bls12_final_exp_optimal_cvma(ANS,ANS);
}
