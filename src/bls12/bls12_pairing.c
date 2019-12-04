#include <ELiPS/bls12_pairing.h>

/*----------------------------------------------------------------------------*/
//bls12
void bls12_ate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bls12_miller_algo_for_plain_ate(ANS,P,Q);
    
    //Final Exp.
    bls12_final_exp_optimal(ANS,ANS);
}

void bls12_optate_pairing_basic(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bls12_miller_algo_for_opt_ate(ANS,P,Q);
    
    //Final Exp.
    bls12_final_exp_optimal(ANS,ANS);
}

void bls12_optate_pairing_compress(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bls12_miller_algo_for_opt_ate(ANS,P,Q);
    
    //Final Exp.
    bls12_final_exp_optimal_compress(ANS,ANS);
}

void bls12_optate_pairing_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    //Miller's Algo.
    bls12_miller_algo_for_opt_ate_lazy(ANS,P,Q);
    
    //Final Exp.
    bls12_final_exp_optimal_lazy(ANS,ANS);
}

void bls12_optate_pairing_compress_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q){
#ifdef DEBUG_COST_A
    cost tmp;
#endif
    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    bls12_miller_algo_for_opt_ate_lazy(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT+=timedifference_msec(tv_start,tv_end);
   
#ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_COST, &tmp);
#endif
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    bls12_final_exp_optimal_compress_lazy(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT+=timedifference_msec(tv_start,tv_end);
}
void bls12_optate_pairing_projective_compress_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q){
#ifdef DEBUG_COST_A
    cost tmp;
#endif
    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    bls12_miller_algo_for_opt_ate_projective_lazy(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT_PROJECTIVE+=timedifference_msec(tv_start,tv_end);
   
#ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_PROJECTIVE_COST, &tmp);
#endif
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    bls12_final_exp_optimal_compress_lazy(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_PROJECTIVE+=timedifference_msec(tv_start,tv_end);
}
void bls12_optate_pairing_compress_lazy_montgomery(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    #ifdef DEBUG_COST_A
    cost tmp;
    #endif

    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    bls12_miller_algo_for_opt_ate_lazy_montgomery(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT_MONTGOMERY+=timedifference_msec(tv_start,tv_end);
    
    #ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_MONTGOMERY_COST,&tmp);
    #endif
   
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    bls12_final_exp_optimal_compress_lazy_montgomery(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_MONTGOMERY+=timedifference_msec(tv_start,tv_end);
}
void bls12_optate_pairing_projective_compress_lazy_montgomery(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    #ifdef DEBUG_COST_A
    cost tmp;
    #endif

    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    bls12_miller_algo_for_opt_ate_projective_lazy_montgomery(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT_PROJECTIVE_MONTGOMERY+=timedifference_msec(tv_start,tv_end);
    
    #ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_PROJECTIVE_MONTGOMERY_COST,&tmp);
    #endif
   
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    bls12_final_exp_optimal_compress_lazy_montgomery(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_PROJECTIVE_MONTGOMERY+=timedifference_msec(tv_start,tv_end);
}