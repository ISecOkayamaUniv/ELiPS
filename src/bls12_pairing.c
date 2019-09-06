#include <ELiPS/bls12_pairing.h>

/*----------------------------------------------------------------------------*/
//bls12
void BLS12_Plain_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    BLS12_Miller_algo_for_plain_ate(ANS,P,Q);
    
    //Final Exp.
    BLS12_Final_exp_optimal(ANS,ANS);
}

void BLS12_Opt_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    BLS12_Miller_algo_for_opt_ate(ANS,P,Q);
    
    //Final Exp.
    BLS12_Final_exp_optimal(ANS,ANS);
}

void BLS12_Opt_ate_pairing_compress(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    BLS12_Miller_algo_for_opt_ate(ANS,P,Q);
    
    //Final Exp.
    BLS12_Final_exp_optimal_compress(ANS,ANS);
}

void BLS12_Opt_ate_pairing_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    BLS12_Miller_algo_for_opt_ate_lazy(ANS,P,Q);
    
    //Final Exp.
    BLS12_Final_exp_optimal_lazy(ANS,ANS);
}

void BLS12_Opt_ate_pairing_compress_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    BLS12_Miller_algo_for_opt_ate_lazy(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT+=timedifference_msec(tv_start,tv_end);
   
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    BLS12_Final_exp_optimal_compress_lazy(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT+=timedifference_msec(tv_start,tv_end);
}
void BLS12_Opt_ate_pairing_compress_lazy_montgomery(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    #ifdef DEBUG_COST_A
    cost tmp;
    #endif

    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    BLS12_Miller_algo_for_opt_ate_lazy_montgomery(ANS,P,Q);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT_MONTGOMERY+=timedifference_msec(tv_start,tv_end);
    
    #ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_MONTGOMERY_COST,&tmp);
    #endif
   
    //Final Exp.
    gettimeofday(&tv_start,NULL);
    BLS12_Final_exp_optimal_compress_lazy_montgomery(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_MONTGOMERY+=timedifference_msec(tv_start,tv_end);
}