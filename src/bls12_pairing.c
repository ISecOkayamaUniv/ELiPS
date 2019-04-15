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
    BLS12_Miller_algo_for_opt_ate_lazy(ANS,P,Q);
    
    //Final Exp.
    BLS12_Final_exp_optimal_compress_lazy(ANS,ANS);
}

void BLS12_X_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //Miller's Algo.
    BLS12_Miller_algo_for_x_ate(ANS,P,Q);
    
    //Final Exp.    
    BLS12_Final_exp_optimal(ANS,ANS);
}
