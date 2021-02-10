#include "ELiPS/bls12.h"
int main(void){
    //Declarate variable
    scalar_t s1,s2,s12;
    efp12_t P,Q;
    efp12_t s1P,s2P,s1Q,s2Q;
    fp12_t ans1,ans2,ans3;
    cost tmp,basic_exp,cvma_exp,compless_exp;
    
    //Initialize variable
    scalar_init(s1);
    scalar_init(s2);
    scalar_init(s12);
    efp12_init(&P);
    efp12_init(&Q);
    fp12_init(&ans1);
    fp12_init(&ans2);
    fp12_init(&ans3);
    efp12_init(&s1P);
    efp12_init(&s2P);
    efp12_init(&s1Q);
    efp12_init(&s2Q);
    cost_init(&tmp);
    cost_init(&basic_exp);
    cost_init(&cvma_exp);
    cost_init(&compless_exp);

    //Initialize and print parameter
    bls12_init();
    bls12_print_parameters();
    
    //Generate random scalar
    scalar_random_order(s1);
    scalar_random_order(s2);
    scalar_mul_order(s12,s1,s2);
    
    //Generate rationalpoit P on G1 and Q on G2
    bls12_generate_g1(&P);
    bls12_generate_g2(&Q);
    
    //Scalar mutiplication
    bls12_g1_scm(&s1P,&P,s1);
    bls12_g1_scm(&s2P,&P,s2);
    bls12_g2_scm(&s1Q,&Q,s1);
    bls12_g2_scm(&s2Q,&Q,s2);
    
    fpm_build();
    //matrix_build();
    matrix_build_fp12_and_fp12cv();
    //Calculate pairing
    
    

    bls12_miller_algo_for_optate_basic(&ans1,&s1P,&s2Q);
    cost_zero();
    bls12_final_exp_optimal(&ans2,&ans1);
    cost_check(&tmp);
    cost_addition(&basic_exp,&tmp);
    cost_printf("basic exp",&basic_exp,1);

    cost_zero();
    bls12_final_exp_optimal_cvma(&ans3,&ans1);
    cost_check(&tmp);
    cost_addition(&cvma_exp,&tmp);
    cost_printf("cvma exp",&cvma_exp,1);

    cost_zero();
    bls12_final_exp_optimal_compress(&ans1,&ans1);
    cost_check(&tmp);
    cost_addition(&compless_exp,&tmp);
    cost_printf("compless exp",&compless_exp,1);
    //
    //
    cost_zero();
    cost_init(&tmp);
    bls12_fp12_pow_X(&ans1,&ans1);
    cost_check(&tmp);
    cost_printf("bls12_fp12_pow_X",&tmp,1);

    cost_zero();
    cost_init(&tmp);
    bls12_optate_pairing_projective(&ans1,&P,&Q);
    cost_check(&tmp);
    cost_printf("bls12_fp12 pairing",&tmp,1);

    cost_zero();
    cost_init(&tmp);
    bls12_optate_pairing_cvma(&ans1,&P,&Q);
    cost_check(&tmp);
    cost_printf("bls12_fp12cv pairing",&tmp,1);

    fp12cv_t cv;
    fp12_to_fp12cv(&cv,&ans1);
    cost_zero();
    cost_init(&tmp);
    bls12_fp12cv_pow_X(&cv,&cv);
    
    cost_check(&tmp);
    cost_printf("bls12_fp12cv_pow_X",&tmp,1);

    cost_zero();
    cost_init(&tmp);
    fp12_sqr_GS(&ans1,&ans1);
    cost_check(&tmp);
    cost_printf("fp12_sqr_GS",&tmp,1);

    cost_zero();
    cost_init(&tmp);
    fp12cv_sqr(&cv,&cv);
    cost_check(&tmp);
    cost_printf("fp12cv_sqr",&tmp,1);

    fp4cv_t fp4cv1;
    fp4cv_set_random(&fp4cv1,state);
    cost_zero();
    cost_init(&tmp);
    fp4cv_mul(&fp4cv1,&fp4cv1,&fp4cv1);
    cost_check(&tmp);
    cost_printf("fp4cv_mul",&tmp,1);


    bls12_optate_pairing_cvma(&ans1,&P,&Q);
    bls12_optate_pairing_cvma(&ans2,&s1P,&s2Q);
    bls12_optate_pairing_cvma(&ans3,&s2P,&s1Q);

    //Calculate exponentiation
    bls12_g3_exp(&ans1,&ans1,s12);
    
    //Print answer of pairing
	fp12_println("e(P,Q)^s12=",&ans1);
	fp12_println("e([s1]P,[s2]Q)=",&ans2);
	fp12_println("e([s2]P,[s1]Q)=",&ans3);
	
    //Compare answer
    if(fp12_cmp(&ans1,&ans2)==0)  printf("e(P,Q)^s12=e([s1]P,[s2]Q)\n");
    else    printf("e(P,Q)^s12!=e([s1]P,[s2]Q)\n");
    if(fp12_cmp(&ans1,&ans3)==0)  printf("e(P,Q)^s12=e([s2]P,[s1]Q)\n");
    else    printf("e(P,Q)^s12!=e([s2]P,[s1]Q)\n");
    if(fp12_cmp(&ans2,&ans3)==0)  printf("e([s1]P,[s2]Q)=e([s2]P,[s1]Q)\n\n");
    else    printf("e([s1]P,[s2]Q)!=e([s2]P,[s1]Q)\n\n");
    
    
    //Clear variable
    scalar_clear(s12);
    scalar_clear(s1);
    scalar_clear(s2);
    
    return 0;
}