#include <ELiPS/bls12_init.h>
/*----------------------------------------------------------------------------*/
//init/set/clear
void bls12_init(){
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    bls12_init_parameters();
    bls12_generate_X();
    if(bls12_generate_prime()==1 && bls12_generate_order()==1){
        bls12_generate_trace();
        bls12_weil();
        bls12_get_epsilon();
        bls12_set_basis();
        bls12_set_frobenius_constant();
        bls12_set_curve_parameter();
	    pre_montgomery();
        bls12_power_init();
        fr_order_init();
    }else{
        bls12_clear();
        printf("error : prime\nexit\n");
    }
}

void bls12_init_parameters(){
    int i,j;
    
    mpz_init(X_z);
    mpz_init(prime_z);
    mpz_init(order_z);
    mpz_init(trace_z);
    mpz_init(X_mod_order_z);
    
    mpz_init(efp_total);
    mpz_init(efp12_total);
    mpn_zero(curve_b,FPLIMB);
    
    for(i=0; i<bls12_X_length+1; i++){
        bls12_X_binary[i]=0;
    }
    
    for(i=0; i<bls12_X2_length+1; i++){
        bls12_X2_binary[i]=0;
    }
    
    mpn_zero(epsilon1,FPLIMB);
    mpn_zero(epsilon2,FPLIMB);
    fp2_init(&Alpha_1);
    fp2_init(&Alpha_1_inv);
    
    for(i=0; i<12; i++){
        for(j=0; j<6; j++){
            fp2_init(&frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            fp2_init(&skew_frobenius_constant[i][j]);
        }
    }
    
}

void bls12_print_parameters(){
    printf("====================================================================================\n");
    printf("bls12\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(X_z,2),X_z);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(prime_z,2),prime_z);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(order_z,2),order_z);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(trace_z,2),trace_z);
    
    gmp_printf("\nelliptic curve\n");
    gmp_printf("E:y^2=x^3+%Nu\n", curve_b,FPLIMB);
    
    gmp_printf("\nmodulo polynomial\n");
    gmp_printf("fp2 = fp[alpha]/(alpha^2+1)\n");
    gmp_printf("fp6 = fp2[beta]/(beta^3-(alpha+1))\n");
    gmp_printf("fp12= fp6[gamma]/(gamma^2-beta)\n");
}

void bls12_clear(){
    int i,j;
    
    mpz_clear(X_z);
    mpz_clear(prime_z);
    mpz_clear(order_z);
    mpz_clear(trace_z);
    
    mpz_clear(efp_total);
    mpz_clear(efp12_total);
    
}
