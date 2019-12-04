#include <ELiPS/bn12_init.h>

//init/set/clear
void bn12_init(){
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    bn12_init_parameters();
    bn12_generate_X();
    if(bn12_generate_prime()==1 && bn12_generate_order()==1){
        bn12_generate_trace();
        bn12_weil();
        bn12_get_epsilon();
        bn12_get_Two_inv();
        bn12_set_basis();
        bn12_set_frobenius_constant();
        bn12_set_curve_parameter();
    }else{
        bn12_clear();
        printf("error : prime\nexit\n");
    }
}
void bn12_init_parameters(){
    int i,j;
    mpz_init(X_z);
    mpz_init(prime_z);
    mpz_init(order_z);
    mpz_init(trace_z);
    
    mpz_init(efp_total);
    mpz_init(efp12_total);
    mpn_init(curve_b,FPLIMB);
    
    
    for(i=0; i<bn12_X_length+1; i++){
        bn12_X_binary[i]=0;
    }
    for(i=0; i<bn12_X6_2_length+1; i++){
        bn12_X6_2_binary[i]=0;
    }
    
    mpn_init(epsilon1,FPLIMB);
    mpn_init(epsilon2,FPLIMB);
    mpn_init(Two_inv,FPLIMB);
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
void bn12_print_parameters(){
    mpz_t mod_12;
    mpz_init(mod_12);
    printf("====================================================================================\n");
    printf("bn12 Class 2 (HW=4)\n\n");
    
    gmp_printf("Parameters\n");
    
    gmp_printf("X = 2^{114}+2^{84}-2^{53}-1\n");
    
    mpz_mod_ui(mod_12,X_z,12);
    gmp_printf("X mod 12 = %Zd\n\n",mod_12);
    
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(X_z,2),X_z);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(prime_z,2),prime_z);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(order_z,2),order_z);
    gmp_printf("trace (%dbit length) : %Zd \n\n",(int)mpz_sizeinbase(trace_z,2),trace_z);
    
    gmp_printf("BN curve\n");
    gmp_printf("E:y^2=x^3+%Nu\n", curve_b,FPLIMB);
    
    gmp_printf("Twisted curve\n");
    gmp_printf("E':y^2=x^3+2(alpha+1)^{-1}\n\n");
    
    gmp_printf("Extension field\n");
    gmp_printf("fp2 = fp[alpha]/(alpha^2+1)\n");
    gmp_printf("fp6 = fp2[beta]/(beta^3-(alpha+1))\n");
    gmp_printf("fp12= fp6[gamma]/(gamma^2-beta)\n");
    
    mpz_clear(mod_12);
}

void bn12_clear(){
    mpz_clear(X_z);
    mpz_clear(prime_z);
    mpz_clear(order_z);
    mpz_clear(trace_z);
    
    mpz_clear(efp_total);
    mpz_clear(efp12_total);
}
