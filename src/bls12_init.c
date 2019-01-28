#include <ELiPS/bls12_init.h>
/*----------------------------------------------------------------------------*/
//init/set/clear
void BLS12_init(){
    BLS12_init_parameters();
    BLS12_generate_X();
    if(BLS12_generate_prime()==1 && BLS12_generate_order()==1){
        BLS12_generate_trace();
        BLS12_weil();
        BLS12_get_epsilon();
        BLS12_get_Two_inv();
        BLS12_set_basis();
        BLS12_set_frobenius_constant();
        BLS12_set_curve_parameter();
        BLS12_set_root2();
    }else{
        BLS12_clear();
        printf("error : prime\nexit\n");
    }
}

void BLS12_init_parameters(){
    int i,j;
    
    mpz_init(X_z);
    mpz_init(prime_z);
    mpz_init(order_z);
    mpz_init(trace_z);
    mpz_init(root_X);
    mpz_init(root_2);
    
    mpz_init(EFp_total);
    mpz_init(EFp12_total);
    mpn_zero(curve_b,FPLIMB);
    
    for(i=0; i<BLS12_X_length+1; i++){
        BLS12_X_binary[i]=0;
    }
    
    for(i=0; i<BLS12_X2_length+1; i++){
        BLS12_X2_binary[i]=0;
    }
    
    mpn_zero(epsilon1,FPLIMB);
    mpn_zero(epsilon2,FPLIMB);
    mpz_init(Two_inv_z);
    Fp2_init(&Alpha_1);
    Fp2_init(&Alpha_1_inv);
    
    for(i=0; i<12; i++){
        for(j=0; j<6; j++){
            Fp2_init(&frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            Fp2_init(&skew_frobenius_constant[i][j]);
        }
    }
    
}

void BLS12_print_parameters(){
    printf("====================================================================================\n");
    printf("BLS12\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(X_z,2),X_z);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(prime_z,2),prime_z);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(order_z,2),order_z);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(trace_z,2),trace_z);
    
    gmp_printf("\nelliptic curve\n");
    gmp_printf("E:y^2=x^3+%Nu\n", curve_b,FPLIMB);
    
    gmp_printf("\nmodulo polynomial\n");
    gmp_printf("Fp2 = Fp[alpha]/(alpha^2+1)\n");
    gmp_printf("Fp6 = Fp2[beta]/(beta^3-(alpha+1))\n");
    gmp_printf("Fp12= Fp6[gamma]/(gamma^2-beta)\n");
}

void BLS12_clear(){
    int i,j;
    
    mpz_clear(X_z);
    mpz_clear(prime_z);
    mpz_clear(order_z);
    mpz_clear(trace_z);
    mpz_clear(root_X);
    mpz_clear(root_2);
    
    mpz_clear(EFp_total);
    mpz_clear(EFp12_total);
    
    mpz_clear(Two_inv_z);
}
