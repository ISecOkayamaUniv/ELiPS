#include <ELiPS/bn12_generate_curve.h>

void BN12_generate_X(){
    int i;
    mpz_t buf;
    mpz_init(buf);
    
    //X_binary
    BN12_X_binary[114]=1;
    BN12_X_binary[84]=1;
    BN12_X_binary[53]=-1;
    BN12_X_binary[0]=-1;
    
    //X_binary_opt
    BN12_X6_2_binary[116]=1;
    BN12_X6_2_binary[115]=1;
    BN12_X6_2_binary[86]=1;
    BN12_X6_2_binary[85]=1;
    BN12_X6_2_binary[55]=-1;
    BN12_X6_2_binary[54]=-1;
    BN12_X6_2_binary[2]=-1;
    
    //BN12.X
    mpz_set_ui(X_z,0);
    for(i=BN12_X_length; i>=0; i--){
        if(BN12_X_binary[i]==1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_add(X_z,X_z,buf);
        }else if(BN12_X_binary[i]==-1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_sub(X_z,X_z,buf);
        }
    }
    
    mpz_clear(buf);
}

int BN12_generate_prime(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,X_z,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X_z,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X_z,2);
    mpz_mul_ui(buf,buf,24);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,X_z,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    mpn_set_mpz(prime,result);
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(prime_z,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

int BN12_generate_order(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,X_z,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X_z,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X_z,2);
    mpz_mul_ui(buf,buf,18);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,X_z,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(order_z,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

void BN12_generate_trace(){
    mpz_t buf;
    mpz_init(buf);
    
    mpz_pow_ui(buf,X_z,2);
    mpz_mul_ui(buf,buf,6);
    mpz_add_ui(trace_z,buf,1);
    
    mpz_clear(buf);
}

void BN12_weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,prime_z,1);
    mpz_sub(EFp_total,buf,trace_z);
    
    //t2←α^2+β^2
    mpz_pow_ui(t2,trace_z,2);
    mpz_mul_ui(buf,prime_z,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,prime_z,2);
    
    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);
    
    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    
    //EFp12_232_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(EFp12_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}
