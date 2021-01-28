#include <ELiPS/bls12_generate_curve.h>

void bls12_generate_X(){
    int i;
    mpz_t buf;
    mpz_init(buf);

    //X_binary
    #if defined(PARAM_ELIPS) || defined(PARAM_ELIPS_CHANGE_B)
    bls12_X_binary[77]=1;
    bls12_X_binary[11]=1;
    bls12_X_binary[9]=-1;
    bls12_X_binary[6]=-1;

    bls12_X2_binary[76]=1;
    bls12_X2_binary[10]=1;
    bls12_X2_binary[8]=-1;
    bls12_X2_binary[5]=-1;
    #endif
    //#ifdef PARAM_TAXONOMY
    #if defined(PARAM_TAXONOMY) || defined(PARAM_TAXONOMY_CHANGE_B)
    bls12_X_binary[77]=-1;
    bls12_X_binary[50]=1;
    bls12_X_binary[33]=1;

    bls12_X2_binary[76]=-1;
    bls12_X2_binary[49]=1;
    bls12_X2_binary[32]=1;
    #endif
    #ifdef PARAM_MINUS_TYPE2
    bls12_X_binary[77]=-1;
    bls12_X_binary[71]=1;
    bls12_X_binary[63]=1;
    bls12_X_binary[22]=-1;

    bls12_X2_binary[76]=-1;
    bls12_X2_binary[70]=1;
    bls12_X2_binary[62]=1;
    bls12_X2_binary[21]=-1;
    #endif
    #ifdef PARAM_PLUS_TYPE1
    bls12_X_binary[77]=1;
    bls12_X_binary[59]=1;
    bls12_X_binary[48]=-1;
    bls12_X_binary[22]=1;

    bls12_X2_binary[76]=1;
    bls12_X2_binary[58]=1;
    bls12_X2_binary[47]=-1;
    bls12_X2_binary[21]=1;
    #endif
    //#ifdef PARAM_ELIPS_IMPROVE_CHANGE_B
    #if defined(PARAM_ELIPS_IMPROVE) || defined(PARAM_ELIPS_IMPROVE_CHANGE_B)
    bls12_X_binary[77]=1;
    bls12_X_binary[10]=1;
    bls12_X_binary[9]=1;
    bls12_X_binary[6]=-1;

    bls12_X2_binary[76]=1;
    bls12_X2_binary[9]=1;
    bls12_X2_binary[8]=1;
    bls12_X2_binary[5]=-1;
    #endif
    //tatara param
    #ifdef PARAM_TATARA_CHANGE_B
    bls12_X_binary[77]=-1;
    bls12_X_binary[55]=-1;
    bls12_X_binary[32]=1;
    bls12_X_binary[8]=1;

    bls12_X2_binary[76]=-1;
    bls12_X2_binary[54]=-1;
    bls12_X2_binary[31]=1;
    bls12_X2_binary[7]=1;
    #endif
    //GOD param
    #ifdef PARAM_GOD
    bls12_X_binary[76]=1;
    bls12_X_binary[70]=-1;
    bls12_X_binary[66]=1;

    bls12_X2_binary[75]=1;
    bls12_X2_binary[69]=-1;
    bls12_X2_binary[65]=1;
    #endif

    /*
    //X_binary
    bls12_X_binary[76]=1;
    bls12_X_binary[59]=-1;
    bls12_X_binary[11]=-1;
    bls12_X_binary[5]=-1;
    */
    /*
    //X2_binary
    bls12_X2_binary[75]=1;
    bls12_X2_binary[58]=-1;
    bls12_X2_binary[10]=-1;
    bls12_X2_binary[4]=-1;
    */

    //bls12.X
    mpz_set_ui(X_z,0);
    for(i=bls12_X_length; i>=0; i--){
        if(bls12_X_binary[i]==1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_add(X_z,X_z,buf);
        }else if(bls12_X_binary[i]==-1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_sub(X_z,X_z,buf);
        }
    }
}

int bls12_generate_prime(){
    mpz_t result,buf1,buf2,modtest;
    mpz_init(result);
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(modtest);

    mpn_zero(prime, FPLIMB);
    mpn_zero(prime2, FPLIMB2);


    mpz_sub_ui(result,X_z,1);
    mpz_pow_ui(result,result,2);

    mpz_pow_ui(buf1,X_z,4);
    mpz_pow_ui(buf2,X_z,2);
    mpz_sub(buf1,buf1,buf2);
    mpz_add_ui(buf1,buf1,1);

    mpz_mul(result,result,buf1);

    //check div3
    mpz_mod_ui(modtest,result,3);
    if(mpz_cmp_ui(modtest,0)!=0){
        mpz_clear(result);
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(modtest);
        //printf("cannot devided by 3\n");
        return 0;
    }

    mpz_tdiv_q_ui(result,result,3);
    mpz_add(result,result,X_z);

    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(result);//init->clear
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(modtest);
        //printf("not prime\n");
        return 0;
    }

    mpz_set(prime_z,result);
    mpn_set_mpz(prime,result);
    mpn_mul_n(prime2,prime,prime,FPLIMB);

    mpz_clear(result);//init->clear
    mpz_clear(buf1);
    mpz_clear(buf2);
    mpz_clear(modtest);
    return 1;
}

int bls12_generate_order(){
    mpz_t buf1,buf2,result;
    mpz_init(result);//err result init
    mpz_init(buf1);
    mpz_init(buf2);

    mpz_pow_ui(buf1,X_z,4);
    mpz_pow_ui(buf2,X_z,2);
    mpz_sub(result,buf1,buf2);
    mpz_add_ui(result,result,1);

    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(result);
        mpz_clear(buf1);
        mpz_clear(buf2);
        //printf("not prime\n");
        return 0;
    }
    mpz_set(order_z,result);
    return 1;

    mpz_clear(buf1);
    mpz_clear(buf2);
    mpz_clear(result);
}

void bls12_generate_trace(){
    mpz_add_ui(trace_z,X_z,1);
}

void bls12_weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);

    //efp_total
    mpz_add_ui(buf,prime_z,1);
    mpz_sub(efp_total,buf,trace_z);

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

    //efp12_232_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(efp12_total,buf,1);

    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}
