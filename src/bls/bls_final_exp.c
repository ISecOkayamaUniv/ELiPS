#include <ELiPS/bls_final_exp.h>

void bls_final_exp_optimal(fpm2_t *ANS, fpm2_t *a){
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime_z,DEGREE_MN);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q(exp,exp,order_z);
    fpm2_pow(ANS,a,exp);
    mpz_clear(exp);
}

int Eulers_totient_function(int n){
    int i,count=0;
    for(i=1;i<n;i++){
        if(gcd(n,i)==1)
            count++;
    }
    return count;
}

int gcd(int a, int b){
    int r;
    while((r = a % b)!=0 ){
        a = b;
        b = r;
    }
    return b;
}

void bls_fpm2_pow_X(fpm2_t *ANS, fpm2_t *a){
    int i;
    fpm2_t tmp,a_inv;
    //set
    fpm2_set(&tmp,a);
    fpm2_frobenius_times(&a_inv,a,6);
    for(i=bls_X_length-1;i>=0;i--){
        switch(bls_X_binary[i]){
            case 0:
                fpm2_mul(&tmp,&tmp,&tmp);
                break;
            case 1:
                fpm2_mul(&tmp,&tmp,&tmp);
                fpm2_mul(&tmp,&tmp,a);
                break;
            case -1:
                fpm2_mul(&tmp,&tmp,&tmp);
                fpm2_mul(&tmp,&tmp,&a_inv);
                break;
            default:
                break;
        }
    }
    fpm2_set(ANS,&tmp);
}

void bls_fpm2_pow_X2(fpm2_t *ANS, fpm2_t *a){
    int i;
    fpm2_t tmp,a_inv;
    //set
    fpm2_set(&tmp,a);
    fpm2_frobenius_times(&a_inv,a,6);
    for(i=bls_X_length-2;i>=0;i--){
        switch(bls_X2_binary[i]){
            case 0:
                fpm2_mul(&tmp,&tmp,&tmp);
                break;
            case 1:
                fpm2_mul(&tmp,&tmp,&tmp);
                fpm2_mul(&tmp,&tmp,a);
                break;
            case -1:
                fpm2_mul(&tmp,&tmp,&tmp);
                fpm2_mul(&tmp,&tmp,&a_inv);
                break;
            default:
                break;
        }
    }
    fpm2_set(ANS,&tmp);
}

void bls_final_exp(fpm2_t *ANS, fpm2_t *a){
    static int e[48][24] = {
        {-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},//Φ1
        {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},//Φ2
        {1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},//Φ6
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},//Φ12
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0,0,0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},//Φ24
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0}//Φ48
    };
    int i;
    int euler = Eulers_totient_function(DEGREE_MN);

    fpm2_t f,tmp1,tmp2,tmp3,tmp4;
    //Easy part
    fpm2_frobenius_times(&tmp1,a,DEGREE_MN/2);
    fpm2_inv(&tmp2,a);
    fpm2_mul(&f,&tmp1,&tmp2);//f = f^(p^6-1)
    fpm2_set(&tmp1,&f);
    fpm2_frobenius_times(&tmp1,&tmp1,DEGREE_MN/6);
    fpm2_mul(&f,&tmp1,&f);
    
    //Hard part
    fpm2_t ramda[euler],ramda_inv;
    fpm2_mul(&tmp1,&f,&f);//tmp1 = f^2
    bls_fpm2_pow_X(&tmp2,&tmp1);//tmp2 = f^2X
    bls_fpm2_pow_X2(&tmp3,&tmp2);//tmp3 = f^(X^2)
    fpm2_frobenius_times(&tmp4,&f,6);//tmp4 = f^-1
    fpm2_mul(&tmp2,&tmp2,&tmp4);//tmp2 = f^(2X-1)
    fpm2_frobenius_times(&tmp2,&tmp2,6);//tmp2 = f^(-2X+1)
    fpm2_mul(&ramda[euler-1],&tmp2,&tmp3);//f^(X^2-2X+1)
    fpm2_frobenius_times(&ramda_inv,&ramda[euler-1],6);

    for(i=euler-2;i>=0;i--){
        switch(e[DEGREE_MN-1][i+1]){
            case 0:
                bls_fpm2_pow_X(&ramda[i],&ramda[i+1]);
                break;
            case 1:
                bls_fpm2_pow_X(&ramda[i],&ramda[i+1]);
                fpm2_mul(&ramda[i],&ramda[i],&ramda[euler-1]);
                break;
            case -1:
                bls_fpm2_pow_X(&ramda[i],&ramda[i+1]);
                fpm2_mul(&ramda[i],&ramda[i],&ramda_inv);
                break;
        }
    }
    fpm2_mul(&ramda[0],&ramda[0],&tmp1);
    fpm2_mul(&tmp1,&ramda[0],&f);
    //fpm2_set(&tmp1,&ramda[0]);
    for(i=1;i<euler;i++){
        fpm2_frobenius_times(&tmp2,&ramda[i],i);
        fpm2_mul(&tmp1,&tmp1,&tmp2);
    }
    fpm2_set(ANS,&tmp1);

}
