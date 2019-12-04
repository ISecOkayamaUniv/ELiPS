#include <ELiPS/scalar.h>

void scalar_init(mpz_t a){
	mpz_init(a);
}

void scalar_clear(mpz_t a){
	mpz_clear(a);
}

void scalar_random_prime(mpz_t a){
    mpz_urandomm(a,state,prime_z);
}

void scalar_add_prime(mpz_t ans, mpz_t a, mpz_t b){
    mpz_add(ans,a,b);
    mpz_mod(ans,ans,prime_z);
}

void scalar_sub_prime(mpz_t ans, mpz_t a, mpz_t b){
    mpz_sub(ans,a,b);
    mpz_mod(ans,ans,prime_z);
}

void scalar_mul_prime(mpz_t ans, mpz_t a, mpz_t b){
    mpz_mul(ans,a,b);
    mpz_mod(ans,ans,prime_z);
}

void scalar_inv_prime(mpz_t ans, mpz_t a){
    mpz_invert(ans,a,prime_z);
}


void scalar_random_order(mpz_t a){
    mpz_urandomm(a,state,order_z);
}
void scalar_add_order(mpz_t ans, mpz_t a, mpz_t b){
    mpz_add(ans,a,b);
    mpz_mod(ans,ans,order_z);
}

void scalar_sub_order(mpz_t ans, mpz_t a, mpz_t b){
    mpz_sub(ans,a,b);
    mpz_mod(ans,ans,order_z);
}

void scalar_mul_order(mpz_t ans, mpz_t a, mpz_t b){
    mpz_mul(ans,a,b);
    mpz_mod(ans,ans,order_z);
}

void scalar_inv_order(mpz_t ans, mpz_t a){
    mpz_invert(ans,a,order_z);
}