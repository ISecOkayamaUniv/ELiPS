#ifndef SCALAR_H
#define SCALAR_H

#include <ELiPS/define.h>

extern void scalar_init(mpz_t a);
extern void scalar_clear(mpz_t a);


extern void scalar_add(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_sub(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_mul(mpz_t ans, mpz_t a, mpz_t b);

extern void scalar_random_prime(mpz_t a);
extern void scalar_add_prime(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_sub_prime(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_mul_prime(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_inv_prime(mpz_t ans, mpz_t a);


extern void scalar_random_order(mpz_t a);
extern void scalar_add_order(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_sub_order(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_mul_order(mpz_t ans, mpz_t a, mpz_t b);
extern void scalar_inv_order(mpz_t ans, mpz_t a);

#endif
