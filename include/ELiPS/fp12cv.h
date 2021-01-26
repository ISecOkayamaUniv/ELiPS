#ifndef INCLUDED_fp12cv_h
#define INCLUDED_fp12cv_h
//fp12_cvma
//can use p = 2,3(mod 5) && p = 5 (mod 7)
#include <ELiPS/fp4cv.h>

extern void fp12cv_init(fp12cv_t *a);
extern void fp12cv_printf(char *str, fp12cv_t *a);
extern void fp12cv_println(char *str, fp12cv_t *a);
//set
extern void fp12cv_set(fp12cv_t *b, fp12cv_t *a);
extern void fp12cv_set_ui(fp12cv_t *a, unsigned long int b);
extern void fp12cv_set_random(fp12cv_t *a, gmp_randstate_t state);
//add
extern void fp12cv_add(fp12cv_t *c, fp12cv_t *a, fp12cv_t *b);
//sub
extern void fp12cv_sub(fp12cv_t *c, fp12cv_t *a, fp12cv_t *b);
//mul
extern void fp12cv_mul(fp12cv_t *c, fp12cv_t *a, fp12cv_t *b);
//sqr
extern void fp12cv_sqr(fp12cv_t *b, fp12cv_t *a);
//frobenius
extern void fp12cv_frobenius(fp12cv_t *b, fp12cv_t *a);
//inv
extern void fp12cv_inv(fp12cv_t *b, fp12cv_t *a);
extern void fp12cv_pow(fp12cv_t *b, fp12cv_t *a, mpz_t scalar);
extern void fp12cv_sqrt(fp12cv_t *b, fp12cv_t *a);
extern int fp12cv_legendre(fp12cv_t *a);
//cmp
extern int fp12cv_cmp(fp12cv_t *a, fp12cv_t *b);
extern int fp12cv_cmp_ui(fp12cv_t *a, unsigned long int b);
extern int fp12cv_cmp_zero(fp12cv_t *a);
extern int fp12cv_cmp_one(fp12cv_t *a);

#endif