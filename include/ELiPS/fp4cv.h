#ifndef INCLUDED_fp4cv_h
#define INCLUDED_fp4cv_h
//fp4_cvma
//can use p=2,3(mod 5)
#include <ELiPS/fp.h>

extern void fp4cv_init(fp4cv_t *a);
extern void fp4cv_printf(char *str, fp4cv_t *a);
extern void fp4cv_println(char *str, fp4cv_t *a);
extern void fpd4cv_println(char *str, fpd4cv_t *a);
//set
extern void fp4cv_set(fp4cv_t *b, fp4cv_t *a);
extern void fpd4cv_set(fpd4cv_t *b, fpd4cv_t *a);
extern void fp4cv_set_ui(fp4cv_t *a, unsigned long int b);
extern void fp4cv_set_mpn(fp4cv_t *b, mp_limb_t *a);
extern void fp4cv_set_neg(fp4cv_t *b, fp4cv_t *a);
extern void fp4cv_set_random(fp4cv_t *a, gmp_randstate_t state);
//montgomery
extern void fp4cv_to_montgomery(fp4cv_t *b, fp4cv_t *a);
//add
extern void fp4cv_add(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b);
extern void fp4cv_add_nonmod_single(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b);
extern void fp4cv_add_nonmod_double(fpd4cv_t *c, fpd4cv_t *a, fpd4cv_t *b);
extern void fp4cv_add_ui(fp4cv_t *c, fp4cv_t *a, unsigned long int b);
//sub
extern void fp4cv_sub(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b);
extern void fp4cv_sub_nonmod_single(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b);
extern void fp4cv_sub_nonmod_double(fpd4cv_t *c, fpd4cv_t *a, fpd4cv_t *b);
extern void fp4cv_sub_ui(fp4cv_t *c, fp4cv_t *a, unsigned long int b);
//mul
extern void fp4cv_mul(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b);
extern void fp4cv_mul_lazy(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b);
//extern void fp4cv_mul_lazy_montgomery(fp4cv_t *c, fp4cv_t *a, fp4cv_t *b);//bug
//sqr
extern void fp4cv_sqr(fp4cv_t *b, fp4cv_t *a);
extern void fp4cv_sqr_lazy(fp4cv_t *b, fp4cv_t *a);
//extern void fp4cv_sqr_lazy_montgomery(fp4cv_t *b, fp4cv_t *a);
//inv
extern void fp4cv_frobenius(fp4cv_t *b, fp4cv_t *a);
extern void fp4cv_frobenius_times(fp4cv_t *b, fp4cv_t *a, unsigned long int number);
extern void fp4cv_inv(fp4cv_t *b, fp4cv_t *a);
extern void fp4cv_pow(fp4cv_t *b, fp4cv_t *a, mpz_t scalar);
extern void fp4cv_sqrt(fp4cv_t *b, fp4cv_t *a);
extern int fp4cv_legendre(fp4cv_t *a);
//cmp
extern int fp4cv_cmp(fp4cv_t *a, fp4cv_t *b);
extern int fp4cv_cmp_ui(fp4cv_t *a, unsigned long int b);
extern int fp4cv_cmp_zero(fp4cv_t *a);
extern int fp4cv_cmp_one(fp4cv_t *a);

#endif
