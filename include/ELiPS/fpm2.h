#ifndef INCLUDED_fpm2_h
#define INCLUDED_fpm2_h

#include <ELiPS/fpm.h>

extern void fpm2_make_cvma();
extern void fpm2_preparate_mul();
extern void fpm2_print_paramerter();
extern void fpm2_build();
extern void fpm2_init(fpm2_t *a);
extern void fpm2_printf(char *str,fpm2_t *a);
extern void fpm2_println(char *str,fpm2_t *a);
extern void fpm2_set_random(fpm2_t *a,gmp_randstate_t state);
extern void fpm2_set(fpm2_t *b, fpm2_t *a);
extern void fpm2_set_ui(fpm2_t *a, unsigned long int b);
extern void fpm2_set_fp(fpm2_t *a, fp_t *b);
extern void fpm2_set_mpn(fpm2_t *a, mp_limb_t *b);
extern void fpm2_set_fpm(fpm2_t *b, fpm_t *a);
extern void fpm2_set_fp2cv(fpm2_t *b, fp2cv_t *a);
extern void fpm2_set_neg(fpm2_t *b, fpm2_t *a);
extern int fpm2_cmp(fpm2_t *a, fpm2_t *b);
extern int fpm2_cmp_ui(fpm2_t *a, unsigned long int b);
//fpm_cmp_one,zeroの実装
extern void fpm2_add(fpm2_t *c, fpm2_t *a, fpm2_t *b);
extern void fpm2_sub(fpm2_t *c, fpm2_t *a, fpm2_t *b);
extern void fpm2_mul(fpm2_t *c, fpm2_t *a, fpm2_t *b);
extern void fpm2_sqr(fpm2_t *b, fpm2_t *a);
extern void fpm2_frobenius(fpm2_t *b, fpm2_t *a);
extern void fpm2_frobenius_times(fpm2_t *b, fpm2_t *a, unsigned long int number);
extern void fpm2_inv(fpm2_t *b, fpm2_t *a);
extern void fpm2_pow(fpm2_t *b,fpm2_t *a,mpz_t scalar);
extern int  fpm2_legendre(fpm2_t *a);
extern int  fpm2_isCNR(fpm2_t *a);
extern void fpm2_sqrt(fpm2_t *b,fpm2_t *a);
extern void fpm2_cbrt(fpm2_t *b,fpm2_t *a);

#endif
