#ifndef INCLUDED_fpm_h
#define INCLUDED_fpm_h

#include <ELiPS/define_cvma.h>
#include <ELiPS/fp.h>

extern void fpm_make_cvma();
extern void fpm_preparate_mul();
extern void fpm_print_paramerter();
extern void fpm_build();
extern void fpm_init(fpm_t *a);
extern void fpm_printf(char *str,fpm_t *a);
extern void fpm_println(char *str,fpm_t *a);
//set
extern void fpm_set(fpm_t *b, fpm_t *a);
extern void fpm_set_ui(fpm_t *a, unsigned long int b);
extern void fpm_set_mpn(fpm_t *a, mp_limb_t *b);
extern void fpm_set_fp(fpm_t *a, fp_t *b);
extern void fpm_set_neg(fpm_t *b, fpm_t *a);
extern void fpm_set_random(fpm_t *a,gmp_randstate_t state);
//cmp
extern int fpm_cmp(fpm_t *a, fpm_t *b);
extern int fpm_cmp_ui(fpm_t *a, unsigned long int b);
//add
extern void fpm_add(fpm_t *c, fpm_t *a, fpm_t *b);
//sub
extern void fpm_sub(fpm_t *c, fpm_t *a, fpm_t *b);
//mul
extern void fpm_mul(fpm_t *c, fpm_t *a, fpm_t *b);
//sqr
extern void fpm_sqr(fpm_t *b, fpm_t *a);
extern void fpm_frobenius(fpm_t *b, fpm_t *a);
extern void fpm_inv(fpm_t *b, fpm_t *a);
extern void fpm_pow(fpm_t *b,fpm_t *a,mpz_t scalar);
extern int fpm_legendre(fpm_t *a);
extern void fpm_sqrt(fpm_t *b,fpm_t *a);

#endif
