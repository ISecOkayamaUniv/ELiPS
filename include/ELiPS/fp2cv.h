#ifndef INCLUDED_fp2cv_h
#define INCLUDED_fp2cv_h

#include <ELiPS/define_cvma.h>
#include <ELiPS/fp.h>

extern void fp2cv_make_cvma();
extern void fp2cv_preparate_mul();
extern void fp2cv_print_paramerter();
extern void fp2cv_build();
extern void fp2cv_init(fp2cv_t *a);
extern void fp2cv_printf(char *str, fp2cv_t *a);
extern void fp2cv_println(char *str, fp2cv_t *a);
//set
extern void fp2cv_set(fp2cv_t *b, fp2cv_t *a);
extern void fp2cv_set_ui(fp2cv_t *a, unsigned long int b);
extern void fp2cv_set_mpn(fp2cv_t *a, mp_limb_t *b);
extern void fp2cv_set_fp(fp2cv_t *a, fp_t *b);
extern void fp2cv_set_neg(fp2cv_t *b, fp2cv_t *a);
extern void fp2cv_set_random(fp2cv_t *a, gmp_randstate_t state);
//cmp
extern int fp2cv_cmp(fp2cv_t *a, fp2cv_t *b);
extern int fp2cv_cmp_ui(fp2cv_t *a, unsigned long int b);
//add
extern void fp2cv_add(fp2cv_t *c, fp2cv_t *a, fp2cv_t *b);
//sub
extern void fp2cv_sub(fp2cv_t *c, fp2cv_t *a, fp2cv_t *b);
//mul
extern void fp2cv_mul(fp2cv_t *c, fp2cv_t *a, fp2cv_t *b);
//sqr
extern void fp2cv_sqr(fp2cv_t *b, fp2cv_t *a);
extern void fp2cv_frobenius(fp2cv_t *b, fp2cv_t *a);
extern void fp2cv_inv(fp2cv_t *b, fp2cv_t *a);
extern void fp2cv_pow(fp2cv_t *b, fp2cv_t *a, mpz_t scalar);
extern int fp2cv_legendre(fp2cv_t *a);
extern int fp2cv_legendre3(fp2cv_t *a);
extern int fp2cv_legendre6(fp2cv_t *a);
extern void fp2cv_sqrt(fp2cv_t *b, fp2cv_t *a);

#endif