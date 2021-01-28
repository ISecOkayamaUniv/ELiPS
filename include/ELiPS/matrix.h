#ifndef INCLUDED_matrix_h
#define INCLUDED_matrix_h

#include <ELiPS/fp12.h>
#include <ELiPS/fp12cv.h>
#include <ELiPS/fpm2.h>

extern void matrix_init(matrix_t *a);
extern void matrix_set(matrix_t *b, matrix_t *a);
extern void matrix_set_unit(matrix_t *a);
extern void matrix_set_random(matrix_t *a, gmp_randstate_t state);
extern int matrix_cmp(matrix_t *a, matrix_t *b);
extern void matrix_printf(char *str, matrix_t *a);
extern void matrix_mul(matrix_t *c, matrix_t *a, matrix_t *b);
extern void matrix_inv(matrix_t *b, matrix_t *a);

extern void matrix_build_fp12_and_fp12cv();
extern void matrix_build_fp12_and_fpm2();
extern void fp12_to_fp12cv(fp12cv_t *b, fp12_t *a);
extern void fp12cv_to_fp12(fp12_t *b ,fp12cv_t *a);
extern void fp12_to_fpm2(fpm2_t *b, fp12_t *a);
extern void fpm2_to_fp12(fp12_t *b, fpm2_t *a);
/*----------test----------*/
extern void matrix_build_fp12_and_fpm();
extern void fp12_to_fpm(fpm_t *b, fp12_t *a);
extern void fpm_to_fp12(fp12_t *a, fpm_t *b);
extern void matrix_build_fpm_and_fp12cv();
extern void fpm_to_fp12cv(fp12cv_t *b, fpm_t *a);
extern void fp12cv_to_fpm(fpm_t *b, fp12cv_t *a);
extern void matrix_build();
#endif