#ifndef INCLUDED_define_cvma_h
#define INCLUDED_define_cvma_h

#include <ELiPS/define.h>

/*============================================================================*/
/* cvma                                                                       */
/*============================================================================*/

#define DEGREE_MN 12
#define DEGREE_EXTENTION_FIELD_M 12
#define DEGREE_EXTENTION_FIELD_M2 3
#define DEGREE_EXTENTION_FIELD_2 2

//とりあえず定義している，修正の必要あり
#define SIZE_ALL 2048
#define SIZE_BLS 256

typedef struct {
  fp_t x[DEGREE_EXTENTION_FIELD_M];
} fpm_t;

typedef struct {
  fpm_t x[DEGREE_EXTENTION_FIELD_M2];
} fpm2_t;

typedef struct {
  fp_t x[DEGREE_EXTENTION_FIELD_2];
} fp2cv_t;

typedef struct {
  fpm_t x;
  fpm_t y;
  int infinity;
} efpm_t;

typedef struct {
  fpm2_t x;
  fpm2_t y;
  int infinity;
} efpm2_t;

typedef struct {
  fp2cv_t x;
  fp2cv_t y;
  int infinity;
} efp2cv_t;

/*=================================================*/
/*bls 2^i*3*/
/*=================================================*/
extern int bls_X_binary[SIZE_BLS];
extern int bls_X2_binary[SIZE_BLS];
extern int bls_X_length;
extern mpz_t efpm2_total;
extern fpm2_t twist_g_sqrt, twist_g_sqrt_inv;
extern fpm2_t twist_g_cbrt, twist_g_cbrt_inv;
extern fpm2_t twist_g_sxrt, twist_g_sxrt_inv;

/*=================================================*/
/*cvma TypeX*/
/*=================================================*/

//1次拡大
extern int h_degree_m;
extern int r_degree_m;
extern int n_ijk_degree_m[SIZE_ALL];

//2次拡大
extern int h_degree_m2;
extern int r_degree_m2;
extern int n_ijk_degree_m2[SIZE_ALL];

//for twist
extern int h_degree_2;
extern int r_degree_2;
extern int n_ijk_degree_2[SIZE_ALL];

#endif
