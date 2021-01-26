#ifndef INCLUDED_define_cvma_h
#define INCLUDED_define_cvma_h

#include <ELiPS/define.h>

/*============================================================================*/
/* cvma                                                                       */
/*============================================================================*/

#define DEGREE_EXTENTION_FIELD_M 4
#define DEGREE_EXTENTION_FIELD_M2 3

//とりあえず定義している，修正の必要あり
#define SIZE_ALL 2048
#define SIZE_BLS 256


typedef struct{
  fp_t x[DEGREE_EXTENTION_FIELD_M];
}fpm_t;

typedef struct{
  fpm_t x[DEGREE_EXTENTION_FIELD_M2];
}fpm2_t;

typedef struct{
  fpm_t x;
  fpm_t y;
  int infinity;
}efpm_t;

typedef struct{
  fpm2_t x;
  fpm2_t y;
  int infinity;
}efpm2_t;

/*=================================================*/
/*bls */
/*=================================================*/
extern int bls_X_binary[SIZE_BLS];
extern mpz_t efpm2_total;
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


#endif
