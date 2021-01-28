#include <ELiPS/define_cvma.h>
/*=================================================*/
/*bls */
/*=================================================*/

int bls_X_binary[SIZE_BLS];
int bls_X2_binary[SIZE_BLS];
int bls_X_length;
mpz_t efpm2_total;
fpm2_t twist_g_sqrt,twist_g_sqrt_inv;
fpm2_t twist_g_cbrt,twist_g_cbrt_inv;
fpm2_t twist_g_sxrt,twist_g_sxrt_inv;
/*=================================================*/
/*cvma TypeX*/
/*=================================================*/

//1次拡大
int h_degree_m;
int r_degree_m;
int n_ijk_degree_m[SIZE_ALL];

//2次拡大
int h_degree_m2;
int r_degree_m2;
int n_ijk_degree_m2[SIZE_ALL];

//2次拡大
int h_degree_2;
int r_degree_2;
int n_ijk_degree_2[SIZE_ALL];
