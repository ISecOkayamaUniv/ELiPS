#ifndef DEFINE_H
#define DEFINE_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

/*************OPTION*************/
//bit
#define X64

//debug
//#define DEBUG_COST_A
//#define DEBUG_ASSERT

/************************************/

#ifdef X64
#define FPLIMB_BITS 512
#define FPLIMB 8
#define FPLIMB2 16
#endif

#ifdef X32
#define FPLIMB_BITS 480
#define FPLIMB 15
#define FPLIMB2 30
#endif

#define scalar_t mpz_t


#define bn12_X_length 114
#define bn12_X6_2_length 116

#define bls12_X_length 77
#define bls12_X2_length 76

//#define bls12_X_length 76
//#define bls12_X2_length 75


extern int bn12_X_binary[bn12_X_length+1];
extern int bn12_X6_2_binary[bn12_X6_2_length+1];

extern int bls12_X_binary[bls12_X_length+1];
extern int bls12_X2_binary[bls12_X2_length+1];

//extern int cost_add,cost_add_ui,cost_sub,cost_sub_ui,cost_mul,cost_mul_ui,cost_sqr,cost_inv,cost_mod;

extern int cost_add,cost_add_ui,cost_sub,cost_sub_ui,cost_mul,cost_set_neg,cost_div2,cost_sqr,cost_inv,cost_mod;
extern int cost_add_nonmod,cost_add_nonmod_double,cost_sub_nonmod,cost_sub_nonmod_double,cost_mod_nomal;

// typedef struct{
// 	int add;
//     int add_ui;
//     int sub;
//     int sub_ui;
//     int mul;
//     int mul_ui;
//     int sqr;
//     int inv;
//     int mod;
// }cost;
typedef struct{
	int add;
	int add_ui;
	int add_nonmod;
	int add_nonmod_double;
	int sub;
	int sub_ui;
	int sub_nonmod;
	int sub_nonmod_double;
	int mul;
	int set_neg;
	int div2;
	int sqr;
	int inv;
	int mod;
	int mod_nomal;
}cost;


/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/

typedef struct{
	mp_limb_t x0[FPLIMB];
}fp_t;
typedef struct{
	fp_t x0;
	fp_t x1;
}fp2_t;
typedef struct{
    fp2_t x0,x1,x2;
}fp6_t;
typedef struct{
    fp6_t x0,x1;
}fp12_t;

//tmp finite field
extern mp_limb_t buf[FPLIMB];

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/

typedef struct{
	mp_limb_t x0[FPLIMB2];
}fpd_t;
typedef struct{
	fpd_t x0;
	fpd_t x1;
}fpd2_t;
typedef struct{
    fpd2_t x0,x1,x2;
}fpd6_t;
typedef struct{
    fpd6_t x0,x1;
}fpd12_t;

/*============================================================================*/
/* Elliptic Curve                                                             */
/*============================================================================*/
typedef struct{
	fp_t x,y;
	int infinity;
}efp_t;

typedef struct{
    fp2_t x,y;
	int infinity;
}efp2_t;

typedef struct{
    fp6_t x,y;
	int infinity;
}efp6_t;

typedef struct{
    fp12_t x,y;
	int infinity;
}efp12_t;

/*============================================================================*/
/* Jacobian Elliptic Curve                                                   */
/*============================================================================*/
typedef struct{
	fp_t x,y,z;
	int infinity;
}efp_projective_t;

typedef struct{
	fp2_t x,y,z;
	int infinity;
}efp2_projective_t;

typedef struct{
	fp6_t x,y,z;
	int infinity;
}efp6_projective_t;

typedef struct{
	fp12_t x,y,z;
	int infinity;
}efp12_projective_t;

/*============================================================================*/
/* Jacobian Elliptic Curve                                                   */
/*============================================================================*/
typedef struct{
	fp_t x,y,z;
	int infinity;
}efp_jacobian_t;

typedef struct{
	fp2_t x,y,z;
	int infinity;
}efp2_jacobian_t;

typedef struct{
	fp6_t x,y,z;
	int infinity;
}efp6_jacobian_t;

typedef struct{
	fp12_t x,y,z;
	int infinity;
}efp12_jacobian_t;

/*============================================================================*/
/* rational point for symmetric pairing                                                             */
/*============================================================================*/
typedef struct{
	efp12_t p;
	efp12_t q;
	int infinity;
}sym_t;

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/
enum f_state{
	f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12
};
extern gmp_randstate_t state;
extern mpz_t X_z,prime_z,order_z,trace_z;
extern mp_limb_t X,prime[FPLIMB],order[FPLIMB],trace[FPLIMB];
extern mp_limb_t prime2[FPLIMB2];
extern fp2_t Alpha_1,Alpha_1_inv;
extern mp_limb_t epsilon1[FPLIMB],epsilon2[FPLIMB];
extern mp_limb_t Two_inv[FPLIMB];
extern mpz_t Two_inv_z;
extern mpz_t root_2,root_X;
extern mpz_t efp_total,efp12_total;
extern fp2_t frobenius_constant[12][6];
extern fp2_t skew_frobenius_constant[12][2];
extern mp_limb_t curve_b[FPLIMB];

//montgomery
extern mp_limb_t R[FPLIMB],Ri[FPLIMB],R1[FPLIMB],RR[FPLIMB],Ni[FPLIMB];
extern int m;

extern mp_limb_t u[FPLIMB+1];
extern mp_limb_t N[FPLIMB2],R2[FPLIMB],R3[FPLIMB],RmodP[FPLIMB];
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
extern struct timeval tv_start,tv_end;
extern float MILLER_OPT_AFFINE,FINALEXP_OPT_AFFINE;
extern float MILLER_OPT_PROJECTIVE,FINALEXP_OPT_PROJECTIVE;

extern cost MILLER_OPT_AFFINE_COST,FINALEXP_OPT_AFFINE_COST;
extern cost MILLER_OPT_PROJECTIVE_COST,FINALEXP_OPT_PROJECTIVE_COST;

#endif
