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

#define X64

#ifdef X64
#define FPLIMB_BITS 512
#define FPLIMB 8
#define FPLIMB2 16
#endif

#ifdef X32
#define FPLIMB_BITS 256
#define FPLIMB 15
#define FPLIMB2 30
#endif


#define BN12_X_length 114
#define BN12_X6_2_length 116

#define BLS12_X_length 77
#define BLS12_X2_length 76

//#define BLS12_X_length 76
//#define BLS12_X2_length 75

#define BLS12_G1_SCM BLS12_EFp12_G1_SCM_2split_JSF_Jacobian_lazy
#define BLS12_G2_SCM BLS12_EFp12_G2_SCM_4split_Jacobian_lazy
#define BLS12_G3_EXP BLS12_Fp12_G3_EXP_4split_GS_lazy
#define BLS12_PAIRING BLS12_Opt_ate_pairing_compress_lazy


extern int BN12_X_binary[BN12_X_length+1];
extern int BN12_X6_2_binary[BN12_X6_2_length+1];


extern int BLS12_X_binary[BLS12_X_length+1];
extern int BLS12_X2_binary[BLS12_X2_length+1];

extern int add,mul;
/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/

typedef struct{
	mp_limb_t x0[FPLIMB];
}Fp;
typedef struct{
	Fp x0;
	Fp x1;
}Fp2;
typedef struct{
    Fp2 x0,x1,x2;
}Fp6;
typedef struct{
    Fp6 x0,x1;
}Fp12;

//tmp finite field
extern mp_limb_t buf[FPLIMB],tmp_mul[FPLIMB2],tmp1[FPLIMB],tmp2[FPLIMB];


/*============================================================================*/
/* Elliptic Curve                                                             */
/*============================================================================*/
typedef struct{
	Fp x,y;
	int infinity;
}EFp;

typedef struct{
    Fp2 x,y;
	int infinity;
}EFp2;

typedef struct{
    Fp6 x,y;
	int infinity;
}EFp6;

typedef struct{
    Fp12 x,y;
	int infinity;
}EFp12;

/*============================================================================*/
/* Jacobian Elliptic Curve                                                   */
/*============================================================================*/
typedef struct{
	Fp x,y,z;
	int infinity;
}EFpP;

typedef struct{
	Fp2 x,y,z;
	int infinity;
}EFpP2;

typedef struct{
	Fp6 x,y,z;
	int infinity;
}EFpP6;

typedef struct{
	Fp12 x,y,z;
	int infinity;
}EFpP12;

/*============================================================================*/
/* Jacobian Elliptic Curve                                                   */
/*============================================================================*/
typedef struct{
	Fp x,y,z;
	int infinity;
}EFpJ;

typedef struct{
	Fp2 x,y,z;
	int infinity;
}EFpJ2;

typedef struct{
	Fp6 x,y,z;
	int infinity;
}EFpJ6;

typedef struct{
	Fp12 x,y,z;
	int infinity;
}EFpJ12;


/*============================================================================*/
/* Jacobian Elliptic Curve Temp                                             */
/*============================================================================*/
typedef struct{
	Fp x,y,z,zz,zzz;
	int infinity;
}EFpJT;

typedef struct{
	Fp2 x,y,z,zz,zzz;
	int infinity;
}EFpJT2;

typedef struct{
	Fp6 x,y,z,zz,zzz;
	int infinity;
}EFpJT6;

typedef struct{
	Fp12 x,y,z,zz,zzz;
	int infinity;
}EFpJT12;

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
extern Fp2 Alpha_1,Alpha_1_inv;
extern mp_limb_t epsilon1[FPLIMB],epsilon2[FPLIMB];
extern mp_limb_t Two_inv[FPLIMB];
extern mpz_t Two_inv_z;
extern mpz_t root_2,root_X;
extern mpz_t EFp_total,EFp12_total;
extern Fp2 frobenius_constant[12][6];
extern Fp2 skew_frobenius_constant[12][2];
extern mp_limb_t curve_b[FPLIMB];

//montgomery
extern mp_limb_t R[FPLIMB],Ri[FPLIMB],R1[FPLIMB],RR[FPLIMB],Ni[FPLIMB];
extern int m;
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
extern struct timeval tv_start,tv_end;
extern float MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
extern float FINALEXP_PLAIN,FINALEXP_OPT,FINALEXP_OPT_EASY,FINALEXP_OPT_HARD;

extern float G1SCM_PLAIN,G1SCM_2SPLIT,G1SCM_2SPLIT_JSF;
extern float G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_2SPLIT_JSF,G2SCM_4SPLIT;
extern float G3EXP_PLAIN,G3EXP_2SPLIT,G3EXP_2SPLIT_JSF,G3EXP_4SPLIT;
#endif
