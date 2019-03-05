#ifndef DEFINE_H
#define DEFINE_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

#define FPLIMB_BITS 512
#define FPLIMB 8
#define FPLIMB2 16

#define BN12_X_length 114
#define BN12_X6_2_length 116

#define BLS12_X_length 77
#define BLS12_X2_length 76

#define BLS12_G1_SCM BLS12_EFp12_G1_SCM_2split_JSF_Jacobian_lazy
#define BLS12_G2_SCM BLS12_EFp12_G2_SCM_4split_Jacobian_lazy
#define BLS12_G3_EXP BLS12_Fp12_G3_EXP_4split_lazy
#define BLS12_PAIRING BLS12_Opt_ate_pairing_lazy


int BN12_X_binary[BN12_X_length+1];
int BN12_X6_2_binary[BN12_X6_2_length+1];


int BLS12_X_binary[BLS12_X_length+1];
int BLS12_X2_binary[BLS12_X2_length+1];

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
mp_limb_t buf[FPLIMB],tmp_mul[FPLIMB2],tmp1[FPLIMB],tmp2[FPLIMB];


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
/* Prjective Elliptic Curve                                                   */
/*============================================================================*/
typedef struct{
	Fp x,y,z;
	int infinity;
}EFpZ;

typedef struct{
	Fp2 x,y,z;
	int infinity;
}EFpZ2;

typedef struct{
	Fp6 x,y,z;
	int infinity;
}EFpZ6;

typedef struct{
	Fp12 x,y,z;
	int infinity;
}EFpZ12;


/*============================================================================*/
/* Prjective Elliptic Curve Table                                             */
/*============================================================================*/
typedef struct{
	Fp x,y,z,zz,zzz;
	int infinity;
}EFpZT;

typedef struct{
	Fp2 x,y,z,zz,zzz;
	int infinity;
}EFpZT2;

typedef struct{
	Fp6 x,y,z,zz,zzz;
	int infinity;
}EFpZT6;

typedef struct{
	Fp12 x,y,z,zz,zzz;
	int infinity;
}EFpZT12;

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/
enum f_state{
	f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12
};
gmp_randstate_t state;
mpz_t X_z,prime_z,order_z,trace_z;
mp_limb_t X,prime[FPLIMB],order[FPLIMB],trace[FPLIMB];
Fp2 Alpha_1,Alpha_1_inv;
mp_limb_t epsilon1[FPLIMB],epsilon2[FPLIMB];
mp_limb_t Two_inv[FPLIMB];
mpz_t Two_inv_z;
mpz_t root_2,root_X;
mpz_t EFp_total,EFp12_total;
Fp2 frobenius_constant[12][6];
Fp2 skew_frobenius_constant[12][2];
mp_limb_t curve_b[FPLIMB];

//montgomery
mp_limb_t R[FPLIMB],Ri[FPLIMB],R1[FPLIMB],RR[FPLIMB],Ni[FPLIMB];
int m;
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
struct timeval tv_start,tv_end;
float MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
float FINALEXP_PLAIN,FINALEXP_OPT_EASY,FINALEXP_OPT_HARD;

float G1SCM_PLAIN,G1SCM_2SPLIT,G1SCM_2SPLIT_JSF;
float G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_2SPLIT_JSF,G2SCM_4SPLIT;
float G3EXP_PLAIN,G3EXP_2SPLIT,G3EXP_2SPLIT_JSF,G3EXP_4SPLIT;
#endif
