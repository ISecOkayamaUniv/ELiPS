#ifndef MACRO_H
#define MACRO_H

#define FPLIMB 8
#define FPLIMB2 16

/***Calculation type***/
#define Optimaization_type
#define Basic_type
#define Custom_type

/***Field***/
#define Field_basic
//#define Field_lazy

/***Eliptiic Curve***/
#define EC_basic
//#define EC_lazy

#define EC_Affine
//#define EC_Jacobian
//#define EC_Mixture

/***Scalar***/
#define scalar_basic
//#define scalar_JSF
//#define scalar_naf

/***GLV***/
#define GLV_1
//#define GLV_2
//#define GLV_4

/***Reduction***/
#define Reduction_basic
#define Reduction_montgomery

/***Field Macro***/
#ifdef Field_basic
#define Fp_add Fp_add_basic
#define Fp_sub Fp_sub_basic
#define Fp_mul Fp_mul_basic
#define Fp_sqr Fp_sqr_basic

#define Fp2_add Fp2_add_basic
#define Fp2_sub Fp2_sub_basic
#define Fp2_mul Fp2_mul_basic
#define Fp2_sqr Fp2_sqr_basic

#define Fp6_add Fp6_add_basic
#define Fp6_sub Fp6_sub_basic
#define Fp6_mul Fp6_mul_basic
#define Fp6_sqr Fp6_sqr_basic

#define Fp12_add Fp12_add_basic
#define Fp12_sub Fp12_sub_basic
#define Fp12_mul Fp12_mul_basic
#define Fp12_sqr Fp12_sqr_basic
#endif

#ifdef Field_lazy
#define Fp_add Fp_add_basic
#define Fp_sub Fp_sub_basic
#define Fp_mul Fp_mul_basic
#define Fp_sqr Fp_sqr_basic

#define Fp2_add Fp2_add_lazy
#define Fp2_sub Fp2_sub_lazy
#define Fp2_mul Fp2_mul_lazy
#define Fp2_sqr Fp2_sqr_lazy

#define Fp6_add Fp6_add_lazy
#define Fp6_sub Fp6_sub_lazy
#define Fp6_mul Fp6_mul_lazy
#define Fp6_sqr Fp6_sqr_lazy

#define Fp12_add Fp12_add_lazy
#define Fp12_sub Fp12_sub_lazy
#define Fp12_mul Fp12_mul_lazy
#define Fp12_sqr Fp12_sqr_lazy
#endif

#define BLS12_EFp12_G1_SCM BLS12_EFp12_G1_SCM_2split_7NAF_interleaving_Mixture_lazy
#define BLS12_G2_SCM BLS12_EFp12_G2_SCM_4split_5NAF_interleaving_Mixture_lazy_montgomery
#define BLS12_G3_EXP BLS12_Fp12_G3_EXP_4split_5NAF_interleaving_GS_lazy_montgomery
#define BLS12_Opt_ate_pairing BLS12_Opt_ate_pairing_compress_lazy_montgomery

#endif
