#ifndef MACRO_H
#define MACRO_H

#define FPLIMB 8
#define FPLIMB2 16


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
#define scalar_JSF
//#define scalar_naf

/***GLV***/
#define GLV_1
#define GLV_2
#define GLV_4




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


#endif
