#ifndef LAZY_TEST_H
#define LAZY_TEST_H


#include <ELiPS/bn12_pairing.h>
#include <ELiPS/bls12_pairing.h>
#include <ELiPS/G1_SCM.h>
#include <ELiPS/G2_SCM.h>
#include <ELiPS/G3_exp.h>
//#include <ELiPS/time.h>

/**
 * @brief test field (Lazy Reduction)
 *
 * @param[in]fp2 --a number of trials on Fp2 test.
 * @param[in]fp6 --a number of trials on Fp6 test.
 * @param[in]fp12 --a number of trials on Fp12 test.
 * 
 * @return int --(All success 0 or failed 1)
 */
extern int test_Field_Lazy(int fp2,int fp6,int fp12);
extern int test_EFp_Lazy(int ecd,int eca,int scm);
extern int test_EFp2_Lazy(int ecd,int eca,int scm);
extern int test_EFp12_Lazy(int ecd,int eca,int scm);
extern void test_BN12_G1_SCM_Lazy();
extern int test_BLS12_G1_SCM_Lazy(int scm);
extern void test_BN12_G2_SCM_Lazy();
extern int test_BLS12_G2_SCM_Lazy(int scm);
extern void test_BN12_G3_exp_Lazy();
extern int test_BLS12_G3_exp_Lazy(int exp);
extern void test_BN12_opt_ate_pairing_Lazy();
extern int test_BLS12_opt_ate_pairing_Lazy(int pairing);

#endif
