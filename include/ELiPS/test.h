#ifndef TEST_H
#define TEST_H

#include <ELiPS/time.h>
#include <ELiPS/bls12_test.h>

extern int test_efp(int ecd,int eca);
extern int test_efp2(int ecd,int eca);
extern int test_efp12(int ecd,int eca,int scm);
extern void test_Frobenius_map();
extern void test_skew_frobenius_map();
extern void test_twist();
//extern void test_fp(int fp);
extern int test_mod(int mod);
extern void test_All();
extern int Fast_test_BLS12(int g1, int g2,int g3,int pairing);
#endif
