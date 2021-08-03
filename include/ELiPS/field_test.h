#ifndef FIELD_TEST_H
#define FIELD_TEST_H

#include <ELiPS/bls12_test.h>
#include <ELiPS/time.h>

extern int test_fp(int fp_n);
extern int test_fp2(int fp2_n);
extern int test_fp6(int fp6_n);
extern int test_field(int fp, int fp2, int fp6, int fp12, int sqr);
extern int test_fp_montgomery(int fp_n);
#endif
