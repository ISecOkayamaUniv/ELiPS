#ifndef FAST_TEST_H
#define FAST_TEST_H

//#include <ELiPS/time.h>
#include <ELiPS/bls12_pairing.h>
#include <ELiPS/G1_SCM.h>
#include <ELiPS/G2_SCM.h>
#include <ELiPS/G3_exp.h>

extern int Fast_test_BLS12_G1_SCM(int scm);
extern int Fast_test_BLS12_G2_SCM(int scm);
extern int Fast_test_BLS12_G3_EXP(int exp);
extern int Fast_test_BLS12_pairing(int pairing);
#endif
