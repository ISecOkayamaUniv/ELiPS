#ifndef FAST_TEST_H
#define FAST_TEST_H

#include <ELiPS/time.h>
#include <ELiPS/Fp12.h>
#include <ELiPS/EFp12.h>
#include <ELiPS/twist.h>
#include <ELiPS/bn12_pairing.h>
#include <ELiPS/bls12_pairing.h>
#include <ELiPS/G1_SCM.h>
#include <ELiPS/G2_SCM.h>
#include <ELiPS/G3_exp.h>

int Fast_test_BLS12_G1_SCM(int scm);
int Fast_test_BLS12_G2_SCM(int scm);
int Fast_test_BLS12_G3_EXP(int exp);
int Fast_test_BLS12_pairing(int pairing);
#endif
