#ifndef BN12_TEST_H
#define BN12_TEST_H

#include <ELiPS/bn12_pairing.h>
#include <ELiPS/G1_SCM.h>
#include <ELiPS/G2_SCM.h>
#include <ELiPS/G3_exp.h>

//test
void BN12_test_rational_point();
void BN12_test_plain_ate_pairing();
void BN12_test_opt_ate_pairing();
void BN12_test_x_ate_pairing();
void BN12_test_G1_SCM();
void BN12_test_G2_SCM();
void BN12_test_G3_EXP();
void BN12_compare_pairings();

#endif
