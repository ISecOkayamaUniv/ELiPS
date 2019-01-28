#ifndef BLS12_TEST_H
#define BLS12_TEST_H

#include <ELiPS/bls12_pairing.h>
#include <ELiPS/G1_SCM.h>
#include <ELiPS/G2_SCM.h>
#include <ELiPS/G3_exp.h>


int BLS12_test_rational_point();
void BLS12_test_plain_ate_pairing();
int BLS12_test_opt_ate_pairing();
void BLS12_test_x_ate_pairing();
int BLS12_test_G1_SCM();
int BLS12_test_G2_SCM();
int BLS12_test_G3_EXP();
void BLS12_compare_pairings();

#endif
