#ifndef BN12_TEST_H
#define BN12_TEST_H

#include <ELiPS/bn12_pairing.h>
#include <ELiPS/bn12_g1_scm.h>
#include <ELiPS/bn12_g2_scm.h>
#include <ELiPS/bn12_g3_exp.h>

//test
extern int bn12_test_rational_point();
extern void bn12_test_ate_pairing();
extern int bn12_test_opt_ate_pairing();
extern void bn12_test_x_ate_pairing();
extern int bn12_test_g1_scm(int scm);
extern int bn12_test_g2_scm(int scm);
extern int bn12_test_g3_exp(int exp);
#endif
