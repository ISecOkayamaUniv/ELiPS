#ifndef BLS12_TEST_H
#define BLS12_TEST_H

#include <ELiPS/bls12_pairing.h>
#include <ELiPS/bls12_g1_scm.h>
#include <ELiPS/bls12_g2_scm.h>
#include <ELiPS/bls12_g3_exp.h>


extern int bls12_test_rational_point();
extern void bls12_test_plain_ate_pairing();
extern int bls12_test_opt_ate_pairing();
extern int bls12_test_g1_scm(int scm);
extern int bls12_test_g2_scm(int scm);
extern int bls12_test_g3_exp(int exp);
#endif
