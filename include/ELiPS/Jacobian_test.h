#ifndef JACOBIAN_TEST_H
#define JACOBIAN_TEST_H


#include <ELiPS/bn12_pairing.h>
#include <ELiPS/bls12_pairing.h>
#include <ELiPS/G1_SCM.h>
#include <ELiPS/G2_SCM.h>
#include <ELiPS/G3_exp.h>
//#include <ELiPS/time.h>


extern int test_EFp_Jacobian(int ecd,int eca,int scm);
extern int test_EFp2_Jacobian(int ecd,int eca,int scm);
extern int test_BLS12_G1_SCM_Jacobian(int scm);
extern int test_BLS12_G2_SCM_Jacobian(int scm);

#endif
