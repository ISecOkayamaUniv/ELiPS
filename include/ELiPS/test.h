#ifndef TEST_H
#define TEST_H

#include <ELiPS/time.h>
#include <ELiPS/Fp12.h>
#include <ELiPS/EFp12.h>
#include <ELiPS/twist.h>
#include <ELiPS/bls12_test.h>
#include <ELiPS/lazy_test.h>
#include <ELiPS/Jacobian_test.h>


void test_Field();
void test_Frobenius_map();
void test_skew_frobenius_map();
void test_twist();
void test_All();
#endif
