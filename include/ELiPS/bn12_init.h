#ifndef BN12_INIT_H
#define BN12_INIT_H

#include <ELiPS/bn12_generate_parameter.h>
#include <ELiPS/bn12_generate_curve.h>
#include <ELiPS/Fp12.h>
#include <ELiPS/EFp12.h>

//init/set/clear
void BN12_init();
void BN12_init_parameters();
void BN12_clear();
void BN12_print_parameters();
#endif
