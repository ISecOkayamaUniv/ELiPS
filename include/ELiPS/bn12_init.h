#ifndef BN12_INIT_H
#define BN12_INIT_H

#include <ELiPS/bn12_generate_parameter.h>
#include <ELiPS/bn12_generate_curve.h>
#include <ELiPS/Fp12.h>
#include <ELiPS/EFp12.h>

//init/set/clear
extern void BN12_init();
extern void BN12_init_parameters();
extern void BN12_clear();
extern void BN12_print_parameters();
#endif
