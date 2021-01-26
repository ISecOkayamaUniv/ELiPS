#ifndef BN12_INIT_H
#define BN12_INIT_H

#include <ELiPS/bn12_generate_parameter.h>
#include <ELiPS/bn12_generate_curve.h>
#include <ELiPS/fp12.h>
#include <ELiPS/efp12.h>

//init/set/clear
extern void bn12_init();
extern void bn12_init_parameters();
extern void bn12_clear();
extern void bn12_print_parameters();
#endif
