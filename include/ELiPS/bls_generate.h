#ifndef INCLUDED_bls_generate_h
#define INCLUDED_bls_generate_h

#include <ELiPS/define_cvma.h>
#include <ELiPS/efp.h>

extern void bls_init();
extern void bls_clear();
extern void bls_build();
extern void bls_print_parameter();
extern void bls_generate_X();
extern void bls_generate_parameter_of_bls12();
extern void bls_search_parameter_all();
extern void bls_weil_of_bls12();
extern void bls_weil_of_bls48();
extern void bls_weil();
extern int combination(int n, int r);
extern void bls_search_curve();
#endif