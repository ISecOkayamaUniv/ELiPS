#ifndef INCLUDED_bls_miller_h
#define INCLUDED_bls_miller_h

#include <ELiPS/bls_line_calculate.h>
#include <ELiPS/bls_twist.h>
#include <ELiPS/efpm2.h>

extern void bls_miller_for_optate_basic(fpm2_t *f, efpm2_t *p, efpm2_t *q);
extern void bls_2i3_miller_for_optate(fpm2_t *ANS, efpm2_t *P, efpm2_t *Q);

#endif