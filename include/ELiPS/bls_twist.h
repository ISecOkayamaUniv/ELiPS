#ifndef INCLUDED_bls_twist_h
#define INCLUDED_bls_twist_h

#include <ELiPS/efp2cv.h>
#include <ELiPS/efpm2.h>
#include <ELiPS/fpm2.h>

extern void set_twist_g();
extern void efpm2_to_efp2cv_for_g2(efp2cv_t *ANS, efpm2_t *p);
extern void efp2cv_to_efpm2_for_g2(efpm2_t *ANS, efp2cv_t *p);
extern void efpm2_to_efp2cv_for_g1(efp2cv_t *ANS, efpm2_t *P);
#endif