#ifndef INCLUDED_bls_line_calculate_h
#define INCLUDED_bls_line_calculate_h

#include <ELiPS/efpm2.h>
#include <ELiPS/efp2cv.h>

extern void bls_f_ltp(fpm2_t *f, efpm2_t *p, efpm2_t *q, efpm2_t *t);
extern void bls_f_ltt(fpm2_t *ANS, efpm2_t *q, efpm2_t *t);
extern void bls_2i3_ff_ltt(fpm2_t *f, efp2cv_t *P, efp2cv_t *T);
extern void bls_2i3_ff_ltq(fpm2_t *f, efp2cv_t *P, efp2cv_t *Q, efp2cv_t *T);

#endif