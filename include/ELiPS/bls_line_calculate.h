#ifndef INCLUDED_bls_line_calculate_h
#define INCLUDED_bls_line_calculate_h

#include <ELiPS/efpm2.h>

extern void bls_f_ltp(fpm2_t *f, efpm2_t *p, efpm2_t *q, efpm2_t *t);
extern void bls_f_ltt(fpm2_t *ANS, efpm2_t *q, efpm2_t *t);

#endif