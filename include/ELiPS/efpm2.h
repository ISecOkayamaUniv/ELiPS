#ifndef INCLUDED_efpm2_h
#define INCLUDED_efpm2_h

#include <ELiPS/fpm2.h>

extern void efpm2_init(efpm2_t *p);
extern void efpm2_set(efpm2_t *q, efpm2_t *p);
extern void efpm2_printf(char *str, efpm2_t *p);
extern void efpm2_println(char *str, efpm2_t *p);
extern void efpm2_rational_point(efpm2_t *p);
extern void efpm2_eca(efpm2_t *p3, efpm2_t *p1, efpm2_t *p2);
extern void efpm2_ecd(efpm2_t *p3, efpm2_t *p1);
extern void efpm2_scm(efpm2_t *q, efpm2_t *p, mpz_t s);

#endif
