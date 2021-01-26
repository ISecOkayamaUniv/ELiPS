#ifndef INCLUDED_efpm_h
#define INCLUDED_efpm_h

#include <ELiPS/fpm.h>

extern void efpm_init(efpm_t *p);
extern void efpm_set(efpm_t *q, efpm_t *p);
extern void efpm_printf(char *str, efpm_t *p);
extern void efpm_rational_point(efpm_t *p);
extern void efpm_eca(efpm_t *p3, efpm_t *p1, efpm_t *p2);
extern void efpm_ecd(efpm_t *p3, efpm_t *p1);
extern void efpm_scm(efpm_t *q, efpm_t *p, mpz_t s);

#endif
