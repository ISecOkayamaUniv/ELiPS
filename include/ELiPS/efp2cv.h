#ifndef INCLUDED_efp2cv_h
#define INCLUDED_efp2cv_h

#include <ELiPS/fp2cv.h>

extern void efp2cv_init(efp2cv_t *p);
extern void efp2cv_set(efp2cv_t *q, efp2cv_t *p);
extern void efp2cv_printf(char *str, efp2cv_t *p);
extern void efp2cv_println(char *str, efp2cv_t *p);
extern void efp2cv_rational_point(efp2cv_t *p);
extern void efp2cv_eca(efp2cv_t *p3, efp2cv_t *p1, efp2cv_t *p2);
extern void efp2cv_ecd(efp2cv_t *p3, efp2cv_t *p1);
extern void efp2cv_scm(efp2cv_t *q, efp2cv_t *p, mpz_t s);

#endif