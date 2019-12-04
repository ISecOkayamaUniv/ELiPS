#ifndef BN12_MILLER_H
#define BN12_MILLER_H

#include <ELiPS/line_calculate.h>
#include <ELiPS/twist.h>
//miller
extern void bn12_miller_algo_for_plain_ate(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bn12_miller_algo_for_opt_ate(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bn12_miller_algo_for_opt_ate_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bn12_miller_algo_for_x_ate(fp12_t *ANS,efp12_t *P,efp12_t *Q);

#endif
