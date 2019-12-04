#ifndef BN12_PAIRING_H
#define BN12_PAIRING_H

#include <ELiPS/bn12_miller.h>
#include <ELiPS/bn12_final_exp.h>
#include <ELiPS/test.h>

//pairing
extern void bn12_ate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bn12_optate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bn12_optate_pairing_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q);
extern void bn12_X_ate_pairing(fp12_t *ANS,efp12_t *P,efp12_t *Q);

#endif
