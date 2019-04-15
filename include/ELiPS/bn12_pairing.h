#ifndef BN12_PAIRING_H
#define BN12_PAIRING_H

#include <ELiPS/bn12_miller.h>
#include <ELiPS/bn12_final_exp.h>
#include <ELiPS/test.h>

//pairing
extern void BN12_Plain_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
extern void BN12_Opt_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);
extern void BN12_Opt_ate_pairing_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q);
extern void BN12_X_ate_pairing(Fp12 *ANS,EFp12 *P,EFp12 *Q);

#endif
