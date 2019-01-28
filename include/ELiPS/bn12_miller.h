#ifndef BN12_MILLER_H
#define BN12_MILLER_H

#include <ELiPS/sparse.h>
#include <ELiPS/twist.h>
//miller
void BN12_Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);
void BN12_Miller_algo_for_opt_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);
void BN12_Miller_algo_for_opt_ate_lazy(Fp12 *ANS,EFp12 *Q,EFp12 *P);
void BN12_Miller_algo_for_x_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);

#endif
