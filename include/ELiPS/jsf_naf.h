#ifndef JSF_H
#define JSF_H

#include <ELiPS/mpn.h>

extern void Joint_sparse_form(int **binary, mpz_t S[2], int *loop_length);
extern int w_naf(int *dw, mpz_t d, int w);
#endif
