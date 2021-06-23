#ifndef FILE_H
#define FILE_H
#include "ELiPS/fr.h"

void mpn_write(FILE *file, mp_limb_t *A, mp_size_t size);
void mpn_read(mp_limb_t *A, FILE *file, mp_size_t size);
void fp_write(FILE *file, fp_t *A);
void fp_read(fp_t *A, FILE *file);
void fp2_write(FILE *file, fp2_t *A);
void fp2_read(fp2_t *A, FILE *file);
void fp6_write(FILE *file, fp6_t *A);
void fp6_read(fp6_t *A, FILE *file);
void fp12_write(FILE *file, fp12_t *A);
void fp12_read(fp12_t *A, FILE *file);
void efp_write(FILE *file, efp_t *A);
void efp_read(efp_t *A, FILE *file);
void efp2_write(FILE *file, efp2_t *A);
void efp2_read(efp2_t *A, FILE *file);
void efp6_write(FILE *file, efp6_t *A);
void efp6_read(efp6_t *A, FILE *file);
void efp12_write(FILE *file, efp12_t *A);
void efp12_read(efp12_t *A, FILE *file);
void g1_write(FILE *file, g1_t *A);
void g1_read(g1_t *A, FILE *file);
void g2_write(FILE *file, g2_t *A);
void g2_read(g2_t *A, FILE *file);
void g3_write(FILE *file, g3_t *A);
void g3_read(g3_t *A, FILE *file);
void fr_write(FILE *file, fr_t *A);
void fr_read(fr_t *A, FILE *file);
void fr_write(FILE *file,fr_t *A);
void fr_read(fr_t *A,FILE *file);
#endif