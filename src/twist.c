#include <ELiPS/twist.h>
//twist
#ifdef TWIST_PHI_INV
void efp12_to_efp2(efp2_t *ANS, efp12_t *A) {
  fp2_set(&ANS->x, &A->x.x0.x1);
  fp2_set(&ANS->y, &A->y.x1.x1);
  ANS->infinity = A->infinity;
}

void efp2_to_efp12(efp12_t *ANS, efp2_t *A) {
  fp12_set_ui_ui(&ANS->x, 0);
  fp12_set_ui_ui(&ANS->y, 0);
  fp2_set(&ANS->x.x0.x1, &A->x);
  fp2_set(&ANS->y.x1.x1, &A->y);
  ANS->infinity = A->infinity;
}
#endif
#ifdef TWIST_PHI
void efp12_to_efp2(efp2_t *ANS, efp12_t *A) {
  fp2_mul_basis(&ANS->x, &A->x.x0.x2);
  fp2_mul_basis(&ANS->y, &A->y.x1.x1);
  ANS->infinity = A->infinity;
}

void efp2_to_efp12(efp12_t *ANS, efp2_t *A) {
  fp2_t tmp;
  fp2_init(&tmp);
  fp12_set_ui_ui(&ANS->x, 0);
  fp12_set_ui_ui(&ANS->y, 0);
  //*(1-i)/2
  fp_set_ui(&tmp.x0, 1);
  fp_set_ui(&tmp.x1, 0);
  fp_sub_ui(&tmp.x1, &tmp.x1, 1);
  fp2_r1shift(&tmp, &tmp);
  fp2_mul(&ANS->x.x0.x2, &A->x, &tmp);
  fp2_mul(&ANS->y.x1.x1, &A->y, &tmp);
  ANS->infinity = A->infinity;
}
#endif
void efp12_to_efp(efp_t *ANS, efp12_t *A) {
  fp_set(&ANS->x, &A->x.x0.x0.x0);
  fp_set(&ANS->y, &A->y.x0.x0.x0);
  ANS->infinity = A->infinity;
}

void efp_to_efp12(efp12_t *ANS, efp_t *A) {
  fp12_set_ui_ui(&ANS->x, 0);
  fp12_set_ui_ui(&ANS->y, 0);
  fp_set(&ANS->x.x0.x0.x0, &A->x);
  fp_set(&ANS->y.x0.x0.x0, &A->y);
  ANS->infinity = A->infinity;
}
