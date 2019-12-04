#include <ELiPS/twist.h>
//twist
void efp12_to_efp2(efp2_t *ANS,efp12_t *A){
    fp2_set(&ANS->x,&A->x.x0.x1);
    fp2_set(&ANS->y,&A->y.x1.x1);
    ANS->infinity=A->infinity;
}

void efp2_to_efp12(efp12_t *ANS,efp2_t *A){
    fp12_set_ui_ui(&ANS->x,0);
    fp12_set_ui_ui(&ANS->y,0);
    fp2_set(&ANS->x.x0.x1,&A->x);
    fp2_set(&ANS->y.x1.x1,&A->y);
    ANS->infinity=A->infinity;
}

void efp12_to_efp(efp_t *ANS,efp12_t *A){
    fp_set(&ANS->x,&A->x.x0.x0.x0);
    fp_set(&ANS->y,&A->y.x0.x0.x0);
    ANS->infinity=A->infinity;
}

void efp_to_efp12(efp12_t *ANS,efp_t *A){
    fp12_set_ui_ui(&ANS->x,0);
    fp12_set_ui_ui(&ANS->y,0);
    fp_set(&ANS->x.x0.x0.x0,&A->x);
    fp_set(&ANS->y.x0.x0.x0,&A->y);
    ANS->infinity=A->infinity;
}
