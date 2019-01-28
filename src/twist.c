#include <ELiPS/twist.h>
//twist
void EFp12_to_EFp2(EFp2 *ANS,EFp12 *A){
    Fp2_set(&ANS->x,&A->x.x0.x1);
    Fp2_set(&ANS->y,&A->y.x1.x1);
    ANS->infinity=A->infinity;
}

void EFp2_to_EFp12(EFp12 *ANS,EFp2 *A){
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp2_set(&ANS->x.x0.x1,&A->x);
    Fp2_set(&ANS->y.x1.x1,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_to_EFp(EFp *ANS,EFp12 *A){
    Fp_set(&ANS->x,&A->x.x0.x0.x0);
    Fp_set(&ANS->y,&A->y.x0.x0.x0);
    ANS->infinity=A->infinity;
}

void EFp_to_EFp12(EFp12 *ANS,EFp *A){
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp_set(&ANS->x.x0.x0.x0,&A->x);
    Fp_set(&ANS->y.x0.x0.x0,&A->y);
    ANS->infinity=A->infinity;
}
