#include <ELiPS/bls_twist.h>

void set_twist_g(){
    fp2cv_t tmp1_fp2cv;
    fpm2_t twist_g;
    while(1){
        fp2cv_set_random(&tmp1_fp2cv,state);
        if(fp2cv_legendre(&tmp1_fp2cv)==-1 && fp2cv_legendre3(&tmp1_fp2cv)==-1){
            break;
        }
    }
    fpm2_set_fp2cv(&twist_g,&tmp1_fp2cv);
    fpm2_sqrt(&twist_g_sqrt,&twist_g);
    fpm2_cbrt(&twist_g_cbrt,&twist_g);
    fpm2_sqrt(&twist_g_sxrt,&twist_g_cbrt);
    fpm2_inv(&twist_g_sqrt_inv,&twist_g_sqrt);
    fpm2_inv(&twist_g_cbrt_inv,&twist_g_cbrt);
    fpm2_inv(&twist_g_sxrt_inv,&twist_g_sxrt);

}
void efpm2_to_efp2cv_for_g2(efp2cv_t *ANS, efpm2_t *p){
    int i;
    fpm2_t x,y;
    fpm2_mul(&x,&p->x,&twist_g_cbrt);
    fpm2_mul(&y,&p->y,&twist_g_sqrt);
    for(i=0;i<DEGREE_EXTENTION_FIELD_2;i++){
        fp_set_neg(&ANS->y.x[i],&y.x[0].x[i]);
        fp_set_neg(&ANS->x.x[i],&x.x[0].x[i]);
    }
    ANS->infinity = p->infinity;
}
void efp2cv_to_efpm2_for_g2(efpm2_t *ANS, efp2cv_t *p){
    fpm2_t x,y;
    fpm2_set_fp2cv(&x,&p->x);
    fpm2_set_fp2cv(&y,&p->y);
    fpm2_mul(&ANS->x,&x,&twist_g_cbrt_inv);
    fpm2_mul(&ANS->y,&y,&twist_g_sqrt_inv);
    ANS->infinity = p->infinity;
}

void efpm2_to_efp2cv_for_g1(efp2cv_t *ANS, efpm2_t *P){
    int i;
    for(i=0;i<DEGREE_EXTENTION_FIELD_2;i++){
        fp_set_neg(&ANS->y.x[i],&P->y.x[0].x[i]);
        fp_set_neg(&ANS->x.x[i],&P->x.x[0].x[i]);
    }
    ANS->infinity = P->infinity;
}