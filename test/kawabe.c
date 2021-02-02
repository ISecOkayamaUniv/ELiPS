#include "ELiPS/bls12.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
void mul_montgomery_test(){
    cost tmp;
    cost fp12_mul_cost,fp4_mul_cost;
    cost fp12_inv_cost,fp4_inv_cost;
    cost_init(&tmp);
    cost_init(&fp12_mul_cost);
    cost_init(&fp12_inv_cost);

    fp12_t A,B,C;
    fp12_init(&A);
    fp12_init(&B);
    fp12_init(&C);

    fp12_set_random(&A,state);
    fp12_to_montgomery(&A,&A);

    fp12_set_random(&B,state);
    fp12_to_montgomery(&B,&B);
    
    cost_zero();
    fp12_mul_lazy_montgomery(&C,&A,&B);
    cost_check(&tmp);
    cost_addition(&fp12_mul_cost,&tmp);

    cost_zero();
    fp12_inv_lazy_montgomery(&C,&A);
    cost_check(&tmp);
    cost_addition(&fp12_inv_cost,&tmp);

    #ifdef DEBUG_COST_A
    printf("*********bls12 g2 scm fp COST.********         \n");
    cost_printf("fp12_mul_lazy_montgomery",&fp12_mul_cost,1);
    printf("***************************************         \n");

    printf("*********bls12 g2 scm fp COST.********         \n");
    cost_printf("fp12_inv_lazy_montgomery",&fp12_inv_cost,1);
    printf("***************************************         \n");
    #endif
}
void mul_test(){
    cost tmp;
    cost fp12_mul_cost,fp4_mul_cost;
    cost fp12_inv_cost,fp4_inv_cost;
    cost_init(&tmp);
    cost_init(&fp12_mul_cost);
    cost_init(&fp12_inv_cost);

    fp12_t A,B,C;
    fp12_init(&A);
    fp12_init(&B);
    fp12_init(&C);

    fp12_set_random(&A,state);

    fp12_set_random(&B,state);
    
    cost_zero();
    fp12_mul(&C,&A,&B);
    cost_check(&tmp);
    cost_addition(&fp12_mul_cost,&tmp);

    cost_zero();
    fp12_inv(&C,&A);
    cost_check(&tmp);
    cost_addition(&fp12_inv_cost,&tmp);

    #ifdef DEBUG_COST_A
    cost_printf("fp12_mul",&fp12_mul_cost,1);
    printf("***************************************         \n");

    cost_printf("fp12_inv",&fp12_inv_cost,1);
    printf("***************************************         \n");
    #endif
}
int main(void){
    bls12_init();
    mul_test();
    mul_montgomery_test();
    
    return 0;
}
