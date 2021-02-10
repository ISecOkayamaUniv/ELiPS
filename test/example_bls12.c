#include "ELiPS/bls12.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    //Declarate variable
    scalar_t s1,s2,s12;
    efp12_t P,Q;
    efp12_t s1P,s2P,s1Q,s2Q;
    fp12_t ans1,ans2,ans3;
    
    //Initialize variable
    scalar_init(s1);
    scalar_init(s2);
    scalar_init(s12);
    efp12_init(&P);
    efp12_init(&Q);
    fp12_init(&ans1);
    fp12_init(&ans2);
    fp12_init(&ans3);
    efp12_init(&s1P);
    efp12_init(&s2P);
    efp12_init(&s1Q);
    efp12_init(&s2Q);

    //Initialize and print parameter
    bls12_init();
    bls12_print_parameters();
    
    //Generate random scalar
    scalar_random_order(s1);
    scalar_random_order(s2);
    scalar_mul_order(s12,s1,s2);
    
    //Generate rationalpoit P on G1 and Q on G2
    bls12_generate_g1(&P);
    bls12_generate_g2(&Q);
    
    //Scalar mutiplication
    bls12_g1_scm(&s1P,&P,s1);
    bls12_g1_scm(&s2P,&P,s2);
    bls12_g2_scm(&s1Q,&Q,s1);
    bls12_g2_scm(&s2Q,&Q,s2);
    
    //Calculate pairing
    bls12_optate_pairing_projective(&ans1,&P,&Q);
    bls12_optate_pairing_projective(&ans2,&s1P,&s2Q);
    bls12_optate_pairing_projective(&ans3,&s2P,&s1Q);
    
    //Calculate exponentiation
    bls12_g3_exp(&ans1,&ans1,s12);
    
    //Print answer of pairing
	fp12_println("e(P,Q)^s12=",&ans1);
	fp12_println("e([s1]P,[s2]Q)=",&ans2);
	fp12_println("e([s2]P,[s1]Q)=",&ans3);
	
    //Compare answer
    if(fp12_cmp(&ans1,&ans2)==0)  printf("e(P,Q)^s12=e([s1]P,[s2]Q)\n");
    else    printf("e(P,Q)^s12!=e([s1]P,[s2]Q)\n");
    if(fp12_cmp(&ans1,&ans3)==0)  printf("e(P,Q)^s12=e([s2]P,[s1]Q)\n");
    else    printf("e(P,Q)^s12!=e([s2]P,[s1]Q)\n");
    if(fp12_cmp(&ans2,&ans3)==0)  printf("e([s1]P,[s2]Q)=e([s2]P,[s1]Q)\n\n");
    else    printf("e([s1]P,[s2]Q)!=e([s2]P,[s1]Q)\n\n");
    
    
    //Clear variable
    scalar_clear(s12);
    scalar_clear(s1);
    scalar_clear(s2);
    
    return 0;
}



