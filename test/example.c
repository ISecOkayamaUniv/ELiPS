#include<stdio.h>
#include<ELiPS/bls12.h>
int main(){
    bls12_init();

    g1_t P,s1P,s2P; // rational point of Efp
    g2_t Q,s1Q,s2Q; // rational point of Efp2 (twisted from Efp12) 
    g3_t F,s1s2F,s2s1F; // paring from P and Q
    fr_t s1,s2,s1s2; // scalar (0~order)

    //initialize
    g1_init(&P);
    g1_init(&s1P);
    g1_init(&s2P);

    g2_init(&Q);
    g2_init(&s1Q);
    g2_init(&s2Q);

    g3_init(&F);
    g3_init(&s1s2F);
    g3_init(&s2s1F);


    fr_init(&s1);
    fr_init(&s2);
    fr_init(&s1s2);

    //Generate random scalar
    fr_set_random(&s1,state);
    fr_set_random(&s2,state);

    //Compute scalar
    fr_mul(&s1s2,&s1,&s2);

    //Generate rationalpoint P on G1 and Q on G2
    g1_set_random(&P,state);
    g2_set_random(&Q,state);

    g1_println("P=",&P);
    g2_println("Q=",&Q);


    //Scalar mutiplication
    g1_scm(&s1P,&P,&s1);
    g1_scm(&s2P,&P,&s2);
    g2_scm(&s1Q,&Q,&s1);
    g2_scm(&s2Q,&Q,&s2);
    
    //Calculate pairing
    g1g2_to_g3_pairing(&F,&P,&Q);
    g1g2_to_g3_pairing(&s1s2F,&s1P,&s2Q);
    g1g2_to_g3_pairing(&s2s1F,&s2P,&s1Q);
    
    //Calculate exponentiation
    g3_exp(&F,&F,&s1s2);


    
    //Print answer of pairing
	fp12_printf_montgomery("e(P,Q)^s12=",&F);
    printf("\n");
	fp12_printf_montgomery("e([s1]P,[s2]Q)=",&s1s2F);
    printf("\n");
	fp12_printf_montgomery("e([s2]P,[s1]Q)=",&s2s1F);
    printf("\n");



    //compare
    if(g3_cmp(&F,&s1s2F)==0 && g3_cmp(&F,&s2s1F)==0) printf("\nbillinear success!\n");
    else printf("error!\n");

    

    

}