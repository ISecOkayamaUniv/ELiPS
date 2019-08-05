#include <ELiPS/Fp12.h>
//Fp12
void Fp12_init(Fp12 *A){
    Fp6_init(&A->x0);
    Fp6_init(&A->x1);
}

void Fp12_printf(char *str,Fp12 *A){
    gmp_printf("%s(",str);
    Fp6_printf("",&A->x0);
    gmp_printf(",");
    Fp6_printf("",&A->x1);
    gmp_printf(")");
}

void Fp12_println(char *str,Fp12 *A){
    gmp_printf("%s(",str);
    Fp6_printf("",&A->x0);
    gmp_printf(",");
    Fp6_printf("",&A->x1);
    gmp_printf(")\n");
}
void Fp12_printf_montgomery(char *str,Fp12 *A){
    gmp_printf("%s(",str);
    Fp6_printf_montgomery("",&A->x0);
    gmp_printf(",");
    Fp6_printf_montgomery("",&A->x1);
    gmp_printf(")");
}
void Fp12_set(Fp12 *ANS,Fp12 *A){
    Fp6_set(&ANS->x0,&A->x0);
    Fp6_set(&ANS->x1,&A->x1);
}
void Fp12_set_ui(Fp12 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x0,UI);
    Fp6_set_ui(&ANS->x1,0);
}
void Fp12_set_ui_ui(Fp12 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x0,UI);
    Fp6_set_ui(&ANS->x1,UI);
}

void Fp12_set_mpn(Fp12 *ANS,mp_limb_t *A){
    Fp6_set_mpn(&ANS->x0,A);
    Fp6_set_ui(&ANS->x1,0);    
}

void Fp12_set_neg(Fp12 *ANS,Fp12 *A){
    Fp6_set_neg(&ANS->x0,&A->x0);
    Fp6_set_neg(&ANS->x1,&A->x1);
}
void Fp12_to_montgomery(Fp12 *ANS,Fp12 *A){
    Fp6_to_montgomery(&ANS->x0,&A->x0);
    Fp6_to_montgomery(&ANS->x1,&A->x1);
}
void Fp12_mod_montgomery(Fp12 *ANS,Fp12 *A){
    Fp6_mod_montgomery(&ANS->x0,&A->x0);
    Fp6_mod_montgomery(&ANS->x1,&A->x1);
}
void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state){
    Fp6_set_random(&ANS->x0,state);
    Fp6_set_random(&ANS->x1,state);
}

void Fp12_mul(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    
    //set
    Fp6_mul(&tmp2_Fp6,&A->x1,&B->x1);//b*d
    Fp6_add(&tmp1_Fp6,&A->x0,&A->x1);//a+b
    Fp6_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp6_mul(&ANS->x1,&tmp1_Fp6,&ANS->x1);//(a+b)(c+d)
    Fp6_mul(&tmp1_Fp6,&A->x0,&B->x0);//a*c
    
    //x0
    Fp6_mul_basis(&ANS->x0,&tmp2_Fp6);//b*d*v
    Fp6_add(&ANS->x0,&ANS->x0,&tmp1_Fp6);//a*c+b*d*v
    
    //x1
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp1_Fp6);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
}
void Fp12_mul_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    //set
    Fp6_mul_lazy(&tmp2_Fp6,&A->x1,&B->x1);//b*d
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);//a+b
    Fp6_add_lazy(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp6_mul_lazy(&ANS->x1,&tmp1_Fp6,&ANS->x1);//(a+b)(c+d)
    Fp6_mul_lazy(&tmp1_Fp6,&A->x0,&B->x0);//a*c
    
    //x0
    Fp6_mul_basis(&ANS->x0,&tmp2_Fp6);//b*d*v
    Fp6_add(&ANS->x0,&ANS->x0,&tmp1_Fp6);//a*c+b*d*v
    
    //x1
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp1_Fp6);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
}
void Fp12_mul_lazy_montgomery(Fp12 *ANS,Fp12 *A,Fp12 *B){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    //set
    Fp6_mul_lazy_montgomery(&tmp2_Fp6,&A->x1,&B->x1);//b*d
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);//a+b
    Fp6_add_lazy(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp6_mul_lazy_montgomery(&ANS->x1,&tmp1_Fp6,&ANS->x1);//(a+b)(c+d)
    Fp6_mul_lazy_montgomery(&tmp1_Fp6,&A->x0,&B->x0);//a*c
    
    //x0
    Fp6_mul_basis(&ANS->x0,&tmp2_Fp6);//b*d*v
    Fp6_add(&ANS->x0,&ANS->x0,&tmp1_Fp6);//a*c+b*d*v
    
    //x1
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp1_Fp6);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
}
void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_mul_ui(&ANS->x0,&A->x0,UI);
    Fp6_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_mul_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B){
    Fp6_mul_mpn(&ANS->x0,&A->x0,B);
    Fp6_mul_mpn(&ANS->x1,&A->x1,B);
}
//complex
void Fp12_sqr(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    Fp6_add(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_mul_basis(&tmp2_Fp6,&A->x1);
    Fp6_add(&tmp2_Fp6,&tmp2_Fp6,&A->x0);
    Fp6_mul(&tmp3_Fp6,&A->x0,&A->x1);
	
    //x0
    Fp6_mul(&ANS->x0,&tmp1_Fp6,&tmp2_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3_Fp6);
    Fp6_mul_basis(&tmp1_Fp6,&tmp3_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp1_Fp6);

    //x1
    Fp6_add(&ANS->x1,&tmp3_Fp6,&tmp3_Fp6);
}

/*
//Karat***NG
void Fp12_sqr(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    Fp6_add(&tmp1_Fp6,&A->x0,&A->x1);
    
    Fp6_mul(&tmp2_Fp6,&A->x0,&A->x1);
    Fp6_add(&tmp2_Fp6,&tmp2_Fp6,&tmp2_Fp6);
    Fp6_mul_basis(&tmp3_Fp6,&tmp2_Fp6);
	
    //x0
    Fp6_sqr(&ANS->x0,&tmp1_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3_Fp6);

    //x1
    Fp6_set(&ANS->x1,&tmp2_Fp6);
}
*/
void Fp12_sqr_lazy(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_mul_basis_lazy(&tmp2_Fp6,&A->x1);
    Fp6_add_lazy(&tmp2_Fp6,&tmp2_Fp6,&A->x0);
    Fp6_mul_lazy(&tmp3_Fp6,&A->x0,&A->x1);
	
    //x0
    Fp6_mul_lazy(&ANS->x0,&tmp1_Fp6,&tmp2_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3_Fp6);
    Fp6_mul_basis(&tmp1_Fp6,&tmp3_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp1_Fp6);
    
    //x1
    Fp6_add(&ANS->x1,&tmp3_Fp6,&tmp3_Fp6);
}
void Fp12_sqr_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6;
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_mul_basis_lazy(&tmp2_Fp6,&A->x1);
    Fp6_add_lazy(&tmp2_Fp6,&tmp2_Fp6,&A->x0);
    Fp6_mul_lazy_montgomery(&tmp3_Fp6,&A->x0,&A->x1);
	
    //x0
    Fp6_mul_lazy_montgomery(&ANS->x0,&tmp1_Fp6,&tmp2_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3_Fp6);
    Fp6_mul_basis(&tmp1_Fp6,&tmp3_Fp6);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp1_Fp6);
    
    //x1
    Fp6_add(&ANS->x1,&tmp3_Fp6,&tmp3_Fp6);
}
void Fp12_sqr_cyclotomic(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    //A=a+b*gamma in G3
    //A^2=(1+2b^2*beta)+((a+b)^2-1-b^2*beta-b^2)
    Fp6_add(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_sqr(&tmp1_Fp6,&tmp1_Fp6);           //(a+b)^2
    Fp6_sqr(&tmp2_Fp6,&A->x1);          //b^2
    Fp6_mul_basis(&ANS->x1,&tmp2_Fp6);  //b^2*beta
    Fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    Fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);   //1+2b^2*beta
    
    Fp6_sub(&ANS->x1,&tmp1_Fp6,&ANS->x1);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
    Fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);   //(a+b)^2-1-b^2*beta-b^2
}
void Fp12_sqr_cyclotomic_lazy(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6;
    Fp6_add_lazy(&tmp1_Fp6,&A->x0,&A->x1);
    Fp6_sqr_lazy(&tmp1_Fp6,&tmp1_Fp6);           //(a+b)^2
    Fp6_sqr_lazy(&tmp2_Fp6,&A->x1);          //b^2
    Fp6_mul_basis(&ANS->x1,&tmp2_Fp6);  //b^2*beta
    Fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    Fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);   //1+2b^2*beta
    
    Fp6_sub(&ANS->x1,&tmp1_Fp6,&ANS->x1);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2_Fp6);
    Fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);   //(a+b)^2-1-b^2*beta-b^2
}
void Fp12_sqr_Karabina(Fp12 *ANS,Fp12 *A){
    Fp12_sqr_compressed(ANS,A);
    Fp12_sqr_recover_g1(ANS,ANS);
    Fp12_sqr_recover_g0(ANS,ANS);
}
void Fp12_sqr_compressed(Fp12 *ANS,Fp12 *A){
    static Fp2 g2,g3,g4,g5;
    static Fp2 T0,T1,T2,T3;
    static Fp2 t0,t1,t2;
    static Fp2 B1;
    
    //set
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    Fp2_sqr(&T0,&g4);
    Fp2_sqr(&T1,&g5);
    
    Fp2_mul_basis(&T2,&T1);
    
    Fp2_add(&T2,&T2,&T0);
    
    Fp2_set(&t2,&T2);
    
    Fp2_add(&t0,&g4,&g5);
    Fp2_sqr(&T2,&t0);
    
    Fp2_add(&T0,&T0,&T1);
    Fp2_sub(&T2,&T2,&T0);
    
    Fp2_set(&t0,&T2);
    Fp2_add(&t1,&g2,&g3);
    Fp2_sqr(&T3,&t1);
    Fp2_sqr(&T2,&g2);
    
    Fp2_mul_basis(&t1,&t0);
    
    Fp2_add(&g2,&g2,&t1);
    Fp2_add(&g2,&g2,&g2);
    
    Fp2_add(&g2,&g2,&t1);
    Fp2_sub(&t1,&t2,&g3);
    Fp2_add(&t1,&t1,&t1);
    
    Fp2_sqr(&T1,&g3);
    
    Fp2_add(&g3,&t1,&t2);
    
    Fp2_mul_basis(&T0,&T1);
    
    Fp2_add(&T0,&T0,&T2);
    
    Fp2_set(&t0,&T0);
    Fp2_sub(&g4,&t0,&g4);//add->sub
    Fp2_add(&g4,&g4,&g4);
    
    Fp2_add(&g4,&g4,&t0);
    
    Fp2_add(&T2,&T2,&T1);
    Fp2_sub(&T3,&T3,&T2);//sub->add
    
    Fp2_set(&t0,&T3);
    Fp2_add(&g5,&g5,&t0);
    Fp2_add(&g5,&g5,&g5);
    
    Fp2_add(&g5,&g5,&t0);
    
    //set
    Fp2_set_ui_ui(&ANS->x0.x0,0);
    Fp2_set_ui_ui(&ANS->x1.x1,0);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
    
}

void Fp12_sqr_compressed_lazy(Fp12 *ANS,Fp12 *A){
    static Fp2 g2,g3,g4,g5;
    static Fp2 T0,T1,T2,T3;
    static Fp2 t0,t1,t2;
    static Fp2 B1;
    
    //set
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    Fp2_sqr_lazy(&T0,&g4);
    Fp2_sqr_lazy(&T1,&g5);
    
    Fp2_mul_basis(&T2,&T1);
    
    Fp2_add(&T2,&T2,&T0);
    
    Fp2_set(&t2,&T2);
    
    Fp2_add(&t0,&g4,&g5);
    Fp2_sqr_lazy(&T2,&t0);
    
    Fp2_add(&T0,&T0,&T1);
    Fp2_sub(&T2,&T2,&T0);
    
    Fp2_set(&t0,&T2);
    Fp2_add(&t1,&g2,&g3);
    Fp2_sqr_lazy(&T3,&t1);
    Fp2_sqr_lazy(&T2,&g2);
    
    Fp2_mul_basis(&t1,&t0);
    
    Fp2_add(&g2,&g2,&t1);
    Fp2_add(&g2,&g2,&g2);
    
    Fp2_add(&g2,&g2,&t1);
    Fp2_sub(&t1,&t2,&g3);
    Fp2_add(&t1,&t1,&t1);
    
    Fp2_sqr_lazy(&T1,&g3);
    
    Fp2_add(&g3,&t1,&t2);
    
    Fp2_mul_basis(&T0,&T1);
    
    Fp2_add(&T0,&T0,&T2);
    
    Fp2_set(&t0,&T0);
    Fp2_sub(&g4,&t0,&g4);//add->sub
    Fp2_add(&g4,&g4,&g4);
    
    Fp2_add(&g4,&g4,&t0);
    
    Fp2_add(&T2,&T2,&T1);
    Fp2_sub(&T3,&T3,&T2);//sub->add
    
    Fp2_set(&t0,&T3);
    Fp2_add(&g5,&g5,&t0);
    Fp2_add(&g5,&g5,&g5);
    
    Fp2_add(&g5,&g5,&t0);
    
    //set
    Fp2_set_ui_ui(&ANS->x0.x0,0);
    Fp2_set_ui_ui(&ANS->x1.x1,0);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}
void Fp12_sqr_compressed_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static Fp2 g2,g3,g4,g5;
    static Fp2 T0,T1,T2,T3;
    static Fp2 t0,t1,t2;
    static Fp2 B1;
    
    //set
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    Fp2_sqr_lazy_montgomery(&T0,&g4);
    Fp2_sqr_lazy_montgomery(&T1,&g5);
    
    Fp2_mul_basis(&T2,&T1);
    
    Fp2_add(&T2,&T2,&T0);
    
    Fp2_set(&t2,&T2);
    
    Fp2_add(&t0,&g4,&g5);
    Fp2_sqr_lazy_montgomery(&T2,&t0);
    
    Fp2_add(&T0,&T0,&T1);
    Fp2_sub(&T2,&T2,&T0);
    
    Fp2_set(&t0,&T2);
    Fp2_add(&t1,&g2,&g3);
    Fp2_sqr_lazy_montgomery(&T3,&t1);
    Fp2_sqr_lazy_montgomery(&T2,&g2);
    
    Fp2_mul_basis(&t1,&t0);
    
    Fp2_add(&g2,&g2,&t1);
    Fp2_add(&g2,&g2,&g2);
    
    Fp2_add(&g2,&g2,&t1);
    Fp2_sub(&t1,&t2,&g3);
    Fp2_add(&t1,&t1,&t1);
    
    Fp2_sqr_lazy_montgomery(&T1,&g3);
    
    Fp2_add(&g3,&t1,&t2);
    
    Fp2_mul_basis(&T0,&T1);
    
    Fp2_add(&T0,&T0,&T2);
    
    Fp2_set(&t0,&T0);
    Fp2_sub(&g4,&t0,&g4);//add->sub
    Fp2_add(&g4,&g4,&g4);
    
    Fp2_add(&g4,&g4,&t0);
    
    Fp2_add(&T2,&T2,&T1);
    Fp2_sub(&T3,&T3,&T2);//sub->add
    
    Fp2_set(&t0,&T3);
    Fp2_add(&g5,&g5,&t0);
    Fp2_add(&g5,&g5,&g5);
    
    Fp2_add(&g5,&g5,&t0);
    
    //set
    Fp2_set_ui_ui(&ANS->x0.x0,0);
    Fp2_set_ui_ui(&ANS->x1.x1,0);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}
void Fp12_sqr_recover_g1(Fp12 *ANS,Fp12 *A){
    static Fp2 g1,g2,g3,g4,g5;
    static Fp2 tmp,f;
    static Fp2 t0,t1;//g1=t0/t1
    static Fp2 C12,C02;
    
    //set
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    //if
    if(Fp2_cmp_zero(&g2)==1){
        Fp2_sqr(&tmp,&g5);
        Fp2_mul_basis(&C12,&tmp);
        Fp2_sqr(&C02,&g4);
        Fp2_set(&t0,&C02);
        Fp2_add(&C02,&C02,&C02);
        Fp2_add(&C02,&C02,&t0);
        Fp2_add(&t0,&C12,&C02);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_add(&t1,&g2,&g2);
        Fp2_add(&t1,&t1,&t1);
    //else
    }else{
        Fp2_mul(&t0,&g4,&g5);
        Fp2_add(&t0,&t0,&t0);
        Fp2_set(&t1,&g3);
    }
    
    Fp2_inv(&t1,&t1);
    Fp2_mul(&g1,&t0,&t1);
    
    //set
    Fp2_set_ui_ui(&ANS->x0.x0,0);
    Fp2_set(&ANS->x1.x1,&g1);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}
void Fp12_sqr_recover_g1_noninv(Fp12 *ANS,Fp12 *A){
    static Fp2 g1,g2,g3,g4,g5;
    static Fp2 tmp,f;
    static Fp2 t0,t1;//g1=t0/t1
    static Fp2 C12,C02;
    
    //set
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    //if
    if(Fp2_cmp_zero(&g2)==1){
        Fp2_sqr(&tmp,&g5);
        Fp2_mul_basis(&C12,&tmp);
        Fp2_sqr(&C02,&g4);
        Fp2_set(&t0,&C02);
        Fp2_add(&C02,&C02,&C02);
        Fp2_add(&C02,&C02,&t0);
        Fp2_add(&t0,&C12,&C02);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_add(&t1,&g2,&g2);
        Fp2_add(&t1,&t1,&t1);
    //else
    }else{
        Fp2_mul(&t0,&g4,&g5);
        Fp2_add(&t0,&t0,&t0);
        Fp2_set(&t1,&g3);
    }
    
    //set
    Fp2_set(&ANS->x0.x0,&t0);
    Fp2_set(&ANS->x1.x1,&t1);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}

void Fp12_sqr_recover_g1_montrick(Fp12 *A,Fp12 *B,Fp12 *C,Fp12 *D){
    static Fp2 a0,a1,b0,b1,c0,c1,d0,d1,ai,bi,ci,di;
    static Fp2 inv,all,buf,ab,bc,cd,da;
    static Fp2 g1,g2,g3,g4,g5;
    static Fp2 tmp,f;
    static Fp2 t0,t1;//g1=t0/t1
    static Fp2 C12,C02;
    
    Fp2_set(&a0,&A->x0.x0);
    Fp2_set(&a1,&A->x1.x1);
    Fp2_set(&b0,&B->x0.x0);
    Fp2_set(&b1,&B->x1.x1);
    Fp2_set(&c0,&C->x0.x0);
    Fp2_set(&c1,&C->x1.x1);
    Fp2_set(&d0,&D->x0.x0);
    Fp2_set(&d1,&D->x1.x1);
    
    //TODO:mul cut back
    Fp2_mul(&ab,&a1,&b1);
    Fp2_mul(&bc,&b1,&c1);
    Fp2_mul(&cd,&c1,&d1);
    Fp2_mul(&da,&d1,&a1);
    
    Fp2_mul(&all,&ab,&cd);
    Fp2_inv(&inv,&all);
    
    Fp2_mul(&buf,&b1,&cd);
    Fp2_mul(&ai,&inv,&buf);
    Fp2_mul(&buf,&c1,&da);
    Fp2_mul(&bi,&inv,&buf);
    Fp2_mul(&buf,&d1,&ab);
    Fp2_mul(&ci,&inv,&buf);
    Fp2_mul(&buf,&a1,&bc);
    Fp2_mul(&di,&inv,&buf);

    
    Fp2_mul(&A->x1.x1,&ai,&a0);
    Fp2_mul(&B->x1.x1,&bi,&b0);
    Fp2_mul(&C->x1.x1,&ci,&c0);
    Fp2_mul(&D->x1.x1,&di,&d0);
}

void Fp12_sqr_recover_g1_lazy(Fp12 *ANS,Fp12 *A){
    static Fp2 g1,g2,g3,g4,g5;
    static Fp2 tmp,f;
    static Fp2 t0,t1;//g1=t0/t1
    static Fp2 C12,C02;
    
    //set
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    //if
    if(Fp2_cmp_zero(&g2)==1){
        Fp2_sqr_lazy(&tmp,&g5);
        Fp2_mul_basis(&C12,&tmp);
        Fp2_sqr_lazy(&C02,&g4);
        Fp2_set(&t0,&C02);
        Fp2_add(&C02,&C02,&C02);
        Fp2_add(&C02,&C02,&t0);
        Fp2_add(&t0,&C12,&C02);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_add(&t1,&g2,&g2);
        Fp2_add(&t1,&t1,&t1);
    //else
    }else{
        Fp2_mul_lazy(&t0,&g4,&g5);
        Fp2_add(&t0,&t0,&t0);
        Fp2_set(&t1,&g3);
    }
    
    Fp2_inv_lazy(&t1,&t1);
    Fp2_mul_lazy(&g1,&t0,&t1);
    
    //set
    Fp2_set_ui_ui(&ANS->x0.x0,0);
    Fp2_set(&ANS->x1.x1,&g1);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}
void Fp12_sqr_recover_g1_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static Fp2 g1,g2,g3,g4,g5;
    static Fp2 tmp,f;
    static Fp2 t0,t1;//g1=t0/t1
    static Fp2 C12,C02;
    
    //set
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    //if
    if(Fp2_cmp_zero(&g2)==1){
        Fp2_sqr_lazy_montgomery(&tmp,&g5);
        Fp2_mul_basis(&C12,&tmp);
        Fp2_sqr_lazy_montgomery(&C02,&g4);
        Fp2_set(&t0,&C02);
        Fp2_add(&C02,&C02,&C02);
        Fp2_add(&C02,&C02,&t0);
        Fp2_add(&t0,&C12,&C02);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_sub(&t0,&t0,&g3);
        Fp2_add(&t1,&g2,&g2);
        Fp2_add(&t1,&t1,&t1);
    //else
    }else{
        Fp2_mul_lazy_montgomery(&t0,&g4,&g5);
        Fp2_add(&t0,&t0,&t0);
        Fp2_set(&t1,&g3);
    }
    
    Fp2_inv_lazy_montgomery(&t1,&t1);
    Fp2_mul_lazy_montgomery(&g1,&t0,&t1);
    
    //set
    Fp2_set_ui_ui(&ANS->x0.x0,0);
    Fp2_set(&ANS->x1.x1,&g1);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}
void Fp12_sqr_recover_g0(Fp12 *ANS,Fp12 *A){
    static Fp2 g0,g1,g2,g3,g4,g5;
    static Fp2 one;
    static Fp2 t0,t1;
    static Fp2 C12,C02;
    
    Fp2_set_ui(&one,1);
    
    //set
    Fp2_set(&g0,&A->x0.x0);
    Fp2_set(&g1,&A->x1.x1);
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    Fp2_sqr(&t0,&g1);
    Fp2_mul(&t1,&g3,&g4);
    Fp2_sub(&t0,&t0,&t1);
    Fp2_add(&t0,&t0,&t0);
    Fp2_sub(&t0,&t0,&t1);
    Fp2_mul(&t1,&g2,&g5);
    Fp2_add(&t0,&t0,&t1);
    Fp2_mul_basis(&g0,&t0);
    Fp2_add(&g0,&g0,&one);
    
    //set
    Fp2_set(&ANS->x0.x0,&g0);
    Fp2_set(&ANS->x1.x1,&g1);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}

void Fp12_sqr_recover_g0_lazy(Fp12 *ANS,Fp12 *A){
    static Fp2 g0,g1,g2,g3,g4,g5;
    static Fp2 one;
    static Fp2 t0,t1;
    static Fp2 C12,C02;
    
    Fp2_set_ui(&one,1);
    
    //set
    Fp2_set(&g0,&A->x0.x0);
    Fp2_set(&g1,&A->x1.x1);
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    Fp2_sqr_lazy(&t0,&g1);
    Fp2_mul_lazy(&t1,&g3,&g4);
    Fp2_sub(&t0,&t0,&t1);
    Fp2_add(&t0,&t0,&t0);
    Fp2_sub(&t0,&t0,&t1);
    Fp2_mul_lazy(&t1,&g2,&g5);
    Fp2_add(&t0,&t0,&t1);
    Fp2_mul_basis(&g0,&t0);
    Fp2_add(&g0,&g0,&one);
    
    //set
    Fp2_set(&ANS->x0.x0,&g0);
    Fp2_set(&ANS->x1.x1,&g1);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}
void Fp12_sqr_recover_g0_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static Fp2 g0,g1,g2,g3,g4,g5;
    static Fp2 one;
    static Fp2 t0,t1;
    static Fp2 C12,C02;
    
    Fp2_set_ui(&one,1);
    
    //set
    Fp2_set(&g0,&A->x0.x0);
    Fp2_set(&g1,&A->x1.x1);
    Fp2_set(&g2,&A->x1.x0);
    Fp2_set(&g3,&A->x0.x2);
    Fp2_set(&g4,&A->x0.x1);
    Fp2_set(&g5,&A->x1.x2);
    
    Fp2_sqr_lazy_montgomery(&t0,&g1);
    Fp2_mul_lazy_montgomery(&t1,&g3,&g4);
    Fp2_sub(&t0,&t0,&t1);
    Fp2_add(&t0,&t0,&t0);
    Fp2_sub(&t0,&t0,&t1);
    Fp2_mul_lazy_montgomery(&t1,&g2,&g5);
    Fp2_add(&t0,&t0,&t1);
    Fp2_mul_basis(&g0,&t0);
    Fp2_add(&g0,&g0,&one);
    
    //set
    Fp2_set(&ANS->x0.x0,&g0);
    Fp2_set(&ANS->x1.x1,&g1);
    Fp2_set(&ANS->x1.x0,&g2);
    Fp2_set(&ANS->x0.x2,&g3);
    Fp2_set(&ANS->x0.x1,&g4);
    Fp2_set(&ANS->x1.x2,&g5);
}
void Fp12_sqr_GS(Fp12 *ANS,Fp12 *A){
    static Fp2 z0,z1,z2,z3,z4,z5;
    static Fp2 t0,t1,t2,t3;
    static Fp2 tmp0,tmp1;
    
    //set
    Fp2_set(&z0,&A->x0.x0);
    Fp2_set(&z1,&A->x1.x1);
    Fp2_set(&z2,&A->x1.x0);
    Fp2_set(&z3,&A->x0.x2);
    Fp2_set(&z4,&A->x0.x1);
    Fp2_set(&z5,&A->x1.x2);
    
    Fp4_sqr(&t0,&t1,&z0,&z1);
    
    Fp2_sub(&z0,&t0,&z0);
    Fp2_add(&z0,&z0,&z0);
    Fp2_add(&z0,&z0,&t0);
    
    Fp2_add(&z1,&t1,&z1);
    Fp2_add(&z1,&z1,&z1);
    Fp2_add(&z1,&z1,&t1);
    
    Fp4_sqr(&t0,&t1,&z2,&z3);
    Fp4_sqr(&t2,&t3,&z4,&z5);
    
    Fp2_sub(&z4,&t0,&z4);
    Fp2_add(&z4,&z4,&z4);
    Fp2_add(&z4,&z4,&t0);
    
    Fp2_add(&z5,&t1,&z5);
    Fp2_add(&z5,&z5,&z5);
    Fp2_add(&z5,&z5,&t1);
    
    Fp2_mul_basis(&t0,&t3);
    Fp2_add(&z2,&t0,&z2);
    Fp2_add(&z2,&z2,&z2);
    Fp2_add(&z2,&z2,&t0);
    
    Fp2_sub(&z3,&t2,&z3);
    Fp2_add(&z3,&z3,&z3);
    Fp2_add(&z3,&z3,&t2);
    //set
    Fp2_set(&ANS->x0.x0,&z0);
    Fp2_set(&ANS->x1.x1,&z1);
    Fp2_set(&ANS->x1.x0,&z2);
    Fp2_set(&ANS->x0.x2,&z3);
    Fp2_set(&ANS->x0.x1,&z4);
    Fp2_set(&ANS->x1.x2,&z5);
}
void Fp4_sqr(Fp2 *t0,Fp2 *t1,Fp2 *g0,Fp2 *g1){
    static Fp2 buf0,buf1;
    
    Fp2_sqr(&buf0,g0);
    Fp2_sqr(&buf1,g1);
    Fp2_add(t1,g0,g1);
    Fp2_mul_basis(t0,&buf1);
    Fp2_add(t0,t0,&buf0);
    Fp2_sqr(t1,t1);
    Fp2_sub(t1,t1,&buf0);
    Fp2_sub(t1,t1,&buf1);
}

void Fp12_sqr_GS_lazy(Fp12 *ANS,Fp12 *A){
    static Fp2 z0,z1,z2,z3,z4,z5;
    static Fp2 t0,t1,t2,t3;
    static Fp2 tmp0,tmp1;
    
    //set
    Fp2_set(&z0,&A->x0.x0);
    Fp2_set(&z1,&A->x1.x1);
    Fp2_set(&z2,&A->x1.x0);
    Fp2_set(&z3,&A->x0.x2);
    Fp2_set(&z4,&A->x0.x1);
    Fp2_set(&z5,&A->x1.x2);
    
    Fp4_sqr_lazy(&t0,&t1,&z0,&z1);
    
    Fp2_sub(&z0,&t0,&z0);
    Fp2_add(&z0,&z0,&z0);
    Fp2_add(&z0,&z0,&t0);
    
    Fp2_add(&z1,&t1,&z1);
    Fp2_add(&z1,&z1,&z1);
    Fp2_add(&z1,&z1,&t1);
    
    Fp4_sqr_lazy(&t0,&t1,&z2,&z3);
    Fp4_sqr_lazy(&t2,&t3,&z4,&z5);
    
    Fp2_sub(&z4,&t0,&z4);
    Fp2_add(&z4,&z4,&z4);
    Fp2_add(&z4,&z4,&t0);
    
    Fp2_add(&z5,&t1,&z5);
    Fp2_add(&z5,&z5,&z5);
    Fp2_add(&z5,&z5,&t1);
    
    Fp2_mul_basis(&t0,&t3);
    Fp2_add(&z2,&t0,&z2);
    Fp2_add(&z2,&z2,&z2);
    Fp2_add(&z2,&z2,&t0);
    
    Fp2_sub(&z3,&t2,&z3);
    Fp2_add(&z3,&z3,&z3);
    Fp2_add(&z3,&z3,&t2);
    //set
    Fp2_set(&ANS->x0.x0,&z0);
    Fp2_set(&ANS->x1.x1,&z1);
    Fp2_set(&ANS->x1.x0,&z2);
    Fp2_set(&ANS->x0.x2,&z3);
    Fp2_set(&ANS->x0.x1,&z4);
    Fp2_set(&ANS->x1.x2,&z5);
}
void Fp12_sqr_GS_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static Fp2 z0,z1,z2,z3,z4,z5;
    static Fp2 t0,t1,t2,t3;
    static Fp2 tmp0,tmp1;
    
    //set
    Fp2_set(&z0,&A->x0.x0);
    Fp2_set(&z1,&A->x1.x1);
    Fp2_set(&z2,&A->x1.x0);
    Fp2_set(&z3,&A->x0.x2);
    Fp2_set(&z4,&A->x0.x1);
    Fp2_set(&z5,&A->x1.x2);
    
    Fp4_sqr_lazy_montgomery(&t0,&t1,&z0,&z1);
    
    Fp2_sub(&z0,&t0,&z0);
    Fp2_add(&z0,&z0,&z0);
    Fp2_add(&z0,&z0,&t0);
    
    Fp2_add(&z1,&t1,&z1);
    Fp2_add(&z1,&z1,&z1);
    Fp2_add(&z1,&z1,&t1);
    
    Fp4_sqr_lazy_montgomery(&t0,&t1,&z2,&z3);
    Fp4_sqr_lazy_montgomery(&t2,&t3,&z4,&z5);
    
    Fp2_sub(&z4,&t0,&z4);
    Fp2_add(&z4,&z4,&z4);
    Fp2_add(&z4,&z4,&t0);
    
    Fp2_add(&z5,&t1,&z5);
    Fp2_add(&z5,&z5,&z5);
    Fp2_add(&z5,&z5,&t1);
    
    Fp2_mul_basis(&t0,&t3);
    Fp2_add(&z2,&t0,&z2);
    Fp2_add(&z2,&z2,&z2);
    Fp2_add(&z2,&z2,&t0);
    
    Fp2_sub(&z3,&t2,&z3);
    Fp2_add(&z3,&z3,&z3);
    Fp2_add(&z3,&z3,&t2);
    //set
    Fp2_set(&ANS->x0.x0,&z0);
    Fp2_set(&ANS->x1.x1,&z1);
    Fp2_set(&ANS->x1.x0,&z2);
    Fp2_set(&ANS->x0.x2,&z3);
    Fp2_set(&ANS->x0.x1,&z4);
    Fp2_set(&ANS->x1.x2,&z5);
}
void Fp4_sqr_lazy(Fp2 *t0,Fp2 *t1,Fp2 *g0,Fp2 *g1){
    static Fp2 buf0,buf1;
    
    Fp2_sqr_lazy(&buf0,g0);
    Fp2_sqr_lazy(&buf1,g1);
    Fp2_add_lazy(t1,g0,g1);
    Fp2_mul_basis(t0,&buf1);
    Fp2_add(t0,t0,&buf0);
    Fp2_sqr_lazy(t1,t1);
    Fp2_sub(t1,t1,&buf0);
    Fp2_sub(t1,t1,&buf1);
}
void Fp4_sqr_lazy_montgomery(Fp2 *t0,Fp2 *t1,Fp2 *g0,Fp2 *g1){
    static Fp2 buf0,buf1;
    
    Fp2_sqr_lazy_montgomery(&buf0,g0);
    Fp2_sqr_lazy_montgomery(&buf1,g1);
    Fp2_add_lazy(t1,g0,g1);
    Fp2_mul_basis(t0,&buf1);
    Fp2_add(t0,t0,&buf0);
    Fp2_sqr_lazy_montgomery(t1,t1);
    Fp2_sub(t1,t1,&buf0);
    Fp2_sub(t1,t1,&buf1);
}
void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_add(&ANS->x0,&A->x0,&B->x0);
    Fp6_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp12_add_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_add_lazy(&ANS->x0,&A->x0,&B->x0);
    Fp6_add_lazy(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_add_ui(&ANS->x0,&A->x0,UI);
    Fp6_add_ui(&ANS->x1,&A->x1,0);
}

void Fp12_add_ui_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_add_ui_ui(&ANS->x0,&A->x0,UI);
    Fp6_add_ui_ui(&ANS->x1,&A->x1,UI);
}
void Fp12_add_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B){
    Fp6_add_mpn(&ANS->x0,&ANS->x0,B);
    Fp6_add_mpn(&ANS->x1,&ANS->x1,B);
}

void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_sub(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub(&ANS->x1,&A->x1,&B->x1);
}
void Fp12_sub_lazy(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_sub_lazy(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub_lazy(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_sub_ui(&ANS->x0,&ANS->x0,UI);
    Fp6_sub_ui(&ANS->x1,&ANS->x1,0);
}

void Fp12_sub_ui_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_sub_ui_ui(&ANS->x0,&ANS->x0,UI);
    Fp6_sub_ui_ui(&ANS->x1,&ANS->x1,UI);
}
void Fp12_sub_mpn(Fp12 *ANS,Fp12 *A,mp_limb_t *B){
    Fp6_sub_mpn(&ANS->x0,&ANS->x0,B);
    Fp6_sub_mpn(&ANS->x1,&ANS->x1,B);
}

void Fp12_inv(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6,tmp4_Fp6;
    Fp6_set(&tmp1_Fp6,&A->x0);
    Fp6_set_neg(&tmp2_Fp6,&A->x1);
    
    Fp6_mul(&tmp3_Fp6,&tmp1_Fp6,&A->x0);
    Fp6_mul(&tmp4_Fp6,&tmp2_Fp6,&A->x1);
    Fp6_mul_basis(&tmp4_Fp6,&tmp4_Fp6);
    Fp6_add(&tmp3_Fp6,&tmp3_Fp6,&tmp4_Fp6);
    Fp6_inv(&tmp3_Fp6,&tmp3_Fp6);
    Fp6_mul(&ANS->x0,&tmp1_Fp6,&tmp3_Fp6);
    Fp6_mul(&ANS->x1,&tmp2_Fp6,&tmp3_Fp6);
}
void Fp12_inv_lazy(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6,tmp4_Fp6;
    Fp6_set(&tmp1_Fp6,&A->x0);
    Fp6_set_neg(&tmp2_Fp6,&A->x1);
    
    Fp6_mul_lazy(&tmp3_Fp6,&tmp1_Fp6,&A->x0);
    Fp6_mul_lazy(&tmp4_Fp6,&tmp2_Fp6,&A->x1);
    Fp6_mul_basis(&tmp4_Fp6,&tmp4_Fp6);
    Fp6_add(&tmp3_Fp6,&tmp3_Fp6,&tmp4_Fp6);
    Fp6_inv_lazy(&tmp3_Fp6,&tmp3_Fp6);
    Fp6_mul_lazy(&ANS->x0,&tmp1_Fp6,&tmp3_Fp6);
    Fp6_mul_lazy(&ANS->x1,&tmp2_Fp6,&tmp3_Fp6);
}
void Fp12_inv_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static Fp6 tmp1_Fp6,tmp2_Fp6,tmp3_Fp6,tmp4_Fp6;
    Fp6_set(&tmp1_Fp6,&A->x0);
    Fp6_set_neg(&tmp2_Fp6,&A->x1);
    
    Fp6_mul_lazy_montgomery(&tmp3_Fp6,&tmp1_Fp6,&A->x0);
    Fp6_mul_lazy_montgomery(&tmp4_Fp6,&tmp2_Fp6,&A->x1);
    Fp6_mul_basis(&tmp4_Fp6,&tmp4_Fp6);
    Fp6_add(&tmp3_Fp6,&tmp3_Fp6,&tmp4_Fp6);
    Fp6_inv_lazy_montgomery(&tmp3_Fp6,&tmp3_Fp6);
    Fp6_mul_lazy_montgomery(&ANS->x0,&tmp1_Fp6,&tmp3_Fp6);
    Fp6_mul_lazy_montgomery(&ANS->x1,&tmp2_Fp6,&tmp3_Fp6);
}
int  Fp12_legendre(Fp12 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp12 tmp;
    Fp12_init(&tmp);
    
    mpz_pow_ui(exp,prime_z,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&tmp,A,exp);
    
    if(Fp12_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

void Fp12_sqrt(Fp12 *ANS,Fp12 *A){
    Fp12 x,y,t,k,n,tmp;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    //gmp_randstate_t state;
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp12_set_random(&n,state);
    while(Fp12_legendre(&n)!=-1){
        Fp12_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp12_pow(&y,&n,q);
    mpz_set_ui(z,e);    
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&x,A,exp);
    Fp12_mul(&tmp,&x,&x);
    Fp12_mul(&k,&tmp,A);
    Fp12_mul(&x,&x,A);
    while(Fp12_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp12_pow(&tmp,&k,exp);
        while(Fp12_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp12_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp12_pow(&t,&y,result);
        Fp12_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp12_mul(&x,&x,&t); 
        Fp12_mul(&k,&k,&y);
    }
    Fp12_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
}

void Fp12_pow(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    Fp12 tmp;
    Fp12_init(&tmp);
    Fp12_set(&tmp,A);
    
    for(i=1;i<length; i++){
        Fp12_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp12_mul(&tmp,A,&tmp);
        }
    }
    
    Fp12_set(ANS,&tmp);
}

int  Fp12_cmp(Fp12 *A,Fp12 *B){
    if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI){
    if(Fp6_cmp_ui(&A->x0,UI)==0 && Fp6_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_mpn(Fp12 *A,mp_limb_t *B){
    if(Fp6_cmp_mpn(&A->x0,B)==0 && Fp6_cmp_mpn(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_zero(Fp12 *A){
    if(Fp6_cmp_zero(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_one(Fp12 *A){
    if(Fp6_cmp_one(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

void Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp2_mul_mpn(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant[f_p1][1].x1.x0);
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    Fp2_mul_mpn(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant[f_p1][2].x0.x0);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p1][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p1][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p1][5]);
}
void Fp12_frobenius_map_p1_lazy(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp2_mul_mpn(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant[f_p1][1].x1.x0);
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    Fp2_mul_mpn(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant[f_p1][2].x0.x0);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul_lazy(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p1][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p1][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul_lazy(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p1][5]);
}
void Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p2][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p2][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p2][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p2][5].x0.x0);
}

void Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p3][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p3][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p3][5]);
}
void Fp12_frobenius_map_p3_lazy(Fp12 *ANS,Fp12 *A){
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul_lazy(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p3][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p3][4].x0.x0);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul_lazy(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p3][5]);
}

void Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p4][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p4][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p4][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p4][5].x0.x0);
}

void Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A){
    //x0
    Fp6_set(&ANS->x0,&A->x0);
    //x1
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p8][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p8][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p8][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p8][5].x0.x0);
}

void Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p10][1].x0.x0);
    Fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p10][2].x0.x0);
    //x1
    Fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p10][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p10][5].x0.x0);
}
void Fp12_frobenius_map_p1_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static mp_limb_t buf[FPLIMB];
    static Fp2 buf2;
    //TODO:global
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    mpn_to_montgomery(buf,frobenius_constant[f_p1][1].x1.x0);
    Fp2_mul_mpn_montgomery(&ANS->x0.x1,&ANS->x0.x1,buf);
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    mpn_to_montgomery(buf,frobenius_constant[f_p1][2].x0.x0);
    Fp2_mul_mpn_montgomery(&ANS->x0.x2,&ANS->x0.x2,buf);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_to_montgomery(&buf2,&frobenius_constant[f_p1][3]);
    Fp2_mul_lazy_montgomery(&ANS->x1.x0,&ANS->x1.x0,&buf2);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    mpn_to_montgomery(buf,frobenius_constant[f_p1][4].x0.x0);
    Fp2_mul_mpn_montgomery(&ANS->x1.x1,&ANS->x1.x1,buf);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_to_montgomery(&buf2,&frobenius_constant[f_p1][5]);
    Fp2_mul_lazy_montgomery(&ANS->x1.x2,&ANS->x1.x2,&buf2);
}
void Fp12_frobenius_map_p2_montgomery(Fp12 *ANS,Fp12 *A){
    static mp_limb_t buf[FPLIMB];
    //TODO:global

    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    mpn_to_montgomery(buf,frobenius_constant[f_p2][1].x0.x0);
    Fp2_mul_mpn_montgomery(&ANS->x0.x1,&A->x0.x1,buf);
    mpn_to_montgomery(buf,frobenius_constant[f_p2][2].x0.x0);
    Fp2_mul_mpn_montgomery(&ANS->x0.x2,&A->x0.x2,buf);
    //x1
    mpn_to_montgomery(buf,frobenius_constant[f_p2][3].x0.x0);
    Fp2_mul_mpn_montgomery(&ANS->x1.x0,&A->x1.x0,buf);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    mpn_to_montgomery(buf,frobenius_constant[f_p2][5].x0.x0);
    Fp2_mul_mpn_montgomery(&ANS->x1.x2,&A->x1.x2,buf);
}
void Fp12_frobenius_map_p3_lazy_montgomery(Fp12 *ANS,Fp12 *A){
    static mp_limb_t buf[FPLIMB];
    static Fp2 buf2;
    //TODO:global
    static Fp tmp1_Fp;
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp1_Fp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp1_Fp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_to_montgomery(&buf2,&frobenius_constant[f_p3][3]);
    Fp2_mul_lazy_montgomery(&ANS->x1.x0,&ANS->x1.x0,&buf2);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    mpn_to_montgomery(buf,frobenius_constant[f_p3][4].x0.x0);
    Fp2_mul_mpn_montgomery(&ANS->x1.x1,&ANS->x1.x1,buf);
    Fp_add(&tmp1_Fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp1_Fp);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_to_montgomery(&buf2,&frobenius_constant[f_p3][5]);
    Fp2_mul_lazy_montgomery(&ANS->x1.x2,&ANS->x1.x2,&buf2);
}

void Fp12_frobenius_map_p6_montgomery(Fp12 *ANS,Fp12 *A){
    //x0
    Fp6_set(&ANS->x0,&A->x0);
    //x1
    Fp6_set_neg(&ANS->x1,&A->x1);
}

int Fp12_montgomery_trick(Fp12 *A_inv,Fp12 *A,int n){
    int i;
    Fp12 ANS[n],ALL_inv;
	Fp12_set(&ANS[0],&A[0]);
	Fp12 check;
	
	for(i=1;i<n;i++){
	Fp12_mul_lazy(&ANS[i],&ANS[i-1],&A[i]);
	}
	Fp12_inv_lazy(&ALL_inv,&ANS[n-1]);	
	for(i=n-1;i>0;i--){
    Fp12_mul_lazy(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    Fp12_mul_lazy(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    Fp12_set(&A_inv[0],&ALL_inv);
    /*
    for(i=0;i<n;i++){
    Fp12_mul(&check,&A[i],&A_inv[i]);
    printf("check:%d",i);	
	Fp12_println("=",&check);
    }
    */
    return 0;
}
int Fp12_montgomery_trick_montgomery(Fp12 *A_inv,Fp12 *A,int n){
    int i;
    Fp12 ANS[n],ALL_inv;
	Fp12_set(&ANS[0],&A[0]);
	Fp12 check;
	
	for(i=1;i<n;i++){
	Fp12_mul_lazy_montgomery(&ANS[i],&ANS[i-1],&A[i]);
	}
	Fp12_inv_lazy_montgomery(&ALL_inv,&ANS[n-1]);	
	for(i=n-1;i>0;i--){
    Fp12_mul_lazy_montgomery(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    Fp12_mul_lazy_montgomery(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    Fp12_set(&A_inv[0],&ALL_inv);
    return 0;
}
