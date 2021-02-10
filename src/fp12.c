#include <ELiPS/fp12.h>

void fp12_init(fp12_t *A){
    fp6_init(&A->x0);
    fp6_init(&A->x1);
}

void fp12_printf(char *str,fp12_t *A){
    gmp_printf("%s(",str);
    fp6_printf("",&A->x0);
    gmp_printf(",");
    fp6_printf("",&A->x1);
    gmp_printf(")");
}

void fp12_println(char *str,fp12_t *A){
    gmp_printf("%s(",str);
    fp6_printf("",&A->x0);
    gmp_printf(",");
    fp6_printf("",&A->x1);
    gmp_printf(")\n");
}

void fp12_printf_montgomery(char *str,fp12_t *A){
    gmp_printf("%s(",str);
    fp6_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp6_printf_montgomery("",&A->x1);
    gmp_printf(")");
}

void fp12_println_montgomery(char *str,fp12_t *A){
    gmp_printf("%s(",str);
    fp6_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp6_printf_montgomery("",&A->x1);
    gmp_printf(")\n");
}

void fp12_set(fp12_t *ANS,fp12_t *A){
    fp6_set(&ANS->x0,&A->x0);
    fp6_set(&ANS->x1,&A->x1);
}

void fp12_set_ui(fp12_t *ANS,unsigned long int UI){
    fp6_set_ui(&ANS->x0,UI);
    fp6_set_ui(&ANS->x1,0);
}

void fp12_set_ui_ui(fp12_t *ANS,unsigned long int UI){
    fp6_set_ui_ui(&ANS->x0,UI);
    fp6_set_ui_ui(&ANS->x1,UI);
}

void fp12_set_mpn(fp12_t *ANS,mp_limb_t *A){
    fp6_set_mpn(&ANS->x0,A);
    fp6_set_ui(&ANS->x1,0);
}

void fp12_set_neg(fp12_t *ANS,fp12_t *A){
    fp6_set_neg(&ANS->x0,&A->x0);
    fp6_set_neg(&ANS->x1,&A->x1);
}

void fp12_to_montgomery(fp12_t *ANS,fp12_t *A){
    fp6_to_montgomery(&ANS->x0,&A->x0);
    fp6_to_montgomery(&ANS->x1,&A->x1);
}

void fp12_mod_montgomery(fp12_t *ANS,fp12_t *A){
    fp6_mod_montgomery(&ANS->x0,&A->x0);
    fp6_mod_montgomery(&ANS->x1,&A->x1);
}

void fp12_set_random(fp12_t *ANS,gmp_randstate_t state){
    fp6_set_random(&ANS->x0,state);
    fp6_set_random(&ANS->x1,state);
}

void fp12_mul(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp6_t tmp1_fp6,tmp2_fp6;

    //set
    fp6_mul(&tmp2_fp6,&A->x1,&B->x1);//b*d
    fp6_add(&tmp1_fp6,&A->x0,&A->x1);//a+b
    fp6_add(&ANS->x1,&B->x0,&B->x1);//c+d
    fp6_mul(&ANS->x1,&tmp1_fp6,&ANS->x1);//(a+b)(c+d)
    fp6_mul(&tmp1_fp6,&A->x0,&B->x0);//a*c

    //x0
    fp6_mul_basis(&ANS->x0,&tmp2_fp6);//b*d*v
    fp6_add(&ANS->x0,&ANS->x0,&tmp1_fp6);//a*c+b*d*v

    //x1
    fp6_sub(&ANS->x1,&ANS->x1,&tmp1_fp6);
    fp6_sub(&ANS->x1,&ANS->x1,&tmp2_fp6);
}

void fp12_mul_lazy(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp6_t tmp1_fp6,tmp2_fp6;
    //set
    fp6_mul_lazy(&tmp2_fp6,&A->x1,&B->x1);//b*d
    fp6_add_nonmod_single(&tmp1_fp6,&A->x0,&A->x1);//a+b
    fp6_add_nonmod_single(&ANS->x1,&B->x0,&B->x1);//c+d
    fp6_mul_lazy(&ANS->x1,&tmp1_fp6,&ANS->x1);//(a+b)(c+d)
    fp6_mul_lazy(&tmp1_fp6,&A->x0,&B->x0);//a*c

    //x0
    fp6_mul_basis(&ANS->x0,&tmp2_fp6);//b*d*v
    fp6_add(&ANS->x0,&ANS->x0,&tmp1_fp6);//a*c+b*d*v

    //x1
    fp6_sub(&ANS->x1,&ANS->x1,&tmp1_fp6);
    fp6_sub(&ANS->x1,&ANS->x1,&tmp2_fp6);
}

void fp12_mul_lazy_montgomery(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp6_t tmp1_fp6,tmp2_fp6;
    //set
    fp6_mul_lazy_montgomery(&tmp2_fp6,&A->x1,&B->x1);//b*d
    fp6_add_nonmod_single(&tmp1_fp6,&A->x0,&A->x1);//a+b
    fp6_add_nonmod_single(&ANS->x1,&B->x0,&B->x1);//c+d
    fp6_mul_lazy_montgomery(&ANS->x1,&tmp1_fp6,&ANS->x1);//(a+b)(c+d)
    fp6_mul_lazy_montgomery(&tmp1_fp6,&A->x0,&B->x0);//a*c

    //x0
    fp6_mul_basis(&ANS->x0,&tmp2_fp6);//b*d*v
    fp6_add(&ANS->x0,&ANS->x0,&tmp1_fp6);//a*c+b*d*v

    //x1
    fp6_sub(&ANS->x1,&ANS->x1,&tmp1_fp6);
    fp6_sub(&ANS->x1,&ANS->x1,&tmp2_fp6);
}

//TODO: bug but fast

// void fp12_mul_lazy_montgomery(fp12_t *ANS,fp12_t *A,fp12_t *B){
//     static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
//     static fpd6_t tmp1_fpd6,tmp2_fpd6,tmp3_fpd6,buf_fpd6;
//     static fp6_t out;
//     //set
//     fp6_mul_nonmod_montgomery(&tmp1_fpd6,&A->x1,&B->x1);//b*d
//     fp6_add_nonmod_single(&tmp2_fp6,&A->x0,&A->x1);//a+b
//     fp6_add_nonmod_single(&tmp3_fp6,&B->x0,&B->x1);//c+d
//     fp6_mul_nonmod_montgomery(&tmp2_fpd6,&tmp2_fp6,&tmp3_fp6);//(a+b)(c+d)
//         fp6_mod_montgomery_double(&out,&tmp2_fpd6);
//         //fp6_println("out_A_0=",&out);printf("\n");
//     fp6_mul_nonmod_montgomery(&tmp3_fpd6,&A->x0,&B->x0);//a*c
//         fp6_mod_montgomery_double(&out,&tmp3_fpd6);
//         //fp6_println("out_A_1=",&out);printf("\n");

//     //x0
//     fp6_mul_basis_nonmod_double(&buf_fpd6,&tmp1_fpd6);//b*d*v
//         fp6_mod_montgomery_double(&out,&buf_fpd6);
//         //fp6_println("out_A_2=",&out);printf("\n");
//     fp6_add_nonmod_double(&buf_fpd6,&buf_fpd6,&tmp3_fpd6);//a*c+b*d*v
//         fp6_mod_montgomery_double(&out,&buf_fpd6);
//         //fp6_println("out_A_3=",&out);printf("\n");
//     fp6_mod_montgomery_double(&ANS->x0,&buf_fpd6);
//         //fp6_println("chk1=",&ANS->x0);printf("\n");

//     //x1
//     fp6_sub_nonmod_double(&buf_fpd6,&tmp2_fpd6,&tmp1_fpd6);
//         fp6_mod_montgomery_double(&out,&buf_fpd6);
//         //fp6_println("out_A_4=",&out);printf("\n");
//     fp6_sub_nonmod_double(&buf_fpd6,&buf_fpd6,&tmp3_fpd6);
//         fp6_mod_montgomery_double(&out,&buf_fpd6);
//         //fp6_println("out_A_5=",&out);printf("\n");
//     fp6_mod_montgomery_double(&ANS->x1,&buf_fpd6);
//         //fp6_println("chk2=",&ANS->x1);printf("\n");

//         //under to debug
//     fp6_mul_lazy_montgomery(&tmp2_fp6,&A->x1,&B->x1);//b*d
//     fp6_add_nonmod_single(&tmp1_fp6,&A->x0,&A->x1);//a+b
//     fp6_add_nonmod_single(&ANS->x1,&B->x0,&B->x1);//c+d
//     fp6_mul_lazy_montgomery(&ANS->x1,&tmp1_fp6,&ANS->x1);//(a+b)(c+d)
//         //fp6_println("out_B_0=",&ANS->x1);printf("\n");
//     fp6_mul_lazy_montgomery(&tmp1_fp6,&A->x0,&B->x0);//a*c
//         //fp6_println("out_B_1=",&tmp1_fp6);printf("\n");

//     //x0
//     fp6_mul_basis(&ANS->x0,&tmp2_fp6);//b*d*v
//         //fp6_println("out_B_2=",&ANS->x0);printf("\n");
//     fp6_add(&ANS->x0,&ANS->x0,&tmp1_fp6);//a*c+b*d*v
//         //fp6_println("out_B_3=",&ANS->x0);printf("\n");
//         //fp6_println("chk1=",&ANS->x0);printf("\n");

//     //x1
//     fp6_sub(&ANS->x1,&ANS->x1,&tmp1_fp6);
//         //fp6_println("out_B_4=",&ANS->x1);printf("\n");
//     fp6_sub(&ANS->x1,&ANS->x1,&tmp2_fp6);
//         //fp6_println("out_B_5=",&ANS->x1);printf("\n");
//         //fp6_println("chk2=",&ANS->x1);printf("\n");
//         //getchar();

// }

void fp12_mul_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp6_mul_ui(&ANS->x0,&A->x0,UI);
    fp6_mul_ui(&ANS->x1,&A->x1,UI);
}

void fp12_mul_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B){
    fp6_mul_mpn(&ANS->x0,&A->x0,B);
    fp6_mul_mpn(&ANS->x1,&A->x1,B);
}
//complex
void fp12_sqr(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
    fp6_add(&tmp1_fp6,&A->x0,&A->x1);
    fp6_mul_basis(&tmp2_fp6,&A->x1);
    fp6_add(&tmp2_fp6,&tmp2_fp6,&A->x0);
    fp6_mul(&tmp3_fp6,&A->x0,&A->x1);

    //x0
    fp6_mul(&ANS->x0,&tmp1_fp6,&tmp2_fp6);
    fp6_sub(&ANS->x0,&ANS->x0,&tmp3_fp6);
    fp6_mul_basis(&tmp1_fp6,&tmp3_fp6);
    fp6_sub(&ANS->x0,&ANS->x0,&tmp1_fp6);

    //x1
    fp6_add(&ANS->x1,&tmp3_fp6,&tmp3_fp6);
}

/*
//TODO: Karat -> bug
void fp12_sqr(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
    fp6_add(&tmp1_fp6,&A->x0,&A->x1);

    fp6_mul(&tmp2_fp6,&A->x0,&A->x1);
    fp6_add(&tmp2_fp6,&tmp2_fp6,&tmp2_fp6);
    fp6_mul_basis(&tmp3_fp6,&tmp2_fp6);

    //x0
    fp6_sqr(&ANS->x0,&tmp1_fp6);
    fp6_sub(&ANS->x0,&ANS->x0,&tmp3_fp6);

    //x1
    fp6_set(&ANS->x1,&tmp2_fp6);
}
*/
//Complex
void fp12_sqr_lazy(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
    fp6_add_nonmod_single(&tmp1_fp6,&A->x0,&A->x1);
    fp6_mul_basis_nonmod_single(&tmp2_fp6,&A->x1);
    fp6_add_nonmod_single(&tmp2_fp6,&tmp2_fp6,&A->x0);
    fp6_mul_lazy(&tmp3_fp6,&A->x0,&A->x1);

    //x0
    fp6_mul_lazy(&ANS->x0,&tmp1_fp6,&tmp2_fp6);
    fp6_sub(&ANS->x0,&ANS->x0,&tmp3_fp6);
    fp6_mul_basis(&tmp1_fp6,&tmp3_fp6);
    fp6_sub(&ANS->x0,&ANS->x0,&tmp1_fp6);

    //x1
    fp6_add(&ANS->x1,&tmp3_fp6,&tmp3_fp6);
}
void fp12_sqr_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6;
    fp6_add_nonmod_single(&tmp1_fp6,&A->x0,&A->x1);
    fp6_mul_basis_nonmod_single(&tmp2_fp6,&A->x1);
    fp6_add_nonmod_single(&tmp2_fp6,&tmp2_fp6,&A->x0);
    fp6_mul_lazy_montgomery(&tmp3_fp6,&A->x0,&A->x1);

    //x0
    fp6_mul_lazy_montgomery(&ANS->x0,&tmp1_fp6,&tmp2_fp6);
    fp6_sub(&ANS->x0,&ANS->x0,&tmp3_fp6);
    fp6_mul_basis(&tmp1_fp6,&tmp3_fp6);
    fp6_sub(&ANS->x0,&ANS->x0,&tmp1_fp6);

    //x1
    fp6_add(&ANS->x1,&tmp3_fp6,&tmp3_fp6);
}
void fp12_sqr_cyclotomic(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6;
    //A=a+b*gamma in G3
    //A^2=(1+2b^2*beta)+((a+b)^2-1-b^2*beta-b^2)
    fp6_add(&tmp1_fp6,&A->x0,&A->x1);
    fp6_sqr(&tmp1_fp6,&tmp1_fp6);           //(a+b)^2
    fp6_sqr(&tmp2_fp6,&A->x1);          //b^2
    fp6_mul_basis(&ANS->x1,&tmp2_fp6);  //b^2*beta
    fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);   //1+2b^2*beta

    fp6_sub(&ANS->x1,&tmp1_fp6,&ANS->x1);
    fp6_sub(&ANS->x1,&ANS->x1,&tmp2_fp6);
    fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);   //(a+b)^2-1-b^2*beta-b^2
}
void fp12_sqr_cyclotomic_lazy(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6;
    fp6_add_nonmod_single(&tmp1_fp6,&A->x0,&A->x1);
    fp6_sqr_lazy(&tmp1_fp6,&tmp1_fp6);           //(a+b)^2
    fp6_sqr_lazy(&tmp2_fp6,&A->x1);          //b^2
    fp6_mul_basis(&ANS->x1,&tmp2_fp6);  //b^2*beta
    fp6_add(&ANS->x0,&ANS->x1,&ANS->x1);
    fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);   //1+2b^2*beta

    fp6_sub(&ANS->x1,&tmp1_fp6,&ANS->x1);
    fp6_sub(&ANS->x1,&ANS->x1,&tmp2_fp6);
    fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);   //(a+b)^2-1-b^2*beta-b^2
}
void fp12_sqr_Karabina(fp12_t *ANS,fp12_t *A){
    fp12_sqr_compressed(ANS,A);
    fp12_sqr_recover_g1(ANS,ANS);
    fp12_sqr_recover_g0(ANS,ANS);
}
void fp12_sqr_compressed(fp12_t *ANS,fp12_t *A){
    static fp2_t g2,g3,g4,g5;
    static fp2_t T0,T1,T2,T3;
    static fp2_t t0,t1,t2;
    static fp2_t B1;

    //set
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    //fp2_sqr(&T0,&A->x0.x1);
    fp2_sqr(&T0,&g4);
    //fp2_sqr(&T1,&A->x1.x2);
    fp2_sqr(&T1,&g5);

    fp2_mul_basis(&T2,&T1);

    fp2_add(&T2,&T2,&T0);

    fp2_set(&t2,&T2);

    //fp2_add(&t0,&A->x0.x1,&A->x1.x2);
    fp2_add(&t0,&g4,&g5);
    fp2_sqr(&T2,&t0);

    fp2_add(&T0,&T0,&T1);
    fp2_sub(&T2,&T2,&T0);

    fp2_set(&t0,&T2);
    fp2_add(&t1,&g2,&g3);
    fp2_sqr(&T3,&t1);
    fp2_sqr(&T2,&g2);

    fp2_mul_basis(&t1,&t0);

    fp2_add(&g2,&g2,&t1);
    fp2_add(&g2,&g2,&g2);

    fp2_add(&g2,&g2,&t1);
    fp2_sub(&t1,&t2,&g3);
    fp2_add(&t1,&t1,&t1);

    fp2_sqr(&T1,&g3);

    fp2_add(&g3,&t1,&t2);

    fp2_mul_basis(&T0,&T1);

    fp2_add(&T0,&T0,&T2);

    fp2_set(&t0,&T0);
    fp2_sub(&g4,&t0,&g4);//add->sub
    fp2_add(&g4,&g4,&g4);

    fp2_add(&g4,&g4,&t0);

    fp2_add(&T2,&T2,&T1);
    fp2_sub(&T3,&T3,&T2);//sub->add

    fp2_set(&t0,&T3);
    fp2_add(&g5,&g5,&t0);
    fp2_add(&g5,&g5,&g5);

    fp2_add(&g5,&g5,&t0);

    //set
    fp2_set_ui_ui(&ANS->x0.x0,0);
    fp2_set_ui_ui(&ANS->x1.x1,0);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);

}

void fp12_sqr_compressed_lazy(fp12_t *ANS,fp12_t *A){
    static fp2_t g2,g3,g4,g5;
    static fp2_t T0,T1,T2,T3;
    static fp2_t t0,t1,t2;
    static fp2_t B1;

    //set
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    fp2_sqr_lazy(&T0,&g4);
    fp2_sqr_lazy(&T1,&g5);

    fp2_mul_basis(&T2,&T1);

    fp2_add(&T2,&T2,&T0);

    fp2_set(&t2,&T2);

    fp2_add(&t0,&g4,&g5);
    fp2_sqr_lazy(&T2,&t0);

    fp2_add(&T0,&T0,&T1);
    fp2_sub(&T2,&T2,&T0);

    fp2_set(&t0,&T2);
    fp2_add(&t1,&g2,&g3);
    fp2_sqr_lazy(&T3,&t1);
    fp2_sqr_lazy(&T2,&g2);

    fp2_mul_basis(&t1,&t0);

    fp2_add(&g2,&g2,&t1);
    fp2_add(&g2,&g2,&g2);

    fp2_add(&g2,&g2,&t1);
    fp2_sub(&t1,&t2,&g3);
    fp2_add(&t1,&t1,&t1);

    fp2_sqr_lazy(&T1,&g3);

    fp2_add(&g3,&t1,&t2);

    fp2_mul_basis(&T0,&T1);

    fp2_add(&T0,&T0,&T2);

    fp2_set(&t0,&T0);
    fp2_sub(&g4,&t0,&g4);//add->sub
    fp2_add(&g4,&g4,&g4);

    fp2_add(&g4,&g4,&t0);

    fp2_add(&T2,&T2,&T1);
    fp2_sub(&T3,&T3,&T2);//sub->add

    fp2_set(&t0,&T3);
    fp2_add(&g5,&g5,&t0);
    fp2_add(&g5,&g5,&g5);

    fp2_add(&g5,&g5,&t0);

    //set
    fp2_set_ui_ui(&ANS->x0.x0,0);
    fp2_set_ui_ui(&ANS->x1.x1,0);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}
void fp12_sqr_compressed_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp2_t g2,g3,g4,g5;
    static fp2_t T0,T1,T2,T3;
    static fp2_t t0,t1,t2;
    static fp2_t B1;

    //set
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    fp2_sqr_lazy_montgomery(&T0,&g4);
    fp2_sqr_lazy_montgomery(&T1,&g5);

    fp2_mul_basis(&T2,&T1);

    fp2_add(&T2,&T2,&T0);

    fp2_set(&t2,&T2);

    fp2_add(&t0,&g4,&g5);
    fp2_sqr_lazy_montgomery(&T2,&t0);

    fp2_add(&T0,&T0,&T1);
    fp2_sub(&T2,&T2,&T0);

    fp2_set(&t0,&T2);
    fp2_add(&t1,&g2,&g3);
    fp2_sqr_lazy_montgomery(&T3,&t1);
    fp2_sqr_lazy_montgomery(&T2,&g2);

    fp2_mul_basis(&t1,&t0);

    fp2_add(&g2,&g2,&t1);
    fp2_add(&g2,&g2,&g2);

    fp2_add(&g2,&g2,&t1);
    fp2_sub(&t1,&t2,&g3);
    fp2_add(&t1,&t1,&t1);

    fp2_sqr_lazy_montgomery(&T1,&g3);

    fp2_add(&g3,&t1,&t2);

    fp2_mul_basis(&T0,&T1);

    fp2_add(&T0,&T0,&T2);

    fp2_set(&t0,&T0);
    fp2_sub(&g4,&t0,&g4);//add->sub
    fp2_add(&g4,&g4,&g4);

    fp2_add(&g4,&g4,&t0);

    fp2_add(&T2,&T2,&T1);
    fp2_sub(&T3,&T3,&T2);//sub->add

    fp2_set(&t0,&T3);
    fp2_add(&g5,&g5,&t0);
    fp2_add(&g5,&g5,&g5);

    fp2_add(&g5,&g5,&t0);

    //set
    fp2_set_ui_ui(&ANS->x0.x0,0);
    fp2_set_ui_ui(&ANS->x1.x1,0);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}

void fp12_sqr_recover_g1(fp12_t *ANS,fp12_t *A){
    static fp2_t g1,g2,g3,g4,g5;
    static fp2_t tmp,f;
    static fp2_t t0,t1;//g1=t0/t1
    static fp2_t C12,C02;

    //set
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    //if
    if(fp2_cmp_zero(&g2)==1){
        fp2_sqr(&tmp,&g5);
        fp2_mul_basis(&C12,&tmp);
        fp2_sqr(&C02,&g4);
        fp2_set(&t0,&C02);
        fp2_add(&C02,&C02,&C02);
        fp2_add(&C02,&C02,&t0);
        fp2_add(&t0,&C12,&C02);
        fp2_sub(&t0,&t0,&g3);
        fp2_sub(&t0,&t0,&g3);
        fp2_add(&t1,&g2,&g2);
        fp2_add(&t1,&t1,&t1);
    //else
    }else{
        fp2_mul(&t0,&g4,&g5);
        fp2_add(&t0,&t0,&t0);
        fp2_set(&t1,&g3);
    }

    fp2_inv(&t1,&t1);
    fp2_mul(&g1,&t0,&t1);

    //set
    fp2_set_ui_ui(&ANS->x0.x0,0);
    fp2_set(&ANS->x1.x1,&g1);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}
void fp12_sqr_recover_g1_noninv(fp12_t *ANS,fp12_t *A){
    static fp2_t g1,g2,g3,g4,g5;
    static fp2_t tmp,f;
    static fp2_t t0,t1;//g1=t0/t1
    static fp2_t C12,C02;

    //set
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    //if
    if(fp2_cmp_zero(&g2)==1){
        fp2_sqr_lazy_montgomery(&tmp,&g5);
        fp2_mul_basis(&C12,&tmp);
        fp2_sqr_lazy_montgomery(&C02,&g4);
        fp2_set(&t0,&C02);
        fp2_add(&C02,&C02,&C02);
        fp2_add(&C02,&C02,&t0);
        fp2_add(&t0,&C12,&C02);
        fp2_sub(&t0,&t0,&g3);
        fp2_sub(&t0,&t0,&g3);
        fp2_add(&t1,&g2,&g2);
        fp2_add(&t1,&t1,&t1);
    //else
    }else{
        fp2_mul_lazy_montgomery(&t0,&g4,&g5);
        fp2_add(&t0,&t0,&t0);
        fp2_set(&t1,&g3);
    }

    //set
    fp2_set(&ANS->x0.x0,&t0);
    fp2_set(&ANS->x1.x1,&t1);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}

void fp12_sqr_recover_g1_montrick_montgomery_ham3(fp12_t *A,fp12_t *B,fp12_t *C){
    static fp2_t a0,a1,b0,b1,c0,c1,ai,bi,ci;
    static fp2_t inv,all,buf,ab,bc,ca;
    static fp2_t g1,g2,g3,g4,g5;
    static fp2_t tmp,f;
    static fp2_t t0,t1;//g1=t0/t1
    static fp2_t C12,C02;

    fp2_set(&a0,&A->x0.x0);
    fp2_set(&a1,&A->x1.x1);
    fp2_set(&b0,&B->x0.x0);
    fp2_set(&b1,&B->x1.x1);
    fp2_set(&c0,&C->x0.x0);
    fp2_set(&c1,&C->x1.x1);

    fp2_mul_lazy_montgomery(&ab,&a1,&b1);
    fp2_mul_lazy_montgomery(&bc,&b1,&c1);
    fp2_mul_lazy_montgomery(&ca,&c1,&a1);

    fp2_mul_lazy_montgomery(&all,&ab,&c1);
    fp2_inv_lazy_montgomery(&inv,&all);

    fp2_mul_lazy_montgomery(&ai,&inv,&bc);
    fp2_mul_lazy_montgomery(&bi,&inv,&ca);
    fp2_mul_lazy_montgomery(&ci,&inv,&ab);

    fp2_mul_lazy_montgomery(&A->x1.x1,&ai,&a0);
    fp2_mul_lazy_montgomery(&B->x1.x1,&bi,&b0);
    fp2_mul_lazy_montgomery(&C->x1.x1,&ci,&c0);
}


void fp12_sqr_recover_g1_montrick_montgomery(fp12_t *A,fp12_t *B,fp12_t *C,fp12_t *D){
    static fp2_t a0,a1,b0,b1,c0,c1,d0,d1,ai,bi,ci,di;
    static fp2_t inv,all,buf,ab,bc,cd,da;
    static fp2_t g1,g2,g3,g4,g5;
    static fp2_t tmp,f;
    static fp2_t t0,t1;//g1=t0/t1
    static fp2_t C12,C02;

    fp2_set(&a0,&A->x0.x0);
    fp2_set(&a1,&A->x1.x1);
    fp2_set(&b0,&B->x0.x0);
    fp2_set(&b1,&B->x1.x1);
    fp2_set(&c0,&C->x0.x0);
    fp2_set(&c1,&C->x1.x1);
    fp2_set(&d0,&D->x0.x0);
    fp2_set(&d1,&D->x1.x1);

    fp2_mul_lazy_montgomery(&ab,&a1,&b1);
    fp2_mul_lazy_montgomery(&bc,&b1,&c1);
    fp2_mul_lazy_montgomery(&cd,&c1,&d1);
    fp2_mul_lazy_montgomery(&da,&d1,&a1);

    fp2_mul_lazy_montgomery(&all,&ab,&cd);
    fp2_inv_lazy_montgomery(&inv,&all);

    fp2_mul_lazy_montgomery(&buf,&b1,&cd);
    fp2_mul_lazy_montgomery(&ai,&inv,&buf);
    fp2_mul_lazy_montgomery(&buf,&c1,&da);
    fp2_mul_lazy_montgomery(&bi,&inv,&buf);
    fp2_mul_lazy_montgomery(&buf,&d1,&ab);
    fp2_mul_lazy_montgomery(&ci,&inv,&buf);
    fp2_mul_lazy_montgomery(&buf,&a1,&bc);
    fp2_mul_lazy_montgomery(&di,&inv,&buf);


    fp2_mul_lazy_montgomery(&A->x1.x1,&ai,&a0);
    fp2_mul_lazy_montgomery(&B->x1.x1,&bi,&b0);
    fp2_mul_lazy_montgomery(&C->x1.x1,&ci,&c0);
    fp2_mul_lazy_montgomery(&D->x1.x1,&di,&d0);
}

void fp12_sqr_recover_g1_montrick(fp12_t *A,fp12_t *B,fp12_t *C,fp12_t *D){
    static fp2_t a0,a1,b0,b1,c0,c1,d0,d1,ai,bi,ci,di;
    static fp2_t inv,all,buf,ab,bc,cd,da;
    static fp2_t g1,g2,g3,g4,g5;
    static fp2_t tmp,f;
    static fp2_t t0,t1;//g1=t0/t1
    static fp2_t C12,C02;

    fp2_set(&a0,&A->x0.x0);
    fp2_set(&a1,&A->x1.x1);
    fp2_set(&b0,&B->x0.x0);
    fp2_set(&b1,&B->x1.x1);
    fp2_set(&c0,&C->x0.x0);
    fp2_set(&c1,&C->x1.x1);
    fp2_set(&d0,&D->x0.x0);
    fp2_set(&d1,&D->x1.x1);

    //TODO:mul cut back
    fp2_mul(&ab,&a1,&b1);
    fp2_mul(&bc,&b1,&c1);
    fp2_mul(&cd,&c1,&d1);
    fp2_mul(&da,&d1,&a1);

    fp2_mul(&all,&ab,&cd);
    fp2_inv(&inv,&all);

    fp2_mul(&buf,&b1,&cd);
    fp2_mul(&ai,&inv,&buf);
    fp2_mul(&buf,&c1,&da);
    fp2_mul(&bi,&inv,&buf);
    fp2_mul(&buf,&d1,&ab);
    fp2_mul(&ci,&inv,&buf);
    fp2_mul(&buf,&a1,&bc);
    fp2_mul(&di,&inv,&buf);


    fp2_mul(&A->x1.x1,&ai,&a0);
    fp2_mul(&B->x1.x1,&bi,&b0);
    fp2_mul(&C->x1.x1,&ci,&c0);
    fp2_mul(&D->x1.x1,&di,&d0);
}
void fp12_sqr_recover_g1_lazy(fp12_t *ANS,fp12_t *A){
    static fp2_t g1,g2,g3,g4,g5;
    static fp2_t tmp,f;
    static fp2_t t0,t1;//g1=t0/t1
    static fp2_t C12,C02;

    //set
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    //if
    if(fp2_cmp_zero(&g2)==1){
        fp2_sqr_lazy(&tmp,&g5);
        fp2_mul_basis(&C12,&tmp);
        fp2_sqr_lazy(&C02,&g4);
        fp2_set(&t0,&C02);
        fp2_add(&C02,&C02,&C02);
        fp2_add(&C02,&C02,&t0);
        fp2_add(&t0,&C12,&C02);
        fp2_sub(&t0,&t0,&g3);
        fp2_sub(&t0,&t0,&g3);
        fp2_add(&t1,&g2,&g2);
        fp2_add(&t1,&t1,&t1);
    //else
    }else{
        fp2_mul_lazy(&t0,&g4,&g5);
        fp2_add(&t0,&t0,&t0);
        fp2_set(&t1,&g3);
    }

    fp2_inv_lazy(&t1,&t1);
    fp2_mul_lazy(&g1,&t0,&t1);

    //set
    fp2_set_ui_ui(&ANS->x0.x0,0);
    fp2_set(&ANS->x1.x1,&g1);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}
void fp12_sqr_recover_g1_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp2_t g1,g2,g3,g4,g5;
    static fp2_t tmp,f;
    static fp2_t t0,t1;//g1=t0/t1
    static fp2_t C12,C02;

    //set
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    //if
    if(fp2_cmp_zero(&g2)==1){
        fp2_sqr_lazy_montgomery(&tmp,&g5);
        fp2_mul_basis(&C12,&tmp);
        fp2_sqr_lazy_montgomery(&C02,&g4);
        fp2_set(&t0,&C02);
        fp2_add(&C02,&C02,&C02);
        fp2_add(&C02,&C02,&t0);
        fp2_add(&t0,&C12,&C02);
        fp2_sub(&t0,&t0,&g3);
        fp2_sub(&t0,&t0,&g3);
        fp2_add(&t1,&g2,&g2);
        fp2_add(&t1,&t1,&t1);
    //else
    }else{
        fp2_mul_lazy_montgomery(&t0,&g4,&g5);
        fp2_add(&t0,&t0,&t0);
        fp2_set(&t1,&g3);
    }

    fp2_inv_lazy_montgomery(&t1,&t1);
    fp2_mul_lazy_montgomery(&g1,&t0,&t1);

    //set
    fp2_set_ui_ui(&ANS->x0.x0,0);
    fp2_set(&ANS->x1.x1,&g1);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}
void fp12_sqr_recover_g0(fp12_t *ANS,fp12_t *A){
    static fp2_t g0,g1,g2,g3,g4,g5;
    static fp2_t one;
    static fp2_t t0,t1;
    static fp2_t C12,C02;

    fp2_set_ui(&one,1);

    //set
    fp2_set(&g0,&A->x0.x0);
    fp2_set(&g1,&A->x1.x1);
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    fp2_sqr(&t0,&g1);
    fp2_mul(&t1,&g3,&g4);
    fp2_sub(&t0,&t0,&t1);
    fp2_add(&t0,&t0,&t0);
    fp2_sub(&t0,&t0,&t1);
    fp2_mul(&t1,&g2,&g5);
    fp2_add(&t0,&t0,&t1);
    fp2_mul_basis(&g0,&t0);
    fp2_add(&g0,&g0,&one);

    //set
    fp2_set(&ANS->x0.x0,&g0);
    fp2_set(&ANS->x1.x1,&g1);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}

void fp12_sqr_recover_g0_lazy(fp12_t *ANS,fp12_t *A){
    static fp2_t g0,g1,g2,g3,g4,g5;
    static fp2_t one;
    static fp2_t t0,t1;
    static fp2_t C12,C02;

    fp2_set_ui(&one,1);

    //set
    fp2_set(&g0,&A->x0.x0);
    fp2_set(&g1,&A->x1.x1);
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    fp2_sqr_lazy(&t0,&g1);
    fp2_mul_lazy(&t1,&g3,&g4);
    fp2_sub(&t0,&t0,&t1);
    fp2_add(&t0,&t0,&t0);
    fp2_sub(&t0,&t0,&t1);
    fp2_mul_lazy(&t1,&g2,&g5);
    fp2_add(&t0,&t0,&t1);
    fp2_mul_basis(&g0,&t0);
    fp2_add(&g0,&g0,&one);

    //set
    fp2_set(&ANS->x0.x0,&g0);
    fp2_set(&ANS->x1.x1,&g1);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}
void fp12_sqr_recover_g0_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp2_t g0,g1,g2,g3,g4,g5;
    static fp2_t one;
    static fp2_t t0,t1;
    static fp2_t C12,C02;

    //OPTIMIZE:montomgery 1
    fp_set_mpn(&one.x0,RmodP);
    //fp2_to_montgomery(&one,&one);

    //set
    fp2_set(&g0,&A->x0.x0);
    fp2_set(&g1,&A->x1.x1);
    fp2_set(&g2,&A->x1.x0);
    fp2_set(&g3,&A->x0.x2);
    fp2_set(&g4,&A->x0.x1);
    fp2_set(&g5,&A->x1.x2);

    fp2_sqr_lazy_montgomery(&t0,&g1);
    fp2_mul_lazy_montgomery(&t1,&g3,&g4);
    fp2_sub(&t0,&t0,&t1);
    fp2_add(&t0,&t0,&t0);
    fp2_sub(&t0,&t0,&t1);
    fp2_mul_lazy_montgomery(&t1,&g2,&g5);
    fp2_add(&t0,&t0,&t1);
    fp2_mul_basis(&g0,&t0);
    fp2_add(&g0,&g0,&one);

    //set
    fp2_set(&ANS->x0.x0,&g0);
    fp2_set(&ANS->x1.x1,&g1);
    fp2_set(&ANS->x1.x0,&g2);
    fp2_set(&ANS->x0.x2,&g3);
    fp2_set(&ANS->x0.x1,&g4);
    fp2_set(&ANS->x1.x2,&g5);
}

void fp12_sqr_GS(fp12_t *ANS,fp12_t *A){
    static fp2_t z0,z1,z2,z3,z4,z5;
    static fp2_t t0,t1,t2,t3;
    static fp2_t tmp0,tmp1;

    //set
    fp2_set(&z0,&A->x0.x0);
    fp2_set(&z1,&A->x1.x1);
    fp2_set(&z2,&A->x1.x0);
    fp2_set(&z3,&A->x0.x2);
    fp2_set(&z4,&A->x0.x1);
    fp2_set(&z5,&A->x1.x2);

    fp4_sqr(&t0,&t1,&z0,&z1);

    fp2_sub(&z0,&t0,&z0);
    fp2_add(&z0,&z0,&z0);
    fp2_add(&z0,&z0,&t0);

    fp2_add(&z1,&t1,&z1);
    fp2_add(&z1,&z1,&z1);
    fp2_add(&z1,&z1,&t1);

    fp4_sqr(&t0,&t1,&z2,&z3);
    fp4_sqr(&t2,&t3,&z4,&z5);

    fp2_sub(&z4,&t0,&z4);
    fp2_add(&z4,&z4,&z4);
    fp2_add(&z4,&z4,&t0);

    fp2_add(&z5,&t1,&z5);
    fp2_add(&z5,&z5,&z5);
    fp2_add(&z5,&z5,&t1);

    fp2_mul_basis(&t0,&t3);
    fp2_add(&z2,&t0,&z2);
    fp2_add(&z2,&z2,&z2);
    fp2_add(&z2,&z2,&t0);

    fp2_sub(&z3,&t2,&z3);
    fp2_add(&z3,&z3,&z3);
    fp2_add(&z3,&z3,&t2);
    //set
    fp2_set(&ANS->x0.x0,&z0);
    fp2_set(&ANS->x1.x1,&z1);
    fp2_set(&ANS->x1.x0,&z2);
    fp2_set(&ANS->x0.x2,&z3);
    fp2_set(&ANS->x0.x1,&z4);
    fp2_set(&ANS->x1.x2,&z5);
}

void fp4_sqr(fp2_t *t0,fp2_t *t1,fp2_t *g0,fp2_t *g1){
    static fp2_t buf0,buf1;

    fp2_sqr(&buf0,g0);
    fp2_sqr(&buf1,g1);
    fp2_add(t1,g0,g1);
    fp2_mul_basis(t0,&buf1);
    fp2_add(t0,t0,&buf0);
    fp2_sqr(t1,t1);
    fp2_sub(t1,t1,&buf0);
    fp2_sub(t1,t1,&buf1);
}

void fp12_sqr_GS_lazy(fp12_t *ANS,fp12_t *A){
    static fp2_t z0,z1,z2,z3,z4,z5;
    static fp2_t t0,t1,t2,t3;
    static fp2_t tmp0,tmp1;

    //set
    fp2_set(&z0,&A->x0.x0);
    fp2_set(&z1,&A->x1.x1);
    fp2_set(&z2,&A->x1.x0);
    fp2_set(&z3,&A->x0.x2);
    fp2_set(&z4,&A->x0.x1);
    fp2_set(&z5,&A->x1.x2);

    fp4_sqr_lazy(&t0,&t1,&z0,&z1);

    fp2_sub(&z0,&t0,&z0);
    fp2_add(&z0,&z0,&z0);
    fp2_add(&z0,&z0,&t0);

    fp2_add(&z1,&t1,&z1);
    fp2_add(&z1,&z1,&z1);
    fp2_add(&z1,&z1,&t1);

    fp4_sqr_lazy(&t0,&t1,&z2,&z3);
    fp4_sqr_lazy(&t2,&t3,&z4,&z5);

    fp2_sub(&z4,&t0,&z4);
    fp2_add(&z4,&z4,&z4);
    fp2_add(&z4,&z4,&t0);

    fp2_add(&z5,&t1,&z5);
    fp2_add(&z5,&z5,&z5);
    fp2_add(&z5,&z5,&t1);

    fp2_mul_basis(&t0,&t3);
    fp2_add(&z2,&t0,&z2);
    fp2_add(&z2,&z2,&z2);
    fp2_add(&z2,&z2,&t0);

    fp2_sub(&z3,&t2,&z3);
    fp2_add(&z3,&z3,&z3);
    fp2_add(&z3,&z3,&t2);
    //set
    fp2_set(&ANS->x0.x0,&z0);
    fp2_set(&ANS->x1.x1,&z1);
    fp2_set(&ANS->x1.x0,&z2);
    fp2_set(&ANS->x0.x2,&z3);
    fp2_set(&ANS->x0.x1,&z4);
    fp2_set(&ANS->x1.x2,&z5);
}
void fp12_sqr_GS_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp2_t z0,z1,z2,z3,z4,z5;
    static fp2_t t0,t1,t2,t3;
    static fp2_t tmp0,tmp1;

    //set
    fp2_set(&z0,&A->x0.x0);
    fp2_set(&z1,&A->x1.x1);
    fp2_set(&z2,&A->x1.x0);
    fp2_set(&z3,&A->x0.x2);
    fp2_set(&z4,&A->x0.x1);
    fp2_set(&z5,&A->x1.x2);

    fp4_sqr_lazy_montgomery(&t0,&t1,&z0,&z1);

    fp2_sub(&z0,&t0,&z0);
    fp2_add(&z0,&z0,&z0);
    fp2_add(&z0,&z0,&t0);

    fp2_add(&z1,&t1,&z1);
    fp2_add(&z1,&z1,&z1);
    fp2_add(&z1,&z1,&t1);

    fp4_sqr_lazy_montgomery(&t0,&t1,&z2,&z3);
    fp4_sqr_lazy_montgomery(&t2,&t3,&z4,&z5);

    fp2_sub(&z4,&t0,&z4);
    fp2_add(&z4,&z4,&z4);
    fp2_add(&z4,&z4,&t0);

    fp2_add(&z5,&t1,&z5);
    fp2_add(&z5,&z5,&z5);
    fp2_add(&z5,&z5,&t1);

    fp2_mul_basis(&t0,&t3);
    fp2_add(&z2,&t0,&z2);
    fp2_add(&z2,&z2,&z2);
    fp2_add(&z2,&z2,&t0);

    fp2_sub(&z3,&t2,&z3);
    fp2_add(&z3,&z3,&z3);
    fp2_add(&z3,&z3,&t2);
    //set
    fp2_set(&ANS->x0.x0,&z0);
    fp2_set(&ANS->x1.x1,&z1);
    fp2_set(&ANS->x1.x0,&z2);
    fp2_set(&ANS->x0.x2,&z3);
    fp2_set(&ANS->x0.x1,&z4);
    fp2_set(&ANS->x1.x2,&z5);
}
void fp4_sqr_lazy(fp2_t *t0,fp2_t *t1,fp2_t *g0,fp2_t *g1){
    static fp2_t buf0,buf1;

    fp2_sqr_lazy(&buf0,g0);
    fp2_sqr_lazy(&buf1,g1);
    fp2_add_nonmod_single(t1,g0,g1);
    fp2_mul_basis(t0,&buf1);
    fp2_add(t0,t0,&buf0);
    fp2_sqr_lazy(t1,t1);
    fp2_sub(t1,t1,&buf0);
    fp2_sub(t1,t1,&buf1);
}
void fp4_sqr_lazy_montgomery(fp2_t *t0,fp2_t *t1,fp2_t *g0,fp2_t *g1){
    static fpd2_t buf0,buf1,t0d,t1d;

    fp2_sqr_nonmod_montgomery(&buf0,g0);
    fp2_sqr_nonmod_montgomery(&buf1,g1);
    fp2_add_nonmod_single(t1,g0,g1);
    fp2_mul_basis_nonmod_double(&t0d,&buf1);
    fp2_add_nonmod_double(&t0d,&t0d,&buf0);
    fp2_mod_montgomery_double(t0,&t0d);
    fp2_sqr_nonmod_montgomery(&t1d,t1);
    fp2_sub_nonmod_double(&t1d,&t1d,&buf0);
    fp2_sub_nonmod_double(&t1d,&t1d,&buf1);
    fp2_mod_montgomery_double(t1,&t1d);
}
void fp12_add(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp6_add(&ANS->x0,&A->x0,&B->x0);
    fp6_add(&ANS->x1,&A->x1,&B->x1);
}
void fp12_add_nonmod_single(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp6_add_nonmod_single(&ANS->x0,&A->x0,&B->x0);
    fp6_add_nonmod_single(&ANS->x1,&A->x1,&B->x1);
}
void fp12_add_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp6_add_ui(&ANS->x0,&A->x0,UI);
    fp6_add_ui(&ANS->x1,&A->x1,0);
}

void fp12_add_ui_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp6_add_ui_ui(&ANS->x0,&A->x0,UI);
    fp6_add_ui_ui(&ANS->x1,&A->x1,UI);
}
void fp12_add_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B){
    fp6_add_mpn(&ANS->x0,&ANS->x0,B);
    fp6_add_mpn(&ANS->x1,&ANS->x1,B);
}

void fp12_sub(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp6_sub(&ANS->x0,&A->x0,&B->x0);
    fp6_sub(&ANS->x1,&A->x1,&B->x1);
}

void fp12_sub_nonmod_single(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp6_sub_nonmod_single(&ANS->x0,&A->x0,&B->x0);
    fp6_sub_nonmod_single(&ANS->x1,&A->x1,&B->x1);
}
void fp12_sub_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp6_sub_ui(&ANS->x0,&ANS->x0,UI);
    fp6_sub_ui(&ANS->x1,&ANS->x1,0);
}

void fp12_sub_ui_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp6_sub_ui_ui(&ANS->x0,&ANS->x0,UI);
    fp6_sub_ui_ui(&ANS->x1,&ANS->x1,UI);
}
void fp12_sub_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B){
    fp6_sub_mpn(&ANS->x0,&ANS->x0,B);
    fp6_sub_mpn(&ANS->x1,&ANS->x1,B);
}

void fp12_inv(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6,tmp4_fp6;
    fp6_set(&tmp1_fp6,&A->x0);
    fp6_set_neg(&tmp2_fp6,&A->x1);

    fp6_sqr(&tmp3_fp6,&tmp1_fp6);
    fp6_mul(&tmp4_fp6,&tmp2_fp6,&A->x1);
    fp6_mul_basis(&tmp4_fp6,&tmp4_fp6);
    fp6_add(&tmp3_fp6,&tmp3_fp6,&tmp4_fp6);
    fp6_inv(&tmp3_fp6,&tmp3_fp6);
    fp6_mul(&ANS->x0,&tmp1_fp6,&tmp3_fp6);
    fp6_mul(&ANS->x1,&tmp2_fp6,&tmp3_fp6);
}
void fp12_inv_lazy(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6,tmp4_fp6;
    fp6_set(&tmp1_fp6,&A->x0);
    fp6_set_neg(&tmp2_fp6,&A->x1);

    fp6_mul_lazy(&tmp4_fp6,&tmp1_fp6,&A->x0);
    fp6_mul_lazy(&tmp4_fp6,&tmp2_fp6,&A->x1);
    fp6_mul_basis(&tmp4_fp6,&tmp4_fp6);
    fp6_add(&tmp3_fp6,&tmp3_fp6,&tmp4_fp6);
    fp6_inv_lazy(&tmp3_fp6,&tmp3_fp6);
    fp6_mul_lazy(&ANS->x0,&tmp1_fp6,&tmp3_fp6);
    fp6_mul_lazy(&ANS->x1,&tmp2_fp6,&tmp3_fp6);
}
void fp12_inv_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp6_t tmp1_fp6,tmp2_fp6,tmp3_fp6,tmp4_fp6;
    fp6_set(&tmp1_fp6,&A->x0);
    fp6_set_neg(&tmp2_fp6,&A->x1);
    //sqr
    fp6_mul_lazy_montgomery(&tmp3_fp6,&tmp1_fp6,&A->x0);
    fp6_mul_lazy_montgomery(&tmp4_fp6,&tmp2_fp6,&A->x1);
    fp6_mul_basis(&tmp4_fp6,&tmp4_fp6);
    fp6_add(&tmp3_fp6,&tmp3_fp6,&tmp4_fp6);
    fp6_inv_lazy_montgomery(&tmp3_fp6,&tmp3_fp6);
    fp6_mul_lazy_montgomery(&ANS->x0,&tmp1_fp6,&tmp3_fp6);
    fp6_mul_lazy_montgomery(&ANS->x1,&tmp2_fp6,&tmp3_fp6);
}
int  fp12_legendre(fp12_t *A){
    mpz_t exp;
    mpz_init(exp);
    fp12_t tmp;
    fp12_init(&tmp);

    mpz_pow_ui(exp,prime_z,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp12_pow(&tmp,A,exp);

    if(fp12_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

void fp12_sqrt(fp12_t *ANS,fp12_t *A){
    fp12_t x,y,t,k,n,tmp;
    fp12_init(&x);
    fp12_init(&y);
    fp12_init(&t);
    fp12_init(&k);
    fp12_init(&n);
    fp12_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);

    fp12_set_random(&n,state);
    while(fp12_legendre(&n)!=-1){
        fp12_set_random(&n,state);
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
    fp12_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp12_pow(&x,A,exp);
    fp12_mul(&tmp,&x,&x);
    fp12_mul(&k,&tmp,A);
    fp12_mul(&x,&x,A);
    while(fp12_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fp12_pow(&tmp,&k,exp);
        while(fp12_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fp12_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fp12_pow(&t,&y,result);
        fp12_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fp12_mul(&x,&x,&t);
        fp12_mul(&k,&k,&y);
    }
    fp12_set(ANS,&x);

    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
}

void fp12_pow(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp12_t tmp;
    fp12_init(&tmp);
    fp12_set(&tmp,A);

    for(i=1;i<length; i++){
        fp12_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            fp12_mul(&tmp,A,&tmp);
        }
    }

    fp12_set(ANS,&tmp);
}

void fp12_pow_lazy_montgomery(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp12_t tmp;
    fp12_init(&tmp);
    fp12_set(&tmp,A);

    for(i=1;i<length; i++){
        fp12_sqr_lazy_montgomery(&tmp,&tmp);
        if(binary[i]=='1'){
            fp12_mul_lazy_montgomery(&tmp,A,&tmp);
        }
    }

    fp12_set(ANS,&tmp);
}

int  fp12_cmp(fp12_t *A,fp12_t *B){
    if(!fp6_cmp(&A->x0,&B->x0) && !fp6_cmp(&A->x1,&B->x1)) return 0;
    else return 1;
}

int  fp12_cmp_ui(fp12_t *A,unsigned long int UI){
    if(!fp6_cmp_ui(&A->x0,UI) && !fp6_cmp_ui(&A->x1,UI)) return 0;
    else return 1;
}

int  fp12_cmp_mpn(fp12_t *A,mp_limb_t *B){
    if(!fp6_cmp_mpn(&A->x0,B) && !fp6_cmp_mpn(&A->x1,B)) return 0;
    else return 1;
}

int  fp12_cmp_zero(fp12_t *A){
    if(!fp6_cmp_zero(&A->x0) && !fp6_cmp_zero(&A->x1)) return 0;
    else return 1;
}

int  fp12_cmp_one(fp12_t *A){
    if(!fp6_cmp_one(&A->x0) && !fp6_cmp_zero(&A->x1)) return 0;
    else return 1;
}

void fp12_frobenius_map_p1(fp12_t *ANS,fp12_t *A){
    static fp_t tmp1_fp;
    //x0
    fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);//
    fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);//
    fp_set(&tmp1_fp,&A->x0.x1.x0);//
    fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    fp_set(&ANS->x0.x1.x1,&tmp1_fp);
    fp2_mul_mpn(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant[f_p1][1].x1.x0);
    fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    fp2_mul_mpn(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant[f_p1][2].x0.x0);
    //x1
    fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p1][3]);
    fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p1][4].x0.x0);
    fp_add(&tmp1_fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_set(&ANS->x1.x1.x1,&tmp1_fp);

    fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p1][5]);
}
void fp12_frobenius_map_p1_lazy(fp12_t *ANS,fp12_t *A){
    static fp_t tmp1_fp;
    //x0
    fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    fp_set(&tmp1_fp,&A->x0.x1.x0);
    fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    fp_set(&ANS->x0.x1.x1,&tmp1_fp);
    fp2_mul_mpn(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant[f_p1][1].x1.x0);
    fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    fp2_mul_mpn(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant[f_p1][2].x0.x0);
    //x1
    fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    fp2_mul_lazy(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p1][3]);
    fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p1][4].x0.x0);
    fp_add(&tmp1_fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_set(&ANS->x1.x1.x1,&tmp1_fp);

    fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    fp2_mul_lazy(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p1][5]);
}
void fp12_frobenius_map_p1_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp_t tmp1_fp;
    //x0
    fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    fp_set(&tmp1_fp,&A->x0.x1.x0);
    fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    fp_set(&ANS->x0.x1.x1,&tmp1_fp);
    fp2_mul_mpn_montgomery(&ANS->x0.x1,&ANS->x0.x1,frobenius_constant_montgomery[f_p1][1].x1.x0);
    fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    fp2_mul_mpn_montgomery(&ANS->x0.x2,&ANS->x0.x2,frobenius_constant_montgomery[f_p1][2].x0.x0);
    //x1
    fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    fp2_mul_lazy_montgomery(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant_montgomery[f_p1][3]);
    fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    fp2_mul_mpn_montgomery(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant_montgomery[f_p1][4].x0.x0);
    fp_add(&tmp1_fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_set(&ANS->x1.x1.x1,&tmp1_fp);

    fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    fp2_mul_lazy_montgomery(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant_montgomery[f_p1][5]);
}
void fp12_frobenius_map_p2(fp12_t *ANS,fp12_t *A){
    //x0
    fp2_set(&ANS->x0.x0,&A->x0.x0);
    fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p2][1].x0.x0);
    fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p2][2].x0.x0);
    //x1
    fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p2][3].x0.x0);
    fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p2][5].x0.x0);
}
void fp12_frobenius_map_p2_montgomery(fp12_t *ANS,fp12_t *A){
    //x0
    fp2_set(&ANS->x0.x0,&A->x0.x0);
    fp2_mul_mpn_montgomery(&ANS->x0.x1,&A->x0.x1,frobenius_constant_montgomery[f_p2][1].x0.x0);
    fp2_mul_mpn_montgomery(&ANS->x0.x2,&A->x0.x2,frobenius_constant_montgomery[f_p2][2].x0.x0);
    //x1
    fp2_mul_mpn_montgomery(&ANS->x1.x0,&A->x1.x0,frobenius_constant_montgomery[f_p2][3].x0.x0);
    fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    fp2_mul_mpn_montgomery(&ANS->x1.x2,&A->x1.x2,frobenius_constant_montgomery[f_p2][5].x0.x0);
}

void fp12_frobenius_map_p3(fp12_t *ANS,fp12_t *A){
    static fp_t tmp1_fp;
    //x0
    fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    fp_set(&tmp1_fp,&A->x0.x1.x0);
    fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    fp_set(&ANS->x0.x1.x1,&tmp1_fp);
    fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p3][3]);
    fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p3][4].x0.x0);
    fp_add(&tmp1_fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_set(&ANS->x1.x1.x1,&tmp1_fp);
    fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p3][5]);
}
void fp12_frobenius_map_p3_lazy(fp12_t *ANS,fp12_t *A){
    static fp_t tmp1_fp;
    //x0
    fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    fp_set(&tmp1_fp,&A->x0.x1.x0);
    fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    fp_set(&ANS->x0.x1.x1,&tmp1_fp);
    fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    fp2_mul_lazy(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant[f_p3][3]);
    fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    fp2_mul_mpn(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant[f_p3][4].x0.x0);
    fp_add(&tmp1_fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_set(&ANS->x1.x1.x1,&tmp1_fp);
    fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    fp2_mul_lazy(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant[f_p3][5]);
}
void fp12_frobenius_map_p3_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp_t tmp1_fp;
    //x0
    fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    fp_set(&tmp1_fp,&A->x0.x1.x0);
    fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    fp_set(&ANS->x0.x1.x1,&tmp1_fp);
    fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    fp2_mul_lazy_montgomery(&ANS->x1.x0,&ANS->x1.x0,&frobenius_constant_montgomery[f_p3][3]);
    fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    fp2_mul_mpn_montgomery(&ANS->x1.x1,&ANS->x1.x1,frobenius_constant_montgomery[f_p3][4].x0.x0);
    fp_add(&tmp1_fp,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_sub(&ANS->x1.x1.x0,&ANS->x1.x1.x0,&ANS->x1.x1.x1);
    fp_set(&ANS->x1.x1.x1,&tmp1_fp);
    fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    fp2_mul_lazy_montgomery(&ANS->x1.x2,&ANS->x1.x2,&frobenius_constant_montgomery[f_p3][5]);
}

void fp12_frobenius_map_p4(fp12_t *ANS,fp12_t *A){
    //x0
    fp2_set(&ANS->x0.x0,&A->x0.x0);
    fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p4][1].x0.x0);
    fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p4][2].x0.x0);
    //x1
    fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p4][3].x0.x0);
    fp2_set(&ANS->x1.x1,&A->x1.x1);
    fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p4][5].x0.x0);
}

void fp12_frobenius_map_p6(fp12_t *ANS,fp12_t *A){
    //x0
    fp6_set(&ANS->x0,&A->x0);
    //x1
    fp6_set_neg(&ANS->x1,&A->x1);
}

void fp12_frobenius_map_p6_montgomery(fp12_t *ANS,fp12_t *A){
    //x0
    fp6_set(&ANS->x0,&A->x0);
    //x1
    fp6_set_neg(&ANS->x1,&A->x1);
}

void fp12_frobenius_map_p8(fp12_t *ANS,fp12_t *A){
    //x0
    fp2_set(&ANS->x0.x0,&A->x0.x0);
    fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p8][1].x0.x0);
    fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p8][2].x0.x0);
    //x1
    fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p8][3].x0.x0);
    fp2_set(&ANS->x1.x1,&A->x1.x1);
    fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p8][5].x0.x0);
}

void fp12_frobenius_map_p10(fp12_t *ANS,fp12_t *A){
    //x0
    fp2_set(&ANS->x0.x0,&A->x0.x0);
    fp2_mul_mpn(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p10][1].x0.x0);
    fp2_mul_mpn(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p10][2].x0.x0);
    //x1
    fp2_mul_mpn(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p10][3].x0.x0);
    fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    fp2_mul_mpn(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p10][5].x0.x0);
}

int fp12_montgomery_trick_montgomery(fp12_t *A_inv,fp12_t *A,int n){
    int i;
    fp12_t ANS[n],ALL_inv;
    fp12_set(ANS,A);
    fp12_t check;

    for(i=1;i<n;i++){
        fp12_mul_lazy_montgomery(&ANS[i],&ANS[i-1],&A[i]);
    }
    fp12_inv_lazy_montgomery(&ALL_inv,&ANS[n-1]);
    for(i=n-1;i>0;i--){
        fp12_mul_lazy_montgomery(&A_inv[i],&ALL_inv,&ANS[i-1]);
        fp12_mul_lazy_montgomery(&ALL_inv,&ALL_inv,&A[i]);
    }

    fp12_set(A_inv,&ALL_inv);
    return 0;
}

// int fp12_legendre_sqrt(fp12_t *ANS,fp12_t *A){
//     fp12_t C,D,A_tmp;
//     int i;

//     //legendre
//     fp12_pow(&C,A,fp12_sqrt_power_z);
//     fp12_mul(&D,&C,&C);
//     fp12_mul(&D,&D,A);

//     if(mpn_cmp_ui(D.x0.x0.x0.x0,FPLIMB,1)==0)        i=1;
//     else if(mpn_cmp_ui(D.x0.x0.x0.x0,FPLIMB,0)==0)    return 0;
//     else                    return -1;

//     //sqrt
//     fp_mul(ANS,&C,A);
//     return 1;
// }