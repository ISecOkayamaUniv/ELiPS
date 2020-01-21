#include <ELiPS/fp6.h>
//fp6_t
void fp6_init(fp6_t *A){
    fp2_init(&A->x0);
    fp2_init(&A->x1);
    fp2_init(&A->x2);
}

void fp6_printf(char *str,fp6_t *A){
    gmp_printf("%s(",str);
    fp2_printf("",&A->x0);
    gmp_printf(",");
    fp2_printf("",&A->x1);
    gmp_printf(",");
    fp2_printf("",&A->x2);
    gmp_printf(")");
}

void fp6_println(char *str,fp6_t *A){
    gmp_printf("%s(",str);
    fp2_printf("",&A->x0);
    gmp_printf(",");
    fp2_printf("",&A->x1);
    gmp_printf(",");
    fp2_printf("",&A->x2);
    gmp_printf(")\n");
}
void fp6_printf_montgomery(char *str,fp6_t *A){
    gmp_printf("%s(",str);
    fp2_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp2_printf_montgomery("",&A->x1);
    gmp_printf(",");
    fp2_printf_montgomery("",&A->x2);
    gmp_printf(")");
}
void fp6_set(fp6_t *ANS,fp6_t *A){
    fp2_set(&ANS->x0,&A->x0);
    fp2_set(&ANS->x1,&A->x1);
    fp2_set(&ANS->x2,&A->x2);
}

void fp6_set_ui(fp6_t *ANS,unsigned long int UI){
    fp2_set_ui(&ANS->x0,UI);
    fp2_set_ui(&ANS->x1,0);
    fp2_set_ui(&ANS->x2,0);
}

void fp6_set_ui_ui(fp6_t *ANS,unsigned long int UI){
    fp2_set_ui(&ANS->x0,UI);
    fp2_set_ui(&ANS->x1,UI);
    fp2_set_ui(&ANS->x2,UI);
}
void fp6_set_mpn(fp6_t *ANS,mp_limb_t *A){
    fp2_set_mpn(&ANS->x0,A);
    fp2_set_ui(&ANS->x1,0);
    fp2_set_ui(&ANS->x2,0);
}

void fp6_set_neg(fp6_t *ANS,fp6_t *A){
    fp2_set_neg(&ANS->x0,&A->x0);
    fp2_set_neg(&ANS->x1,&A->x1);
    fp2_set_neg(&ANS->x2,&A->x2);
}
void fp6_to_montgomery(fp6_t *ANS,fp6_t *A){
    fp2_to_montgomery(&ANS->x0,&A->x0);
    fp2_to_montgomery(&ANS->x1,&A->x1);
    fp2_to_montgomery(&ANS->x2,&A->x2);
}
void fp6_mod_montgomery(fp6_t *ANS,fp6_t *A){
    fp2_mod_montgomery(&ANS->x0,&A->x0);
    fp2_mod_montgomery(&ANS->x1,&A->x1);
    fp2_mod_montgomery(&ANS->x2,&A->x2);
}
void fp6_mod_montgomery_double(fp6_t *ANS,fpd6_t *A){
    fp2_mod_montgomery_double(&ANS->x0,&A->x0);
    fp2_mod_montgomery_double(&ANS->x1,&A->x1);
    fp2_mod_montgomery_double(&ANS->x2,&A->x2);
}
void fp6_set_random(fp6_t *ANS,gmp_randstate_t state){
    fp2_set_random(&ANS->x0,state);
    fp2_set_random(&ANS->x1,state);
    fp2_set_random(&ANS->x2,state);
}

void fp6_mul(fp6_t *ANS,fp6_t *A,fp6_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2;
    //set
    fp2_mul(&tmp1_fp2,&A->x0,&B->x0);//x0*y0
    fp2_mul(&tmp2_fp2,&A->x1,&B->x1);//x1*y1
    fp2_mul(&tmp3_fp2,&A->x2,&B->x2);//x2*y2
    
    fp2_add(&tmp5_fp2,&A->x0,&A->x1);//x0+x1
    fp2_add(&tmp4_fp2,&B->x0,&B->x1);//y0+y1
    fp2_mul(&tmp5_fp2,&tmp5_fp2,&tmp4_fp2);//(x0+x1)(y0+y1)
    
    fp2_add(&tmp6_fp2,&A->x1,&A->x2);//x1+x2
    fp2_add(&tmp4_fp2,&B->x1,&B->x2);//y1+y2
    fp2_mul(&tmp6_fp2,&tmp6_fp2,&tmp4_fp2);//(x1+x2)(y1+y2)
    
    fp2_add(&tmp7_fp2,&B->x0,&B->x2);//y2+y0
    fp2_add(&tmp4_fp2,&A->x0,&A->x2);//x2+x0
    fp2_mul(&tmp7_fp2,&tmp7_fp2,&tmp4_fp2);//(x2+x0)(y2+y0)
    //x0
    fp2_sub(&tmp6_fp2,&tmp6_fp2,&tmp2_fp2);
    fp2_sub(&tmp6_fp2,&tmp6_fp2,&tmp3_fp2);//(x1+x2)(y1+y2)-x1y1-x2y2
    fp2_mul_basis(&tmp4_fp2,&tmp6_fp2);
    fp2_add(&ANS->x0,&tmp1_fp2,&tmp4_fp2);
    //x1
    fp2_sub(&tmp5_fp2,&tmp5_fp2,&tmp1_fp2);
    fp2_sub(&tmp5_fp2,&tmp5_fp2,&tmp2_fp2);
    fp2_mul_basis(&tmp4_fp2,&tmp3_fp2);
    fp2_add(&ANS->x1,&tmp4_fp2,&tmp5_fp2);
    //x2
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp1_fp2);
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp3_fp2);
    fp2_add(&ANS->x2,&tmp2_fp2,&tmp7_fp2);
}
/*
//karat
void fp6_mul(fp6_t *ANS,fp6_t *A,fp6_t *B){
    static fp2_t v0,v1,v2,a0a1,a0a2,a1a2,b0b1,b0b2,b1b2,buf;
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2;
    //set
    fp2_mul(&v0,&A->x0,&B->x0);//x0*y0
    fp2_mul(&v1,&A->x1,&B->x1);//x1*y1
    fp2_mul(&v2,&A->x2,&B->x2);//x2*y2
    
    fp2_add(&a0a1,&A->x0,&A->x1);//x0+x1
    fp2_add(&b0b1,&B->x0,&B->x1);//y0+y1
    fp2_mul(&tmp1_fp2,&a0a1,&b0b1);//(x0+x1)(y0+y1)
    
    fp2_add(&a1a2,&A->x1,&A->x2);//x1+x2
    fp2_add(&b1b2,&B->x1,&B->x2);//y1+y2
    fp2_mul(&tmp2_fp2,a1a2,b1b2);//(x1+x2)(y1+y2)
    
    fp2_add(&a0a2,&B->x0,&B->x2);//y2+y0
    fp2_add(&b0b2,&A->x0,&A->x2);//x2+x0
    fp2_mul(&tmp3_fp2,&a0a2,&b0b2);//(x2+x0)(y2+y0)
    
    //x0
    fp2_sub(&buf,&tmp2_fp2,&v1);
    fp2_sub(&buf,&buf,&v2);//(x1+x2)(y1+y2)-x1y1-x2y2
    fp2_mul_basis(&buf,&buf);
    fp2_add(&ANS->x0,&v0,&buf);
    
    //mada
    //x1
    fp2_sub(&buf,&tmp1_fp2,&tmp1_fp2);
    fp2_sub(&buf,&tmp5_fp2,&tmp2_fp2);
    fp2_mul_basis(&tmp4_fp2,&tmp3_fp2);
    fp2_add(&ANS->x1,&tmp4_fp2,&tmp5_fp2);
    //x2
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp1_fp2);
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp3_fp2);
    fp2_add(&ANS->x2,&tmp2_fp2,&tmp7_fp2);
}
*/
void fp6_mul_lazy(fp6_t *ANS,fp6_t *A,fp6_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2;
    //set
    fp2_mul_lazy(&tmp1_fp2,&A->x0,&B->x0);//tmp1
    fp2_mul_lazy(&tmp2_fp2,&A->x1,&B->x1);//tmp2
    fp2_mul_lazy(&tmp3_fp2,&A->x2,&B->x2);//tmp3
    
    fp2_add_nonmod_single(&tmp5_fp2,&A->x0,&A->x1);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x0,&B->x1);
    fp2_mul_lazy(&tmp5_fp2,&tmp5_fp2,&tmp4_fp2);//tmp5

    fp2_add_nonmod_single(&tmp6_fp2,&A->x1,&A->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x1,&B->x2);
    fp2_mul_lazy(&tmp6_fp2,&tmp6_fp2,&tmp4_fp2);//tmp6
    
    fp2_add_nonmod_single(&tmp7_fp2,&B->x0,&B->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&A->x0,&A->x2);
    fp2_mul_lazy(&tmp7_fp2,&tmp7_fp2,&tmp4_fp2);//tmp7
    //x0
    fp2_sub(&tmp6_fp2,&tmp6_fp2,&tmp2_fp2);
    fp2_sub(&tmp6_fp2,&tmp6_fp2,&tmp3_fp2);
    fp2_mul_basis(&tmp4_fp2,&tmp6_fp2);
    fp2_add(&ANS->x0,&tmp1_fp2,&tmp4_fp2);
    //x1
    fp2_sub(&tmp5_fp2,&tmp5_fp2,&tmp1_fp2);
    fp2_sub(&tmp5_fp2,&tmp5_fp2,&tmp2_fp2);
    fp2_mul_basis(&tmp4_fp2,&tmp3_fp2);
    fp2_add(&ANS->x1,&tmp4_fp2,&tmp5_fp2);
    //x2
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp1_fp2);
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp3_fp2);
    fp2_add(&ANS->x2,&tmp2_fp2,&tmp7_fp2);
}
void fp6_mul_lazy_montgomery2(fp6_t *ANS,fp6_t *A,fp6_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2;
    //set
    	//fp2_println("tmpoldA->x0=",&A->x0);
    	//fp2_println("tmpoldB->x0=",&B->x0);
    fp2_mul_lazy_montgomery(&tmp1_fp2,&A->x0,&B->x0);//tmp1
    	//fp2_println("tmpold1=",&tmp1_fp2);printf("\n\n");
    fp2_mul_lazy_montgomery(&tmp2_fp2,&A->x1,&B->x1);//tmp2
    	//fp2_println("tmpold2=",&tmp2_fp2);printf("\n\n");
    fp2_mul_lazy_montgomery(&tmp3_fp2,&A->x2,&B->x2);//tmp3
    	//fp2_println("tmpold3=",&tmp3_fp2);printf("\n\n");
    
    fp2_add_nonmod_single(&tmp5_fp2,&A->x0,&A->x1);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x0,&B->x1);
    fp2_mul_lazy_montgomery(&tmp5_fp2,&tmp5_fp2,&tmp4_fp2);//tmp5
    	//fp2_println("tmpold5=",&tmp5_fp2);printf("\n\n");

    fp2_add_nonmod_single(&tmp6_fp2,&A->x1,&A->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x1,&B->x2);
    fp2_mul_lazy_montgomery(&tmp6_fp2,&tmp6_fp2,&tmp4_fp2);//tmp6
    	//fp2_println("tmpold6=",&tmp6_fp2);printf("\n\n");

    
    fp2_add_nonmod_single(&tmp7_fp2,&B->x0,&B->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&A->x0,&A->x2);
    fp2_mul_lazy_montgomery(&tmp7_fp2,&tmp7_fp2,&tmp4_fp2);//tmp7
    	//fp2_println("tmpold7=",&tmp7_fp2); printf("\n\n");
//printf("chk\n\n");
    //x0
    fp2_sub(&tmp6_fp2,&tmp6_fp2,&tmp2_fp2);
    fp2_sub(&tmp6_fp2,&tmp6_fp2,&tmp3_fp2);
    fp2_mul_basis(&tmp4_fp2,&tmp6_fp2);
    	//fp2_println("tmpoldx0=",&tmp4_fp2);printf("\n\n");
    	//fp2_println("tmpoldx0'=",&tmp1_fp2);printf("\n\n");
    fp2_add(&ANS->x0,&tmp1_fp2,&tmp4_fp2);
    //x1
    fp2_sub(&tmp5_fp2,&tmp5_fp2,&tmp1_fp2);
    fp2_sub(&tmp5_fp2,&tmp5_fp2,&tmp2_fp2);
    fp2_mul_basis(&tmp4_fp2,&tmp3_fp2);
    	//fp2_println("tmpoldx1=",&tmp4_fp2);printf("\n\n");
    	//fp2_println("tmpoldx1'=",&tmp5_fp2);printf("\n\n");
    fp2_add(&ANS->x1,&tmp4_fp2,&tmp5_fp2);
    //x2
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp1_fp2);
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp3_fp2);
    	//fp2_println("tmpoldx2=",&tmp7_fp2);printf("\n\n");
    	//fp2_println("tmpoldx2'=",&tmp2_fp2);printf("\n\n");
    fp2_add(&ANS->x2,&tmp2_fp2,&tmp7_fp2);
    	//fp2_println("tmpoldx2''=",&ANS->x2);printf("\n\n");
    //getchar();
}
void fp6_mul_lazy_montgomery(fp6_t *ANS,fp6_t *A,fp6_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2;
    static fpd2_t tmp1_fpd2,tmp2_fpd2,tmp3_fpd2,tmp4_fpd2,tmp5_fpd2,tmp6_fpd2,tmp7_fpd2,buf_fpd2;
    fp2_t out;
    
    //set
    fp2_mul_nonmod_montgomery(&tmp1_fpd2,&A->x0,&B->x0);//tmp1
    fp2_mul_nonmod_montgomery(&tmp2_fpd2,&A->x1,&B->x1);//tmp2
    fp2_mul_nonmod_montgomery(&tmp3_fpd2,&A->x2,&B->x2);//tmp3
    
    fp2_add_nonmod_single(&tmp5_fp2,&A->x0,&A->x1);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x0,&B->x1);
    fp2_mul_nonmod_montgomery(&tmp5_fpd2,&tmp5_fp2,&tmp4_fp2);//tmp5

    fp2_add_nonmod_single(&tmp6_fp2,&A->x1,&A->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x1,&B->x2);
    fp2_mul_nonmod_montgomery(&tmp6_fpd2,&tmp6_fp2,&tmp4_fp2);//tmp6
    
    fp2_add_nonmod_single(&tmp7_fp2,&B->x0,&B->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&A->x0,&A->x2);
    fp2_mul_nonmod_montgomery(&tmp7_fpd2,&tmp7_fp2,&tmp4_fp2);//tmp7
    
    //x0
    fp2_sub_nonmod_double(&tmp6_fpd2,&tmp6_fpd2,&tmp2_fpd2);
    fp2_sub_nonmod_double(&tmp6_fpd2,&tmp6_fpd2,&tmp3_fpd2);
    fp2_mul_basis_lazy_double(&tmp4_fpd2,&tmp6_fpd2);
    	
    fp2_add_nonmod_double(&buf_fpd2,&tmp1_fpd2,&tmp4_fpd2);
    fp2_mod_montgomery_double(&ANS->x0,&buf_fpd2);
    
    //x1
    fp2_sub_nonmod_double(&tmp5_fpd2,&tmp5_fpd2,&tmp1_fpd2);
    fp2_sub_nonmod_double(&tmp5_fpd2,&tmp5_fpd2,&tmp2_fpd2);
    fp2_mul_basis_lazy_double(&tmp4_fpd2,&tmp3_fpd2);
    fp2_add_nonmod_double(&buf_fpd2,&tmp4_fpd2,&tmp5_fpd2);
    fp2_mod_montgomery_double(&ANS->x1,&buf_fpd2);
    
    //x2
    fp2_sub_nonmod_double(&tmp7_fpd2,&tmp7_fpd2,&tmp1_fpd2);
    fp2_sub_nonmod_double(&tmp7_fpd2,&tmp7_fpd2,&tmp3_fpd2);
    fp2_add_nonmod_double(&buf_fpd2,&tmp2_fpd2,&tmp7_fpd2);
    fp2_mod_montgomery_double(&ANS->x2,&buf_fpd2);
    
}
void fp6_mul_nonmod_montgomery(fpd6_t *ANS,fp6_t *A,fp6_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2;
    static fpd2_t tmp1_fpd2,tmp2_fpd2,tmp3_fpd2,tmp4_fpd2,tmp5_fpd2,tmp6_fpd2,tmp7_fpd2,buf_fpd2;
    fp2_t out;
    
    //set
    fp2_mul_nonmod_montgomery(&tmp1_fpd2,&A->x0,&B->x0);//tmp1
    fp2_mul_nonmod_montgomery(&tmp2_fpd2,&A->x1,&B->x1);//tmp2
    fp2_mul_nonmod_montgomery(&tmp3_fpd2,&A->x2,&B->x2);//tmp3
    
    fp2_add_nonmod_single(&tmp5_fp2,&A->x0,&A->x1);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x0,&B->x1);
    fp2_mul_nonmod_montgomery(&tmp5_fpd2,&tmp5_fp2,&tmp4_fp2);//tmp5

    fp2_add_nonmod_single(&tmp6_fp2,&A->x1,&A->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&B->x1,&B->x2);
    fp2_mul_nonmod_montgomery(&tmp6_fpd2,&tmp6_fp2,&tmp4_fp2);//tmp6
    
    fp2_add_nonmod_single(&tmp7_fp2,&B->x0,&B->x2);
    fp2_add_nonmod_single(&tmp4_fp2,&A->x0,&A->x2);
    fp2_mul_nonmod_montgomery(&tmp7_fpd2,&tmp7_fp2,&tmp4_fp2);//tmp7
    
    //x0
    fp2_sub_nonmod_double(&tmp6_fpd2,&tmp6_fpd2,&tmp2_fpd2);
    fp2_sub_nonmod_double(&tmp6_fpd2,&tmp6_fpd2,&tmp3_fpd2);
    fp2_mul_basis_lazy_double(&tmp4_fpd2,&tmp6_fpd2);
    	
    fp2_add_nonmod_double(&ANS->x0,&tmp1_fpd2,&tmp4_fpd2);
    
    //x1
    fp2_sub_nonmod_double(&tmp5_fpd2,&tmp5_fpd2,&tmp1_fpd2);
    fp2_sub_nonmod_double(&tmp5_fpd2,&tmp5_fpd2,&tmp2_fpd2);
    fp2_mul_basis_lazy_double(&tmp4_fpd2,&tmp3_fpd2);
    fp2_add_nonmod_double(&ANS->x1,&tmp4_fpd2,&tmp5_fpd2);
    
    //x2
    fp2_sub_nonmod_double(&tmp7_fpd2,&tmp7_fpd2,&tmp1_fpd2);
    fp2_sub_nonmod_double(&tmp7_fpd2,&tmp7_fpd2,&tmp3_fpd2);
    fp2_add_nonmod_double(&ANS->x2,&tmp2_fpd2,&tmp7_fpd2);
    
}
void fp6_mul_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI){
    fp2_mul_ui(&ANS->x0,&A->x0,UI);
    fp2_mul_ui(&ANS->x1,&A->x1,UI);
    fp2_mul_ui(&ANS->x2,&A->x2,UI);
}

void fp6_mul_mpn(fp6_t *ANS,fp6_t *A,mp_limb_t *B){
    fp2_mul_mpn(&ANS->x0,&A->x0,B);
    fp2_mul_mpn(&ANS->x1,&A->x1,B);
    fp2_mul_mpn(&ANS->x2,&A->x2,B);
}

void fp6_mul_basis(fp6_t *ANS,fp6_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    fp2_set(&tmp1_fp2,&A->x0);
    fp2_set(&tmp2_fp2,&A->x1);
    fp2_set(&tmp3_fp2,&A->x2);
	
    fp_sub(&ANS->x0.x0,&tmp3_fp2.x0,&tmp3_fp2.x1);
    fp_add(&ANS->x0.x1,&tmp3_fp2.x0,&tmp3_fp2.x1);
    fp_set(&ANS->x1.x0,&tmp1_fp2.x0);
    fp_set(&ANS->x1.x1,&tmp1_fp2.x1);
    fp_set(&ANS->x2.x0,&tmp2_fp2.x0);
    fp_set(&ANS->x2.x1,&tmp2_fp2.x1);
}
void fp6_mul_basis_double(fpd6_t *ANS,fpd6_t *A){
    static fpd2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    fpd2_set(&tmp1_fp2,&A->x0);
    fpd2_set(&tmp2_fp2,&A->x1);
    fpd2_set(&tmp3_fp2,&A->x2);
	
    fp_sub_nonmod_double(&ANS->x0.x0,&tmp3_fp2.x0,&tmp3_fp2.x1);
    fp_add_nonmod_double(&ANS->x0.x1,&tmp3_fp2.x0,&tmp3_fp2.x1);
    fpd_set(&ANS->x1.x0,&tmp1_fp2.x0);
    fpd_set(&ANS->x1.x1,&tmp1_fp2.x1);
    fpd_set(&ANS->x2.x0,&tmp2_fp2.x0);
    fpd_set(&ANS->x2.x1,&tmp2_fp2.x1);
}
void fp6_mul_basis_lazy(fp6_t *ANS,fp6_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    fp2_set(&tmp1_fp2,&A->x0);
    fp2_set(&tmp2_fp2,&A->x1);
    fp2_set(&tmp3_fp2,&A->x2);
	
    fp_sub(&ANS->x0.x0,&tmp3_fp2.x0,&tmp3_fp2.x1);
    fp_add(&ANS->x0.x1,&tmp3_fp2.x0,&tmp3_fp2.x1);
    fp_set(&ANS->x1.x0,&tmp1_fp2.x0);
    fp_set(&ANS->x1.x1,&tmp1_fp2.x1);
    fp_set(&ANS->x2.x0,&tmp2_fp2.x0);
    fp_set(&ANS->x2.x1,&tmp2_fp2.x1);
}
void fp6_sqr(fp6_t *ANS,fp6_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    fp2_sqr(&tmp1_fp2,&A->x0);        //x0^2
    fp2_sqr(&tmp4_fp2,&A->x2);        //x2^2
    fp2_add(&tmp5_fp2,&A->x1,&A->x1);        //2x1
    fp2_mul(&tmp2_fp2,&tmp5_fp2,&A->x2);  //2x1x2
    fp2_mul(&tmp3_fp2,&A->x0,&tmp5_fp2);  //2x0x1
    fp2_add(&tmp5_fp2,&A->x0,&A->x1);        //x0+x1+x2
    fp2_add(&tmp5_fp2,&tmp5_fp2,&A->x2);
    
    //x0
    fp2_mul_basis(&ANS->x0,&tmp2_fp2);
    fp2_add(&ANS->x0,&ANS->x0,&tmp1_fp2);
    //x1
    fp2_mul_basis(&ANS->x1,&tmp4_fp2);
    fp2_add(&ANS->x1,&ANS->x1,&tmp3_fp2);
    //x2
    fp2_sqr(&ANS->x2,&tmp5_fp2);
    fp2_add(&tmp5_fp2,&tmp1_fp2,&tmp4_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp2_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp3_fp2);
    fp2_sub(&ANS->x2,&ANS->x2,&tmp5_fp2);
}
void fp6_sqr_lazy(fp6_t *ANS,fp6_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    fp2_sqr_lazy(&tmp1_fp2,&A->x0);        //x0^2
    fp2_sqr_lazy(&tmp4_fp2,&A->x2);        //x2^2
    fp2_add_nonmod_single(&tmp5_fp2,&A->x1,&A->x1);        //2x1

    fp2_mul_lazy(&tmp2_fp2,&tmp5_fp2,&A->x2);  //2x1x2
    fp2_mul_lazy(&tmp3_fp2,&A->x0,&tmp5_fp2);  //2x0x1
    fp2_add(&tmp5_fp2,&A->x0,&A->x1);        //x0+x1+x2
    fp2_add(&tmp5_fp2,&tmp5_fp2,&A->x2);
    
    //x0
    fp2_mul_basis(&ANS->x0,&tmp2_fp2);
    fp2_add(&ANS->x0,&ANS->x0,&tmp1_fp2);

    //x1
    fp2_mul_basis(&ANS->x1,&tmp4_fp2);
    fp2_add(&ANS->x1,&ANS->x1,&tmp3_fp2);

    //x2
    fp2_sqr_lazy(&ANS->x2,&tmp5_fp2);
    fp2_add(&tmp5_fp2,&tmp1_fp2,&tmp4_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp2_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp3_fp2);
    fp2_sub(&ANS->x2,&ANS->x2,&tmp5_fp2);
}
void fp6_sqr_lazy_montgomery(fp6_t *ANS,fp6_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    fp2_sqr_lazy(&tmp1_fp2,&A->x0);        //x0^2
    fp2_sqr_lazy(&tmp4_fp2,&A->x2);        //x2^2
    fp2_add_nonmod_single(&tmp5_fp2,&A->x1,&A->x1);        //2x1

    fp2_mul_lazy_montgomery(&tmp2_fp2,&tmp5_fp2,&A->x2);  //2x1x2
    fp2_mul_lazy_montgomery(&tmp3_fp2,&A->x0,&tmp5_fp2);  //2x0x1
    fp2_add(&tmp5_fp2,&A->x0,&A->x1);        //x0+x1+x2
    fp2_add(&tmp5_fp2,&tmp5_fp2,&A->x2);
    
    //x0
    fp2_mul_basis(&ANS->x0,&tmp2_fp2);
    fp2_add(&ANS->x0,&ANS->x0,&tmp1_fp2);

    //x1
    fp2_mul_basis(&ANS->x1,&tmp4_fp2);
    fp2_add(&ANS->x1,&ANS->x1,&tmp3_fp2);

    //x2
    fp2_sqr_lazy(&ANS->x2,&tmp5_fp2);
    fp2_add(&tmp5_fp2,&tmp1_fp2,&tmp4_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp2_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp3_fp2);
    fp2_sub(&ANS->x2,&ANS->x2,&tmp5_fp2);
}
void fp6_add(fp6_t *ANS,fp6_t *A,fp6_t *B){
    fp2_add(&ANS->x0,&A->x0,&B->x0);
    fp2_add(&ANS->x1,&A->x1,&B->x1);
    fp2_add(&ANS->x2,&A->x2,&B->x2);
}
void fp6_add_nonmod_single(fp6_t *ANS,fp6_t *A,fp6_t *B){
    fp2_add_nonmod_single(&ANS->x0,&A->x0,&B->x0);
    fp2_add_nonmod_single(&ANS->x1,&A->x1,&B->x1);
    fp2_add_nonmod_single(&ANS->x2,&A->x2,&B->x2);
}
void fp6_add_nonmod_double(fpd6_t *ANS,fpd6_t *A,fpd6_t *B){
    fp2_add_nonmod_double(&ANS->x0,&A->x0,&B->x0);
    fp2_add_nonmod_double(&ANS->x1,&A->x1,&B->x1);
    fp2_add_nonmod_double(&ANS->x2,&A->x2,&B->x2);
}
void fp6_add_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI){
    fp2_add_ui(&ANS->x0,&A->x0,UI);
    fp2_add_ui(&ANS->x1,&A->x1,0);
    fp2_add_ui(&ANS->x2,&A->x2,0);
}

void fp6_add_ui_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI){
    fp2_add_ui_ui(&ANS->x0,&A->x0,UI);
    fp2_add_ui_ui(&ANS->x1,&A->x1,UI);
    fp2_add_ui_ui(&ANS->x2,&A->x2,UI);
}
void fp6_add_mpn(fp6_t *ANS,fp6_t *A,mp_limb_t *B){
    fp2_add_mpn(&ANS->x0,&A->x0,B);
    fp2_add_mpn(&ANS->x1,&A->x1,B);
    fp2_add_mpn(&ANS->x2,&A->x2,B);
}
void fp6_sub(fp6_t *ANS,fp6_t *A,fp6_t *B){
    fp2_sub(&ANS->x0,&A->x0,&B->x0);
    fp2_sub(&ANS->x1,&A->x1,&B->x1);
    fp2_sub(&ANS->x2,&A->x2,&B->x2);
}
void fp6_sub_nonmod_single(fp6_t *ANS,fp6_t *A,fp6_t *B){
    fp2_sub_nonmod_single(&ANS->x0,&A->x0,&B->x0);
    fp2_sub_nonmod_single(&ANS->x1,&A->x1,&B->x1);
    fp2_sub_nonmod_single(&ANS->x2,&A->x2,&B->x2);
}
void fp6_sub_nonmod_double(fpd6_t *ANS,fpd6_t *A,fpd6_t *B){
    fp2_sub_nonmod_double(&ANS->x0,&A->x0,&B->x0);
    fp2_sub_nonmod_double(&ANS->x1,&A->x1,&B->x1);
    fp2_sub_nonmod_double(&ANS->x2,&A->x2,&B->x2);
}
void fp6_sub_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI){
    fp2_sub_ui(&ANS->x0,&A->x0,UI);
    fp2_sub_ui(&ANS->x1,&A->x1,0);
    fp2_sub_ui(&ANS->x2,&A->x2,0);
}

void fp6_sub_ui_ui(fp6_t *ANS,fp6_t *A,unsigned long int UI){
    fp2_sub_ui_ui(&ANS->x0,&A->x0,UI);
    fp2_sub_ui_ui(&ANS->x1,&A->x1,UI);
    fp2_sub_ui_ui(&ANS->x2,&A->x2,UI);
}
void fp6_sub_mpn(fp6_t *ANS,fp6_t *A,mp_limb_t *B){
    fp2_sub_mpn(&ANS->x0,&A->x0,B);
    fp2_sub_mpn(&ANS->x1,&A->x1,B);
    fp2_sub_mpn(&ANS->x2,&A->x2,B);
}

void fp6_inv(fp6_t *ANS,fp6_t *A){
	static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2,tmp8_fp2;

    fp2_sqr(&tmp1_fp2,&A->x0);
    fp2_sqr(&tmp2_fp2,&A->x1);
    fp2_sqr(&tmp3_fp2,&A->x2);
    
    fp2_mul(&tmp4_fp2,&A->x1,&A->x2);
    fp2_mul_basis(&tmp4_fp2,&tmp4_fp2);
    fp2_sub(&tmp6_fp2,&tmp1_fp2,&tmp4_fp2);
    
    fp2_mul(&tmp4_fp2,&A->x0,&A->x1);
    fp2_mul_basis(&tmp7_fp2,&tmp3_fp2);
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp4_fp2);
    
    fp2_mul(&tmp4_fp2,&A->x0,&A->x2);
    fp2_sub(&tmp8_fp2,&tmp2_fp2,&tmp4_fp2);
    
    fp2_mul(&tmp1_fp2,&tmp1_fp2,&A->x0);
    fp2_mul(&tmp3_fp2,&tmp3_fp2,&A->x2);
    fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);
    
    fp2_add(&tmp5_fp2,&tmp4_fp2,&tmp4_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp4_fp2);
    fp2_sub(&tmp5_fp2,&tmp2_fp2,&tmp5_fp2);
    fp2_mul(&tmp5_fp2,&tmp5_fp2,&A->x1);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp3_fp2);
    fp2_mul_basis(&tmp5_fp2,&tmp5_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp1_fp2);
    
    fp2_inv(&tmp5_fp2,&tmp5_fp2);
    
    fp2_mul(&ANS->x0,&tmp6_fp2,&tmp5_fp2);
    fp2_mul(&ANS->x1,&tmp7_fp2,&tmp5_fp2);
    fp2_mul(&ANS->x2,&tmp8_fp2,&tmp5_fp2);
}
void fp6_inv_lazy(fp6_t *ANS,fp6_t *A){
	static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2,tmp8_fp2;

    fp2_sqr_lazy(&tmp1_fp2,&A->x0);
    fp2_sqr_lazy(&tmp2_fp2,&A->x1);
    fp2_sqr_lazy(&tmp3_fp2,&A->x2);
    
    fp2_mul_lazy(&tmp4_fp2,&A->x1,&A->x2);
    fp2_mul_basis(&tmp4_fp2,&tmp4_fp2);
    fp2_sub(&tmp6_fp2,&tmp1_fp2,&tmp4_fp2);//tmp6
    
    fp2_mul_lazy(&tmp4_fp2,&A->x0,&A->x1);
    fp2_mul_basis(&tmp7_fp2,&tmp3_fp2);
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp4_fp2);//tmp7
    
    fp2_mul_lazy(&tmp4_fp2,&A->x0,&A->x2);//tmp4
    fp2_sub(&tmp8_fp2,&tmp2_fp2,&tmp4_fp2);//tmp8
    
    fp2_mul_lazy(&tmp1_fp2,&tmp1_fp2,&A->x0);//tmp1
    fp2_mul_lazy(&tmp3_fp2,&tmp3_fp2,&A->x2);
    fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);//tmp3
    
    fp2_add_nonmod_single(&tmp5_fp2,&tmp4_fp2,&tmp4_fp2);
    fp2_add_nonmod_single(&tmp5_fp2,&tmp5_fp2,&tmp4_fp2);
    fp2_sub_nonmod_single(&tmp5_fp2,&tmp2_fp2,&tmp5_fp2);
    fp2_mul_lazy(&tmp5_fp2,&tmp5_fp2,&A->x1);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp3_fp2);
    fp2_mul_basis(&tmp5_fp2,&tmp5_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp1_fp2);//mod
    
    fp2_inv_lazy(&tmp5_fp2,&tmp5_fp2);
    
    fp2_mul_lazy(&ANS->x0,&tmp6_fp2,&tmp5_fp2);
    fp2_mul_lazy(&ANS->x1,&tmp7_fp2,&tmp5_fp2);
    fp2_mul_lazy(&ANS->x2,&tmp8_fp2,&tmp5_fp2);
}
void fp6_inv_lazy_montgomery(fp6_t *ANS,fp6_t *A){
	static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,tmp6_fp2,tmp7_fp2,tmp8_fp2;

    fp2_sqr_lazy_montgomery(&tmp1_fp2,&A->x0);
    fp2_sqr_lazy_montgomery(&tmp2_fp2,&A->x1);
    fp2_sqr_lazy_montgomery(&tmp3_fp2,&A->x2);
    
    fp2_mul_lazy_montgomery(&tmp4_fp2,&A->x1,&A->x2);
    fp2_mul_basis(&tmp4_fp2,&tmp4_fp2);
    fp2_sub(&tmp6_fp2,&tmp1_fp2,&tmp4_fp2);//tmp6
    
    fp2_mul_lazy_montgomery(&tmp4_fp2,&A->x0,&A->x1);
    fp2_mul_basis(&tmp7_fp2,&tmp3_fp2);
    fp2_sub(&tmp7_fp2,&tmp7_fp2,&tmp4_fp2);//tmp7
    
    fp2_mul_lazy_montgomery(&tmp4_fp2,&A->x0,&A->x2);//tmp4
    fp2_sub(&tmp8_fp2,&tmp2_fp2,&tmp4_fp2);//tmp8
    
    fp2_mul_lazy_montgomery(&tmp1_fp2,&tmp1_fp2,&A->x0);//tmp1
    fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp3_fp2,&A->x2);
    fp2_mul_basis(&tmp3_fp2,&tmp3_fp2);//tmp3
    
    fp2_add_nonmod_single(&tmp5_fp2,&tmp4_fp2,&tmp4_fp2);
    fp2_add_nonmod_single(&tmp5_fp2,&tmp5_fp2,&tmp4_fp2);
    fp2_sub_nonmod_single(&tmp5_fp2,&tmp2_fp2,&tmp5_fp2);
    fp2_mul_lazy_montgomery(&tmp5_fp2,&tmp5_fp2,&A->x1);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp3_fp2);
    fp2_mul_basis(&tmp5_fp2,&tmp5_fp2);
    fp2_add(&tmp5_fp2,&tmp5_fp2,&tmp1_fp2);//mod
    
    fp2_inv_lazy_montgomery(&tmp5_fp2,&tmp5_fp2);
    
    fp2_mul_lazy_montgomery(&ANS->x0,&tmp6_fp2,&tmp5_fp2);
    fp2_mul_lazy_montgomery(&ANS->x1,&tmp7_fp2,&tmp5_fp2);
    fp2_mul_lazy_montgomery(&ANS->x2,&tmp8_fp2,&tmp5_fp2);
}
int  fp6_legendre(fp6_t *A){
    mpz_t exp;
    mpz_init(exp);
    fp6_t tmp;
    fp6_init(&tmp);
    
    mpz_pow_ui(exp,prime_z,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp6_pow(&tmp,A,exp);
    
    if(fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

int  fp6_isCNR(fp6_t *A){
    fp6_t tmp;
    fp6_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime_z,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    fp6_pow(&tmp,A,exp);
    
    if(fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

void fp6_sqrt(fp6_t *ANS,fp6_t *A){
    fp6_t tmp1,tmp2;
    fp6_init(&tmp1);
    fp6_init(&tmp2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    fp2_set(&tmp1.x0,&A->x0);
    fp2_mul_mpn(&tmp1.x1,&A->x1,frobenius_constant[f_p4][1].x0.x0);
    fp2_mul_mpn(&tmp1.x2,&A->x2,frobenius_constant[f_p4][2].x0.x0);
    
    fp2_set(&tmp2.x0,&A->x0);
    fp2_mul_mpn(&tmp2.x1,&A->x1,frobenius_constant[f_p2][1].x0.x0);
    fp2_mul_mpn(&tmp2.x2,&A->x2,frobenius_constant[f_p2][2].x0.x0);
    
    fp6_mul(&tmp1,&tmp1,&tmp2);
    fp6_mul(&tmp1,&tmp1,A);
    fp6_set_ui(&tmp2,0);
    fp2_sqrt(&tmp2.x0,&tmp1.x0);
    fp2_inv(&tmp2.x0,&tmp2.x0);
    fp2_set(&tmp2.x0,&tmp2.x0);
    mpz_pow_ui(exp,prime_z,8);
    mpz_pow_ui(buf,prime_z,4);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    fp6_pow(&tmp1,A,exp);
    fp6_mul(&tmp1,&tmp1,&tmp2);
    fp6_set(ANS,&tmp1);
    
    mpz_clear(exp);
    mpz_clear(buf);
}

void fp6_pow(fp6_t *ANS,fp6_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp6_t tmp;
    fp6_init(&tmp);
    fp6_set(&tmp,A);
    
    for(i=1;i<length; i++){
        fp6_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            fp6_mul(&tmp,A,&tmp);
        }
    }
    
    fp6_set(ANS,&tmp);
}

int  fp6_cmp(fp6_t *A,fp6_t *B){
    if(fp2_cmp(&A->x0,&B->x0)==0 && fp2_cmp(&A->x1,&B->x1)==0 && fp2_cmp(&A->x2,&B->x2)==0){
        return 0;   
    }
    return 1;
}

int  fp6_cmp_ui(fp6_t *A,unsigned long int UI){
    if(fp2_cmp_ui(&A->x0,UI)==0 && fp2_cmp_ui(&A->x1,UI)==0 && fp2_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  fp6_cmp_mpn(fp6_t *A,mp_limb_t *B){
    if(fp2_cmp_mpn(&A->x0,B)==0 && fp2_cmp_mpn(&A->x1,B)==0 && fp2_cmp_mpn(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  fp6_cmp_zero(fp6_t *A){
    if(fp2_cmp_zero(&A->x0)==0 && fp2_cmp_zero(&A->x1)==0 && fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

int  fp6_cmp_one(fp6_t *A){
    if(fp2_cmp_one(&A->x0)==0 && fp2_cmp_zero(&A->x1)==0 && fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}
