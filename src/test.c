#include <ELiPS/test.h>
/*----------------------------------------------------------------------------*/
//test
int test_Field(int fp2,int fp6,int fp12,int sqr){
    int i,n=0;
    float mul_time=0,mul_lazy_time=0;
    float sqr_time=0,sqr_lazy_time=0;
    float com_time=0,rec_time=0,gs_time=0;
    float com_lazy_time=0,rec_lazy_time=0,gs_lazy_time=0;
    EFp12 P,Q;
    Fp12 tmp,t0,t1,f;
    struct timeval tv_A,tv_B;
    
    printf("====================================================================================\n");
    printf("Field test\n");
    Fp2 A_Fp2,B_Fp2,test1,test2;
    Fp6 A_Fp6,B_Fp6,test1_Fp6,test2_Fp6,test3_Fp6;
    Fp12 A_Fp12,B_Fp12,test1_Fp12,test2_Fp12,test3_Fp12,test4_Fp12;
    Fp12 test0_Fp12,test1l_Fp12,test2l_Fp12,test3l_Fp12;

    Fp2_init(&A_Fp2);
    Fp2_init(&B_Fp2);
    Fp2_init(&test1);
    Fp2_init(&test2);

    
    Fp6_init(&A_Fp6);
    Fp6_init(&B_Fp6);
    Fp6_init(&test1_Fp6);
    Fp6_init(&test2_Fp6);
    Fp6_init(&test3_Fp6);
    
    Fp12_init(&A_Fp12);
    Fp12_init(&B_Fp12);
    Fp12_init(&test0_Fp12);
    Fp12_init(&test1_Fp12);
    Fp12_init(&test2_Fp12);
    Fp12_init(&test3_Fp12);

    Fp12_init(&test1l_Fp12);
    Fp12_init(&test2l_Fp12);
    Fp12_init(&test3l_Fp12);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_mul test\n");
mul_time=0,mul_lazy_time=0;
for(i=0;i<fp2;i++){
    Fp2_set_random(&A_Fp2,state);
    Fp2_set_random(&B_Fp2,state);

    gettimeofday(&tv_A,NULL);
    Fp2_mul(&test1,&A_Fp2,&B_Fp2);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp2_mul_lazy(&test2,&A_Fp2,&B_Fp2);
    gettimeofday(&tv_B,NULL);
    mul_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp2_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	Fp2_printf(&test1,"");
	Fp2_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("Fp2 mul.      : %.4f[ms]\n",mul_time/fp2);
    printf("Fp2 mul lazy. : %.4f[ms]\n",mul_lazy_time/fp2);


printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_sqr test\n");
sqr_time=0,sqr_lazy_time=0;
for(i=0;i<fp2;i++){
    Fp2_set_random(&A_Fp2,state);

    gettimeofday(&tv_A,NULL);
    Fp2_sqr(&test1,&A_Fp2);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp2_sqr_lazy(&test2,&A_Fp2);
    gettimeofday(&tv_B,NULL);
    sqr_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp2_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	Fp2_printf(&test1,"");
	Fp2_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("Fp2 sqr.      : %.4f[ms]\n",sqr_time/fp2);
    printf("Fp2 sqr lazy. : %.4f[ms]\n",sqr_lazy_time/fp2);


printf("------------------------------------------------------------------------------------\n");
    printf("Fp6_mul test\n");
mul_time=0,mul_lazy_time=0;
for(i=0;i<fp6;i++){

    Fp6_set_random(&A_Fp6,state);
    Fp6_set_random(&B_Fp6,state);

    gettimeofday(&tv_A,NULL);
    Fp6_mul(&test1_Fp6,&A_Fp6,&B_Fp6);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp6_mul_lazy(&test2_Fp6,&A_Fp6,&B_Fp6);
    gettimeofday(&tv_B,NULL);
    mul_lazy_time+=timedifference_msec(tv_A,tv_B);



    if(Fp6_cmp(&test1_Fp6,&test2_Fp6)!=0){
        printf("failed!\n\n");
	Fp6_printf(&test1_Fp6,"");
	Fp6_printf(&test2_Fp6,"\n");
	printf("\n\n");
	return 1;
    }
}
    
    printf("Fp6 mul.       : %.4f[ms]\n",mul_time/fp6);
    printf("Fp6 mul lazy.  : %.4f[ms]\n",mul_lazy_time/fp6);

printf("------------------------------------------------------------------------------------\n");
    printf("Fp6_sqr test\n");
sqr_time=0,sqr_lazy_time=0;
for(i=0;i<fp6;i++){
    Fp6_set_random(&A_Fp6,state);

    gettimeofday(&tv_A,NULL);
    Fp6_sqr(&test1_Fp6,&A_Fp6);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp6_sqr_lazy(&test2_Fp6,&A_Fp6);
    gettimeofday(&tv_B,NULL);
    sqr_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp6_cmp(&test1_Fp6,&test2_Fp6)!=0){
        printf("failed!\n\n");
	Fp6_printf(&test1_Fp6,"");
	Fp6_printf(&test2_Fp6,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("Fp6 sqr.      : %.4f[ms]\n",sqr_time/fp6);
    printf("Fp6 sqr lazy. : %.4f[ms]\n",sqr_lazy_time/fp6);

printf("------------------------------------------------------------------------------------\n");
    printf("Fp12_mul test\n");
mul_time=0,mul_lazy_time=0;
for(i=0;i<fp12;i++){
    Fp12_set_random(&A_Fp12,state);
    Fp12_set_random(&B_Fp12,state);

    gettimeofday(&tv_A,NULL);
    Fp12_mul(&test3_Fp12,&A_Fp12,&B_Fp12);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp12_mul_lazy(&test4_Fp12,&A_Fp12,&B_Fp12);
    gettimeofday(&tv_B,NULL);
    mul_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(Fp12_cmp(&test3_Fp12,&test4_Fp12)!=0){
        printf("failed!\n\n");
	Fp12_printf(&test3_Fp12,"");
	Fp12_printf(&test4_Fp12,"\n");
	printf("\n\n");
	return 1;
    }
}
    
    printf("Fp12 mul.        : %.4f[ms]\n",mul_time/fp12);
    printf("Fp12 mul lazy.   : %.4f[ms]\n",mul_lazy_time/fp12);

printf("------------------------------------------------------------------------------------\n");
    printf("Fp12_sqr test\n");
sqr_time=0,sqr_lazy_time=0;
for(i=0;i<fp12;i++){
    Fp12_set_random(&A_Fp12,state);

    gettimeofday(&tv_A,NULL);
    Fp12_sqr(&test1_Fp12,&A_Fp12);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp12_sqr_lazy(&test2_Fp12,&A_Fp12);
    gettimeofday(&tv_B,NULL);
    sqr_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp12_cmp(&test1_Fp12,&test2_Fp12)!=0){
        printf("failed!\n\n");
	Fp12_printf(&test1_Fp12,"");
	Fp12_printf(&test2_Fp12,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("Fp12 sqr.      : %.4f[ms]\n",sqr_time/fp12);
    printf("Fp12 sqr lazy. : %.4f[ms]\n",sqr_lazy_time/fp12);

printf("------------------------------------------------------------------------------------\n");
    printf("Fp12_sqr_cyclotomic test\n");
sqr_time=0,sqr_lazy_time=0;
for(i=0;i<fp12;i++){
    Fp12_set_random(&A_Fp12,state);

    gettimeofday(&tv_A,NULL);
    Fp12_sqr_cyclotomic(&test1_Fp12,&A_Fp12);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp12_sqr_cyclotomic_lazy(&test2_Fp12,&A_Fp12);
    gettimeofday(&tv_B,NULL);
    sqr_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp12_cmp(&test1_Fp12,&test2_Fp12)!=0){
        printf("failed!\n\n");
	Fp12_printf(&test1_Fp12,"");
	Fp12_printf(&test2_Fp12,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("Fp12 sqr_cyclotomic.      : %.4f[ms]\n",sqr_time/fp12);
    printf("Fp12 sqr_cyclotomic lazy. : %.4f[ms]\n",sqr_lazy_time/fp12);

printf("------------------------------------------------------------------------------------\n");
    printf("Fp12_sqr_compressed test\n");

for(i=0;i<sqr;i++){
    BLS12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);
    BLS12_Miller_algo_for_opt_ate(&f,&P,&Q);
    
    //EASY PART
    Fp12_frobenius_map_p6(&t0,&f);//f^(p^6)
    Fp12_inv(&t1,&f);//f^-1
    Fp12_mul(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    Fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    Fp12_mul(&tmp,&t0,&tmp);//f^(p^2)*f
    
    
    Fp12_sqr_cyclotomic(&test0_Fp12,&tmp);
    
    gettimeofday(&tv_A,NULL);
    Fp12_sqr_compressed(&test1_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    com_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    Fp12_sqr_compressed_lazy(&test1l_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    com_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    Fp12_sqr_recover_g1(&test2_Fp12,&test1_Fp12);
    Fp12_sqr_recover_g0(&test2_Fp12,&test2_Fp12);
    gettimeofday(&tv_B,NULL);
    rec_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    Fp12_sqr_recover_g1_lazy(&test2l_Fp12,&test1l_Fp12);
    Fp12_sqr_recover_g0_lazy(&test2l_Fp12,&test2l_Fp12);
    gettimeofday(&tv_B,NULL);
    rec_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    Fp12_sqr_GS(&test3_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    gs_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    Fp12_sqr_GS_lazy(&test3l_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    gs_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    if(Fp12_cmp(&test0_Fp12,&test2_Fp12)!=0 || Fp12_cmp(&test1_Fp12,&test1l_Fp12)!=0 || Fp12_cmp(&test2_Fp12,&test2l_Fp12)!=0 || Fp12_cmp(&test0_Fp12,&test3_Fp12)!=0 || Fp12_cmp(&test3_Fp12,&test3l_Fp12)!=0){
        printf("failed!\n\n");
	    return 1;
    }
}
    
    printf("Sqr compressed.             : %.4f[ms]\n",com_time/sqr);
    printf("Sqr compressed lazy.        : %.4f[ms]\n",com_lazy_time/sqr);
    printf("Sqr compress_recover.       : %.4f[ms]\n",rec_time/sqr);
    printf("Sqr compress_recover lazy.  : %.4f[ms]\n",rec_lazy_time/sqr);
    printf("Sqr Karabina.               : %.4f[ms]\n",(com_time+rec_time)/sqr);
    printf("Sqr Karabina lazy.          : %.4f[ms]\n",(com_lazy_time+rec_lazy_time)/sqr);
    printf("Sqr GS.                     : %.4f[ms]\n",gs_time/sqr);
    printf("Sqr GS lazy.                : %.4f[ms]\n",gs_lazy_time/sqr);

return 0;
}

int test_EFp(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0,ecd_Jacobian_time=0,ecd_Jacobian_lazy_time=0;
    float eca_time=0,eca_lazy_time=0,eca_Jacobian_time=0,eca_Jacobian_lazy_time=0;
    float scm_time=0,scm_lazy_time=0,scm_Jacobian_time=0,scm_Jacobian_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp A_EFp,B_EFp,test0,test1,test2,test3;
    EFpZ A_EFpZ,B_EFpZ,testZ1,testZ2,testZ3;
    EFp_init(&A_EFp);
    EFp_init(&B_EFp);
    EFp_init(&test0);
    EFp_init(&test1);
    EFp_init(&test2);
    EFp_init(&test3);
    EFpZ_init(&A_EFpZ);
    EFpZ_init(&B_EFpZ);
    EFpZ_init(&testZ1);
    EFpZ_init(&testZ2);
    EFpZ_init(&testZ3);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


printf("------------------------------------------------------------------------------------\n");
    printf("EFp_ECD test\n");

for(i=0;i<ecd;i++){

    EFp_rational_point(&B_EFp);
    EFp_to_EFpZ(&B_EFpZ,&B_EFp);

    gettimeofday(&tv_A,NULL);
    EFp_ECD(&test0,&B_EFp);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_ECD_lazy(&test1,&B_EFp);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp_ECD_Jacobian(&testZ2,&B_EFpZ);
    EFp_Jacobian(&test2,&testZ2);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_ECD_Jacobian_lazy(&testZ3,&B_EFpZ);
    EFp_Jacobian(&test3,&testZ3);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(EFp_cmp(&test0,&test1)!=0 || EFp_cmp(&test1,&test2)!=0 || EFp_cmp(&test1,&test3)!=0){
        printf("failed!\n\n");
	    EFp_printf(&test1,"");
	    EFp_printf(&test2,"\n");
	    EFp_printf(&test3,"\n");
        EFpZ_printf(&testZ2,"\ntestZ2=");printf("\n");
        EFpZ_printf(&testZ3,"\ntestZ3=");printf("\n");
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp ECD.               : %.4f[ms]\n",ecd_time/ecd);
    printf("EFp ECD lazy.          : %.4f[ms]\n",ecd_lazy_time/ecd);
    printf("EFp ECD Jacobian.      : %.4f[ms]\n",ecd_Jacobian_time/ecd);
    printf("EFp ECD Jacobian lazy. : %.4f[ms]\n",ecd_Jacobian_lazy_time/ecd);


printf("------------------------------------------------------------------------------------\n");
    printf("EFp_ECA test\n");
    
    EFp_rational_point(&A_EFp);
    EFp_to_EFpZ(&A_EFpZ,&A_EFp);

for(i=0;i<eca;i++){

    EFp_rational_point(&B_EFp);
    EFp_to_EFpZ(&B_EFpZ,&B_EFp);

    gettimeofday(&tv_A,NULL);
    EFp_ECA(&test0,&A_EFp,&B_EFp);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp_ECA_lazy(&test1,&A_EFp,&B_EFp);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp_ECA_Jacobian(&testZ2,&A_EFpZ,&B_EFpZ);
    EFp_Jacobian(&test2,&testZ2);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_ECA_Jacobian_lazy(&testZ3,&A_EFpZ,&B_EFpZ);
    EFp_Jacobian(&test3,&testZ3);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(EFp_cmp(&test0,&test1)!=0 || EFp_cmp(&test1,&test2)!=0 || EFp_cmp(&test1,&test3)!=0){
        printf("failed!\n\n");
	    EFp_printf(&test1,"");
	    EFp_printf(&test2,"\n");
	    EFp_printf(&test3,"\n");
        EFpZ_printf(&testZ2,"\ntestZ2=");printf("\n");
        EFpZ_printf(&testZ3,"\ntestZ3=");printf("\n");
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp ECA.               : %.4f[ms]\n",eca_time/eca);
    printf("EFp ECA lazy.          : %.4f[ms]\n",eca_lazy_time/eca);
    printf("EFp ECA Jacobian.      : %.4f[ms]\n",eca_Jacobian_time/eca);
    printf("EFp ECA Jacobian lazy. : %.4f[ms]\n",eca_Jacobian_lazy_time/eca);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp_SCM test\n");

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    EFp_rational_point(&A_EFp);

    gettimeofday(&tv_A,NULL);
    EFp_SCM(&test0,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_SCM_lazy(&test1,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp_SCM_Jacobian(&test2,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_SCM_Jacobian_lazy(&test3,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp_cmp(&test0,&test1)!=0 || EFp_cmp(&test1,&test2)!=0 || EFp_cmp(&test1,&test3)!=0){
        printf("failed!\n\n");
	    EFp_printf(&test1,"");
	    EFp_printf(&test2,"\n");
	    EFp_printf(&test3,"\n");
	    printf("\n\n");
	    return 1;
    }
}

    printf("EFp SCM.               : %.4f[ms]\n",scm_time/scm);
    printf("EFp SCM lazy.          : %.4f[ms]\n",scm_lazy_time/scm);
    printf("EFp SCM Jacobian.      : %.4f[ms]\n",scm_Jacobian_time/scm);
    printf("EFp SCM Jacobian Lazy. : %.4f[ms]\n",scm_Jacobian_lazy_time/scm);

    return 0;
}
int test_EFp2(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0,ecd_Jacobian_time=0,ecd_Jacobian_lazy_time=0;
    float eca_time=0,eca_lazy_time=0,eca_Jacobian_time=0,eca_Jacobian_lazy_time=0;
    float scm_time=0,scm_lazy_time=0,scm_Jacobian_time=0,scm_Jacobian_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp2 A_EFp2,B_EFp2,test0,test1,test2,test3;
    EFpZ2 A_EFpZ2,B_EFpZ2,testZ1,testZ2,testZ3;
    EFp2_init(&A_EFp2);
    EFp2_init(&B_EFp2);
    EFp2_init(&test0);
    EFp2_init(&test1);
    EFp2_init(&test2);
    EFp2_init(&test3);
    EFpZ2_init(&A_EFpZ2);
    EFpZ2_init(&B_EFpZ2);
    EFpZ2_init(&testZ1);
    EFpZ2_init(&testZ2);
    EFpZ2_init(&testZ3);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_ECD test\n");

for(i=0;i<ecd;i++){

    EFp2_rational_point(&B_EFp2);
    EFp2_to_EFpZ2(&B_EFpZ2,&B_EFp2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD(&test0,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD_lazy(&test1,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECD_Jacobian(&testZ2,&B_EFpZ2);
    EFp2_Jacobian(&test2,&testZ2);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD_Jacobian_lazy(&testZ3,&B_EFpZ2);
    EFp2_Jacobian(&test3,&testZ3);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(EFp2_cmp(&test0,&test1)!=0 || EFp2_cmp(&test1,&test2)!=0 || EFp2_cmp(&test1,&test3)!=0){
        printf("failed!\n\n");
	    EFp2_printf(&test1,"");
	    EFp2_printf(&test2,"\n");
	    EFp2_printf(&test3,"\n");
        EFpZ2_printf(&testZ2,"\ntestZ2=");printf("\n");
        EFpZ2_printf(&testZ3,"\ntestZ3=");printf("\n");
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp2 ECD.               : %.4f[ms]\n",ecd_time/ecd);
    printf("EFp2 ECD lazy.          : %.4f[ms]\n",ecd_lazy_time/ecd);
    printf("EFp2 ECD Jacobian.      : %.4f[ms]\n",ecd_Jacobian_time/ecd);
    printf("EFp2 ECD Jacobian lazy. : %.4f[ms]\n",ecd_Jacobian_lazy_time/ecd);


printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_ECA test\n");
    
    EFp2_rational_point(&A_EFp2);
    EFp2_to_EFpZ2(&A_EFpZ2,&A_EFp2);

for(i=0;i<eca;i++){

    EFp2_rational_point(&B_EFp2);
    EFp2_to_EFpZ2(&B_EFpZ2,&B_EFp2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA(&test0,&A_EFp2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECA_lazy(&test1,&A_EFp2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Jacobian(&testZ2,&A_EFpZ2,&B_EFpZ2);
    EFp2_Jacobian(&test2,&testZ2);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Jacobian_lazy(&testZ3,&A_EFpZ2,&B_EFpZ2);
    EFp2_Jacobian(&test3,&testZ3);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(EFp2_cmp(&test0,&test1)!=0 || EFp2_cmp(&test1,&test2)!=0 || EFp2_cmp(&test1,&test3)!=0){
        printf("failed!\n\n");
	    EFp2_printf(&test1,"");
	    EFp2_printf(&test2,"\n");
	    EFp2_printf(&test3,"\n");
        EFpZ2_printf(&testZ2,"\ntestZ2=");printf("\n");
        EFpZ2_printf(&testZ3,"\ntestZ3=");printf("\n");
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp2 ECA.               : %.4f[ms]\n",eca_time/eca);
    printf("EFp2 ECA lazy.          : %.4f[ms]\n",eca_lazy_time/eca);
    printf("EFp2 ECA Jacobian.      : %.4f[ms]\n",eca_Jacobian_time/eca);
    printf("EFp2 ECA Jacobian lazy. : %.4f[ms]\n",eca_Jacobian_lazy_time/eca);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_SCM test\n");

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    EFp2_rational_point(&A_EFp2);

    gettimeofday(&tv_A,NULL);
    EFp2_SCM(&test0,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_SCM_lazy(&test1,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp2_SCM_Jacobian(&test2,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_SCM_Jacobian_lazy(&test3,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp2_cmp(&test0,&test1)!=0 || EFp2_cmp(&test1,&test2)!=0 || EFp2_cmp(&test1,&test3)!=0){
        printf("failed!\n\n");
	    EFp2_printf(&test1,"");
	    EFp2_printf(&test2,"\n");
	    EFp2_printf(&test3,"\n");
	    printf("\n\n");
	    return 1;
    }
}

    printf("EFp2 SCM.               : %.4f[ms]\n",scm_time/scm);
    printf("EFp2 SCM lazy.          : %.4f[ms]\n",scm_lazy_time/scm);
    printf("EFp2 SCM Jacobian.      : %.4f[ms]\n",scm_Jacobian_time/scm);
    printf("EFp2 SCM Jacobian Lazy. : %.4f[ms]\n",scm_Jacobian_lazy_time/scm);

    return 0;
}
int test_EFp12(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0;
    float eca_time=0,eca_lazy_time=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp12 A_EFp12,B_EFp12,test1,test2;
    EFp12_init(&A_EFp12);
    EFp12_init(&B_EFp12);
    EFp12_init(&test1);
    EFp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


printf("------------------------------------------------------------------------------------\n");
    printf("EFp12_ECD test\n");

for(i=0;i<ecd;i++){

    EFp12_rational_point(&B_EFp12);

    gettimeofday(&tv_A,NULL);
    EFp12_ECD(&test1,&B_EFp12);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp12_ECD_lazy(&test2,&B_EFp12);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(EFp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp12_printf(&test1,"");
	EFp12_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp12 ECD.      : %.4f[ms]\n",ecd_time/ecd);
    printf("EFp12 ECD lazy. : %.4f[ms]\n",ecd_lazy_time/ecd);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp12_ECA test\n");

    EFp12_rational_point(&A_EFp12);
for(i=0;i<eca;i++){

    EFp12_rational_point(&B_EFp12);

    gettimeofday(&tv_A,NULL);
    EFp12_ECA(&test1,&A_EFp12,&B_EFp12);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp12_ECA_lazy(&test2,&A_EFp12,&B_EFp12);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp12_printf(&test1,"");
	EFp12_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp12 ECA.      : %.4f[ms]\n",eca_time/eca);
    printf("EFp12 ECA lazy. : %.4f[ms]\n",eca_lazy_time/eca);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp12_SCM test\n");

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    EFp12_rational_point(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    EFp12_SCM(&test1,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp12_SCM_lazy(&test2,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp12_printf(&test1,"");
	EFp12_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp12 SCM.      : %.4f[ms]\n",scm_time/scm);
    printf("EFp12 SCM lazy. : %.4f[ms]\n",scm_lazy_time/scm);

	return 0;
}
void test_Frobenius_map(){
    printf("====================================================================================\n");
    Fp12 A_Fp12,test1,test2;
    Fp12_init(&A_Fp12);
    Fp12_init(&test1);
    Fp12_init(&test2);
    mpz_t exp;
    mpz_init(exp);
    
    gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp12_set_random(&A_Fp12,state);
    Fp12_printf(&A_Fp12,"");
    printf("\n\n");
    
    printf("frobenius\n");
    Fp12_frobenius_map_p10(&test1,&A_Fp12);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    mpz_pow_ui(exp,prime_z,10);
    Fp12_pow(&test2,&A_Fp12,exp);
    Fp12_printf(&test2,"");
    printf("\n");
    
    if(Fp12_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(exp);
}

void test_skew_frobenius_map(){
    printf("====================================================================================\n");
    EFp12 Q,test1,test2;
    EFp12_init(&Q);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp2 twisted_Q;
    EFp2_init(&twisted_Q);
    
    EFp12_generate_G2(&Q);
    EFp12_to_EFp2(&twisted_Q,&Q);
    
    Fp12_frobenius_map_p10(&test1.x,&Q.x);
    Fp12_frobenius_map_p10(&test1.y,&Q.y);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    EFp2_skew_frobenius_map_p10(&twisted_Q,&twisted_Q);
    EFp2_to_EFp12(&test2,&twisted_Q);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
}


void test_twist(){
    printf("====================================================================================\n");
    EFp12 Q,test1,test2;
    EFp12_init(&Q);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp2 twist_Q;
    EFp2_init(&twist_Q);
    
    
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    
    EFp12_to_EFp2(&twist_Q,&Q);
    EFp2_ECD(&twist_Q,&twist_Q);
    EFp2_to_EFp12(&test1,&twist_Q);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    EFp12_ECD(&test2,&Q);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
}
int test_mod(int mulmod,int mod){
    int i,n=0;
    float add_time=0,dbl_time=0,mul_time=0,sqr_time=0,mod_time=0,mod_mont_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("Field Mod test\n");
    mp_limb_t A[FPLIMB],B[FPLIMB],C[FPLIMB],C1[FPLIMB2];
    mp_limb_t At[FPLIMB2],Bt[FPLIMB2],Ct[FPLIMB2],test_mul[FPLIMB2];
    Fp test1,test2;

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<mulmod;i++){
    mpn_random(A,FPLIMB);
    mpn_random(B,FPLIMB);

    Lazy_mod(A,A,FPLIMB);
    Lazy_mod(B,B,FPLIMB);

    mpn_mul_n(At,A,A,FPLIMB);
    mpn_mul_n(Bt,B,B,FPLIMB);

    gettimeofday(&tv_A,NULL);
    mpn_add_n(test_mul,At,At,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    add_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    mpn_lshift(test_mul,At,FPLIMB2,1);
    gettimeofday(&tv_B,NULL);
    dbl_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    mpn_mul_n(test_mul,A,A,FPLIMB);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    mpn_sqr(test_mul,A,FPLIMB);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    mpn_mod(&test1,test_mul,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    mod_time+=timedifference_msec(tv_A,tv_B);

}
    printf("Fp add.            : %.6f[ms]\n",add_time/mulmod);
    printf("Fp dbl.            : %.6f[ms]\n",dbl_time/mulmod);
    printf("Fp mod.            : %.6f[ms]\n",mod_time/mulmod);
    printf("Fp sqr.            : %.6f[ms]\n",sqr_time/mulmod);
    printf("Fp mul.            : %.6f[ms]\n",mul_time/mulmod);
    printf("Fp add + mul * 2.  : %.6f[ms]\n",(add_time+(mul_time*2))/mulmod);
mod_time=0;
for(i=0;i<mod;i++){
    mpn_random(A,FPLIMB);
    mpn_random(B,FPLIMB);

    Lazy_mod(A,A,FPLIMB);
    Lazy_mod(B,B,FPLIMB);

    mpn_mul_n(C1,A,B,FPLIMB);

    Fp_mul_montgomery(At,FPLIMB2,A,FPLIMB);
    Fp_MR(A,At,FPLIMB2);

    Fp_mul_montgomery(Bt,FPLIMB2,B,FPLIMB);
    Fp_MR(B,Bt,FPLIMB2);

    mpn_mul_n(Ct,A,B,FPLIMB);
    Fp_MR(C,Ct,FPLIMB2);

    gettimeofday(&tv_A,NULL);
    mpn_mod(&test1,C1,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    mod_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp_MR(C,Ct,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    mod_mont_time+=timedifference_msec(tv_A,tv_B);

    Fp_MR(test2.x0,C,FPLIMB);

    if(Fp_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	Fp_printf(&test1,"");
	Fp_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("Fp mod.            : %.6f[ms]\n",mod_time/mod);
    printf("Fp mod montgomery. : %.6f[ms]\n",mod_mont_time/mod);

    return 0;
    
}

void test_All(){
	int test_point,test_G1,test_G2,test_G3,test_pairing;
	int Lazy_Field,Lazy_EFp,Lazy_EFp2,Lazy_EFp12,Lazy_G1,Lazy_G2,Lazy_G3,Lazy_pairing;
	int Jacobian_EFp,Jacobian_EFp2,Jacobian_G1;
	
printf("====================================================================================\n");
    printf("Normal test\n");
	
	test_point = BLS12_test_rational_point();
	test_G1 = BLS12_test_G1_SCM(1);
	test_G2 = BLS12_test_G2_SCM(1);
	test_G3 = BLS12_test_G3_EXP(1);
	test_pairing = BLS12_test_opt_ate_pairing();

printf("====================================================================================\n");
    printf("Lazy test\n");
	
	Lazy_EFp = test_EFp_Lazy(100,100,10);
	Lazy_EFp2 = test_EFp2_Lazy(10,10,10);
	Lazy_EFp12 = test_EFp12_Lazy(10,10,1);

printf("====================================================================================\n");
    printf("Jacobian test\n");

	
	Jacobian_EFp = test_EFp_Jacobian(100,100,10);
	Jacobian_EFp2 = test_EFp2_Jacobian(100,100,10);


printf("====================================================================================\n");
    printf("Test Result\n\n");

printf("------------------------------------------------------------------------------------\n");
    printf("Normal test\n");

	printf("Test Rational Point  :");
	if(test_point==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G1 test         :");
	if(test_G1==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G2 test         :");
	if(test_G2==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G3 test         :");
	if(test_G3==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test pairing test    :");
	if(test_pairing==0) printf("Success\n");
	else printf("Failed\n");

printf("------------------------------------------------------------------------------------\n");
    printf("Lazy Reduction test\n");

	printf("Test Field Lazy      :");
	if(Lazy_Field==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test EFp Lazy        :");
	if(Lazy_EFp==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test EFp2 Lazy       :");
	if(Lazy_EFp2==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test EFp12 Lazy      :");
	if(Lazy_EFp12==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G1 Lazy         :");
	if(Lazy_G1==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G2 Lazy         :");
	if(Lazy_G2==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G3 Lazy         :");
	if(Lazy_G3==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test pairing Lazy    :");
	if(Lazy_pairing==0) printf("Success\n");
	else printf("Failed\n");

printf("------------------------------------------------------------------------------------\n");
    printf("Jacobian test\n");

	printf("Test EFp Jacobian    :");
	if(Jacobian_EFp==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test EFp2 Jacobian   :");
	if(Jacobian_EFp2==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G1 Jacobian     :");
	if(Jacobian_G1==0) printf("Success\n");
	else printf("Failed\n");
}
int Fast_test_BLS12(int g1, int g2,int g3,int pairing){
    int i;
    float scm1_time=0;
    float scm2_time=0;
    float exp_time=0;
    float opt_time=0;
    struct timeval tv_A,tv_B;
    
    EFp12 A_EFp12,test_EFp12;
    EFp12_init(&A_EFp12);
    EFp12_init(&test_EFp12);
    
    Fp12 A_Fp12,B_Fp12,test_Fp12;
    Fp12_init(&A_Fp12);
    Fp12_init(&B_Fp12);
    Fp12_init(&test_Fp12);
    
    mpz_t scalar;
    mpz_init(scalar);

    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);

    Fp12 Z;
    Fp12_init(&Z);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    printf("====================================================================================\n");
    printf("BLS12_G1_SCM test\n");

for(i=0;i<g1;i++){
    
    mpz_urandomm(scalar,state,order_z);
    BLS12_EFp12_generate_G1(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BLS12_G1_SCM(&test_EFp12,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm1_time+=timedifference_msec(tv_A,tv_B);

}

    printf("====================================================================================\n");
    printf("BLS12_G2_SCM test\n");

for(i=0;i<g2;i++){

    mpz_urandomm(scalar,state,order_z);
    EFp12_generate_G2(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BLS12_G2_SCM(&test_EFp12,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm2_time+=timedifference_msec(tv_A,tv_B);

}
    printf("====================================================================================\n");
    printf("BLS12_G3_exp test\n");

for(i=0;i<g3;i++){
    mpz_urandomm(scalar,state,order_z);
    BLS12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);
    BLS12_PAIRING(&A_Fp12,&P,&Q);

    gettimeofday(&tv_A,NULL);
    BLS12_G3_EXP(&test_Fp12,&A_Fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_time+=timedifference_msec(tv_A,tv_B);
}

    printf("====================================================================================\n");
    printf("BLS12_Opt-ate pairing\n\n");

for(i=0;i<pairing;i++){

    BLS12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);

    gettimeofday(&tv_A,NULL);
    BLS12_PAIRING(&test_Fp12,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_time+=timedifference_msec(tv_A,tv_B);
}

    printf("BLS12 G1 scm.    : %.4f[ms]\n",scm1_time/g1);
    printf("BLS12 G2 scm.      : %.4f[ms]\n",scm2_time/g2);
    printf("BLS12 G3 exp.      : %.4f[ms]\n",exp_time/g3);
    printf("BLS12 optimal ate. : %.4f[ms]\n",opt_time/pairing);

	return 0;
}