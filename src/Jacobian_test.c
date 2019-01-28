#include <ELiPS/Jacobian_test.h>

int test_EFp_Jacobian(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_Jacobian_time=0,ecd_Jacobian_lazy_time=0;
    float eca_time=0,eca_Jacobian_time=0,eca_Jacobian_lazy_time=0;
    float scm_time=0,scm_Jacobian_time=0,scm_Jacobian_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp A_EFp,B_EFp,test1,test2,test3;
    EFpZ A_EFpZ,B_EFpZ,testZ1,testZ2,testZ3;
    EFp_init(&A_EFp);
    EFp_init(&B_EFp);
    EFp_init(&test1);
    EFp_init(&test2);
    EFp_init(&test3);
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
    EFp_ECD(&test1,&B_EFp);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

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

    if(EFp_cmp(&test1,&test2)!=0&&EFp_cmp(&test1,&test3)!=0){
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
    EFp_ECA(&test1,&A_EFp,&B_EFp);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);

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

    if(EFp_cmp(&test2,&test3)!=0&&EFp_cmp(&test1,&test3)!=0){
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
    printf("EFp ECA Jacobian.      : %.4f[ms]\n",eca_Jacobian_time/eca);
    printf("EFp ECA Jacobian lazy. : %.4f[ms]\n",eca_Jacobian_lazy_time/eca);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp_SCM test\n");

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    EFp_rational_point(&A_EFp);

    gettimeofday(&tv_A,NULL);
    EFp_SCM(&test1,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_SCM_Jacobian(&test2,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_SCM_Jacobian_lazy(&test3,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp_cmp(&test1,&test3)!=0){
        printf("failed!\n\n");
	EFp_printf(&test1,"");
	EFp_printf(&test3,"\n");
	printf("\n\n");
	return 1;
    }
}

    printf("EFp SCM.               : %.4f[ms]\n",scm_time/scm);
    printf("EFp SCM Jacobian.      : %.4f[ms]\n",scm_Jacobian_time/scm);
    printf("EFp SCM Jacobian Lazy. : %.4f[ms]\n",scm_Jacobian_lazy_time/scm);

    return 0;
}
int test_EFp2_Jacobian(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_Jacobian_time=0,ecd_Jacobian_lazy_time=0;
    float eca_time=0,eca_Jacobian_time=0,eca_Jacobian_lazy_time=0;
    float scm_time=0,scm_Jacobian_time=0,scm_Jacobian_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp2 A_EFp2,B_EFp2,test1,test2,test3;
    EFpZ2 A_EFpZ2,B_EFpZ2,testZ1,testZ2,testZ3;
    EFp2_init(&A_EFp2);
    EFp2_init(&B_EFp2);
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
    EFp2_ECD(&test1,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD_Jacobian(&testZ2,&B_EFpZ2);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    EFp2_Jacobian(&test2,&testZ2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD_Jacobian_lazy(&testZ3,&B_EFpZ2);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);

    EFp2_Jacobian(&test3,&testZ3);

    if(EFp2_cmp(&test1,&test2)!=0&&EFp2_cmp(&test1,&test2)!=0){
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
    printf("EFp2 ECD Jacobian.      : %.4f[ms]\n",ecd_Jacobian_time/ecd);
    printf("EFp2 ECD Jacobian lazy. : %.4f[ms]\n",ecd_Jacobian_lazy_time/ecd);


printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_ECA test\n");
    EFp2_rational_point(&A_EFp2);
    EFp2_to_EFpZ2(&A_EFpZ2,&A_EFp2);

for(i=0;i<ecd;i++){

    EFp2_rational_point(&B_EFp2);
    EFp2_to_EFpZ2(&B_EFpZ2,&B_EFp2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA(&test1,&A_EFp2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Jacobian_scm(&testZ2,&A_EFpZ2,&B_EFpZ2);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    EFp2_Jacobian(&test2,&testZ2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Jacobian_lazy(&testZ3,&A_EFpZ2,&B_EFpZ2);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);

    EFp2_Jacobian(&test3,&testZ3);

    if(EFp2_cmp(&test1,&test2)!=0&&EFp2_cmp(&test1,&test2)!=0){
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

    printf("EFp2 ECA.               : %.4f[ms]\n",eca_time/ecd);
    printf("EFp2 ECA Jacobian.      : %.4f[ms]\n",eca_Jacobian_time/ecd);
    printf("EFp2 ECA Jacobian lazy. : %.4f[ms]\n",eca_Jacobian_lazy_time/ecd);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_SCM test\n");

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    EFp2_rational_point(&A_EFp2);

    gettimeofday(&tv_A,NULL);
    EFp2_SCM(&test1,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_SCM_Jacobian(&test2,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_SCM_Jacobian_lazy(&test3,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp2_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp2_printf(&test1,"");
	EFp2_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}

    printf("EFp2 SCM.               : %.4f[ms]\n",scm_time/scm);
    printf("EFp2 SCM Jacobian.      : %.4f[ms]\n",scm_Jacobian_time/scm);
    printf("EFp2 SCM Jacobian Lazy. : %.4f[ms]\n",scm_Jacobian_lazy_time/scm);

    return 0;
}


int test_BLS12_G1_SCM_Jacobian(int scm){
    int i,n=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BLS12_G1_SCM test\n");
    EFp12 A_EFp12,B_EFp12,test1,test2;
    EFp12_init(&A_EFp12);
    EFp12_init(&B_EFp12);
    EFp12_init(&test1);
    EFp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<scm;i++){
    
    mpz_urandomm(scalar,state,order_z);
    BLS12_EFp12_generate_G1(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BLS12_EFp12_G1_SCM_2split_JSF(&test1,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);


    gettimeofday(&tv_A,NULL);
    BLS12_EFp12_G1_SCM_2split_JSF_Jacobian_lazy(&test2,&A_EFp12,scalar);
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
    printf("BLS12 G1 SCM.               : %.4f[ms]\n",scm_time/scm);
    printf("BLS12 G1 SCM Jacobian lazy. : %.4f[ms]\n",scm_lazy_time/scm);

    return 0;
}


int test_BLS12_G2_SCM_Jacobian(int scm){
    int i,n=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BLS12_G2_SCM test\n");
    EFp12 A_EFp12,B_EFp12,test1,test2;
    EFp12_init(&A_EFp12);
    EFp12_init(&B_EFp12);
    EFp12_init(&test1);
    EFp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<scm;i++){

    mpz_urandomm(scalar,state,order_z);
    EFp12_generate_G2(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BLS12_EFp12_G2_SCM_4split_lazy(&test1,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);


    gettimeofday(&tv_A,NULL);
    BLS12_EFp12_G2_SCM_4split_Jacobian_lazy(&test2,&A_EFp12,scalar);
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
    printf("BLS12 G2 SCM.               : %.4f[ms]\n",scm_time/scm);
    printf("BLS12 G2 SCM Jacobian lazy. : %.4f[ms]\n",scm_lazy_time/scm);

    return 0;
}

