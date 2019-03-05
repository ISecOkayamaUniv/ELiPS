#include <ELiPS/lazy_test.h>


int test_Field_Lazy(int fp2,int fp6,int fp12){
    int i,n=0;
    float mul_time=0,mul_lazy_time=0;
    float sqr_time=0,sqr_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("Field Lazy test\n");
    Fp2 A_Fp2,B_Fp2,test1,test2;
    Fp6 A_Fp6,B_Fp6,test1_Fp6,test2_Fp6,test3_Fp6;
    Fp12 A_Fp12,B_Fp12,test1_Fp12,test2_Fp12,test3_Fp12,test4_Fp12;

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
    Fp12_init(&test1_Fp12);
    Fp12_init(&test2_Fp12);
    Fp12_init(&test3_Fp12);


    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_mul test\n");

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

for(i=0;i<n;i++){
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

return 0;
}


int test_EFp_Lazy(int eca,int ecd,int scm){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0;
    float eca_time=0,eca_lazy_time=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp A_EFp,B_EFp,test1,test2;
    EFp_init(&A_EFp);
    EFp_init(&B_EFp);
    EFp_init(&test1);
    EFp_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

printf("------------------------------------------------------------------------------------\n");
    printf("EFp_ECD test\n");

for(i=0;i<ecd;i++){

    EFp_rational_point(&B_EFp);

    gettimeofday(&tv_A,NULL);
    EFp_ECD(&test1,&B_EFp);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_ECD_lazy(&test2,&B_EFp);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(EFp_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp_printf(&test1,"");
	EFp_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp ECD.      : %.4f[ms]\n",ecd_time/ecd);
    printf("EFp ECD lazy. : %.4f[ms]\n",ecd_lazy_time/ecd);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp_ECA test\n");

    EFp_rational_point(&A_EFp);
for(i=0;i<eca;i++){

    EFp_rational_point(&B_EFp);

    gettimeofday(&tv_A,NULL);
    EFp_ECA(&test1,&A_EFp,&B_EFp);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp_ECA_lazy(&test2,&A_EFp,&B_EFp);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp_printf(&test1,"");
	EFp_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp ECA.      : %.4f[ms]\n",eca_time/eca);
    printf("EFp ECA lazy. : %.4f[ms]\n",eca_lazy_time/eca);

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
    EFp_SCM_lazy(&test2,&A_EFp,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp_printf(&test1,"");
	EFp_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp SCM.      : %.4f[ms]\n",scm_time/scm);
    printf("EFp SCM lazy. : %.4f[ms]\n",scm_lazy_time/scm);

	return 0;
}

int test_EFp2_Lazy(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0;
    float eca_time=0,eca_lazy_time=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp2 A_EFp2,B_EFp2,test1,test2;
    EFp2_init(&A_EFp2);
    EFp2_init(&B_EFp2);
    EFp2_init(&test1);
    EFp2_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_ECD test\n");

for(i=0;i<ecd;i++){

    EFp2_rational_point(&B_EFp2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD(&test1,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD_lazy(&test2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(EFp2_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp2_printf(&test1,"");
	EFp2_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp2 ECD.      : %.4f[ms]\n",ecd_time/ecd);
    printf("EFp2 ECD lazy. : %.4f[ms]\n",ecd_lazy_time/ecd);

printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_ECA test\n");

    EFp2_rational_point(&A_EFp2);
for(i=0;i<eca;i++){

    EFp2_rational_point(&B_EFp2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA(&test1,&A_EFp2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA_lazy(&test2,&A_EFp2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp2_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp2_printf(&test1,"");
	EFp2_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp2 ECA.      : %.4f[ms]\n",eca_time/eca);
    printf("EFp2 ECA lazy. : %.4f[ms]\n",eca_lazy_time/eca);

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
    EFp2_SCM_lazy(&test2,&A_EFp2,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp2_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp2_printf(&test1,"");
	EFp2_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("EFp2 SCM.      : %.4f[ms]\n",scm_time/scm);
    printf("EFp2 SCM lazy. : %.4f[ms]\n",scm_lazy_time/scm);

	return 0;
}

int test_EFp12_Lazy(int ecd,int eca,int scm){
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

void test_BN12_G1_SCM_Lazy(){
    int i,n=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BN12_G1_SCM test\n");
    EFp12 A_EFp12,B_EFp12,test1,test2;
    EFp12_init(&A_EFp12);
    EFp12_init(&B_EFp12);
    EFp12_init(&test1);
    EFp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


n=1000;
for(i=0;i<n;i++){
    
    mpz_urandomm(scalar,state,order_z);
    BN12_EFp12_generate_G1(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BN12_EFp12_G1_SCM_2split_JSF(&test1,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);


    gettimeofday(&tv_A,NULL);
    BN12_EFp12_G1_SCM_2split_JSF_lazy(&test2,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp12_printf(&test1,"");
	EFp12_printf(&test2,"\n");
	printf("\n\n");
	break;
    }
}
    printf("BN12 G1 SCM.      : %.4f[ms]\n",scm_time/n);
    printf("BN12 G1 SCM lazy. : %.4f[ms]\n",scm_lazy_time/n);

}

int test_BLS12_G1_SCM_Lazy(int scm){
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
    BLS12_EFp12_G1_SCM_2split_JSF_lazy(&test2,&A_EFp12,scalar);
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
    printf("BLS12 G1 SCM.      : %.4f[ms]\n",scm_time/scm);
    printf("BLS12 G1 SCM lazy. : %.4f[ms]\n",scm_lazy_time/scm);

	return 0;
}

void test_BN12_G2_SCM_Lazy(){
    int i,n=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BN12_G2_SCM test\n");
    EFp12 A_EFp12,B_EFp12,test1,test2;
    EFp12_init(&A_EFp12);
    EFp12_init(&B_EFp12);
    EFp12_init(&test1);
    EFp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

n=10;

for(i=0;i<n;i++){

    mpz_urandomm(scalar,state,order_z);
    EFp12_generate_G2(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BN12_EFp12_G2_SCM_4split(&test1,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);


    gettimeofday(&tv_A,NULL);
    BN12_EFp12_G2_SCM_4split_lazy(&test2,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(EFp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	EFp12_printf(&test1,"");
	EFp12_printf(&test2,"\n");
	printf("\n\n");
	break;
    }
}
    printf("BN12 G2 SCM.      : %.4f[ms]\n",scm_time/n);
    printf("BN12 G2 SCM lazy. : %.4f[ms]\n",scm_lazy_time/n);

}

int test_BLS12_G2_SCM_Lazy(int scm){
    int i,n=0;
    float scm_time=0,scm_lazy_time=0,scm_lazy_time2=0,scm_lazy_time3=0;
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
    BLS12_EFp12_G2_SCM_4split(&test1,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);


    gettimeofday(&tv_A,NULL);
    BLS12_EFp12_G2_SCM_4split_lazy(&test2,&A_EFp12,scalar);
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
    printf("BLS12 G2 SCM.       : %.4f[ms]\n",scm_time/scm);
    printf("BLS12 G2 SCM lazy.  : %.4f[ms]\n",scm_lazy_time/scm);

	return 0;
}

void test_BN12_G3_exp_Lazy(){
    int i,n=0;
    float exp_time=0,exp_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BN12_G3_exp test\n");
    Fp12 A_Fp12,B_Fp12,test1,test2;
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12_init(&A_Fp12);
    Fp12_init(&B_Fp12);
    Fp12_init(&test1);
    Fp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

n=10;

for(i=0;i<n;i++){

    mpz_urandomm(scalar,state,order_z);
    BN12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);
    BN12_Opt_ate_pairing(&A_Fp12,&Q,&P);

    gettimeofday(&tv_A,NULL);
    BN12_Fp12_G3_EXP_4split(&test1,&A_Fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_time+=timedifference_msec(tv_A,tv_B);


    gettimeofday(&tv_A,NULL);
    BN12_Fp12_G3_EXP_4split_lazy(&test2,&A_Fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	Fp12_printf(&test1,"");
	Fp12_printf(&test2,"\n");
	printf("\n\n");
	break;
    }
}
    printf("BN12 G3 exp.      : %.4f[ms]\n",exp_time/n);
    printf("BN12 G3 exp lazy. : %.4f[ms]\n",exp_lazy_time/n);

}

int test_BLS12_G3_exp_Lazy(int exp){
    int i,n=0;
    float exp_time=0,exp_lazy_time=0,exp_lazy_time2=0,exp_lazy_time3=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BLS12_G3_exp test\n");
    Fp12 A_Fp12,B_Fp12,test1,test2,test3,test4;
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12_init(&A_Fp12);
    Fp12_init(&B_Fp12);
    Fp12_init(&test1);
    Fp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


for(i=0;i<exp;i++){

    mpz_urandomm(scalar,state,order_z);
    BLS12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);
    BLS12_Opt_ate_pairing(&A_Fp12,&P,&Q);

    gettimeofday(&tv_A,NULL);
    BLS12_Fp12_G3_EXP_4split(&test1,&A_Fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    BLS12_Fp12_G3_EXP_4split_lazy(&test2,&A_Fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	Fp12_printf(&test1,"");
	Fp12_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}
    printf("BLS12 G3 exp.       : %.4f[ms]\n",exp_time/exp);
    printf("BLS12 G3 exp lazy.  : %.4f[ms]\n",exp_lazy_time/exp);


	return 0;

}


void test_BN12_opt_ate_pairing_Lazy(){
    int i,n=0;
    float opt_time=0,opt_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BN12_Opt-ate pairing\n\n");

    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);

    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);


    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
n=10;
for(i=0;i<n;i++){

    BN12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);

    gettimeofday(&tv_A,NULL);
    BN12_Opt_ate_pairing(&test1,&Q,&P);
    gettimeofday(&tv_B,NULL);
    opt_time+=timedifference_msec(tv_A,tv_B);


    gettimeofday(&tv_A,NULL);
    BN12_Opt_ate_pairing_lazy(&test2,&Q,&P);
    gettimeofday(&tv_B,NULL);
    opt_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(Fp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	Fp12_printf(&test1,"");
	Fp12_printf(&test2,"\n");
	printf("\n\n");
	break;
    }
}

    printf("BN12 opt ate.      : %.4f[ms]\n",opt_time/n);
    printf("BN12 opt ate lazy. : %.4f[ms]\n",opt_lazy_time/n);
}


int test_BLS12_opt_ate_pairing_Lazy(int pairing){
    int i,n=0;
    float opt_time=0,opt_lazy_time=0,opt_lazy_time2=0,opt_lazy_time3=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("BLS12_Opt-ate pairing\n\n");

    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);

    Fp12 Z,test1,test2,test3,test4;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    Fp12_init(&test4);


    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<pairing;i++){

    BLS12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);

    gettimeofday(&tv_A,NULL);
    BLS12_Opt_ate_pairing(&test1,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    BLS12_Opt_ate_pairing_lazy(&test2,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(Fp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	Fp12_printf(&test1,"");
	Fp12_printf(&test2,"\n");
	printf("\n\n");
	return 1;
    }
}

    printf("BLS12 opt ate.       : %.4f[ms]\n",opt_time/pairing);
    printf("BLS12 opt ate lazy.  : %.4f[ms]\n",opt_lazy_time/pairing);

	return 0;
}
