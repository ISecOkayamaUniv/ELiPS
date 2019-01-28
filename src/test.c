#include <ELiPS/test.h>
/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) * 1000.0f + (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_usec - tv_start.tv_usec);
}

/*----------------------------------------------------------------------------*/
//test
void test_Field(){
    printf("====================================================================================\n");
    Fp6 A_Fp6,test1,test2;
    Fp6_init(&A_Fp6);
    Fp6_init(&test1);
    Fp6_init(&test2);
    mpz_t exp;
    mpz_init(exp);
    
    gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp6_set_random(&A_Fp6,state);
    Fp6_printf(&A_Fp6,"");
    printf("\n\n");
    
    printf("mul/sqr\n");
    Fp6_mul(&test1,&A_Fp6,&A_Fp6);
    Fp6_printf(&test1,"");
    printf("\n");
    
    Fp6_sqr(&test2,&A_Fp6);
    Fp6_printf(&test2,"");
    printf("\n");
    
    if(Fp6_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("pow/inv\n");
    mpz_pow_ui(exp,prime_z,12);
    mpz_sub_ui(exp,exp,2);
    Fp6_pow(&test1,&A_Fp6,exp);
    Fp6_printf(&test1,"");
    printf("\n");
    
    Fp6_inv(&test2,&A_Fp6);
    Fp6_printf(&test2,"");
    printf("\n");
    
    if(Fp6_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
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

void test_All(){
	int test_point,test_G1,test_G2,test_G3,test_pairing;
	int Lazy_Field,Lazy_EFp,Lazy_EFp2,Lazy_EFp12,Lazy_G1,Lazy_G2,Lazy_G3,Lazy_pairing;
	int Jacobian_EFp,Jacobian_EFp2,Jacobian_G1;
	
printf("====================================================================================\n");
    printf("Normal test\n");
	
	test_point = BLS12_test_rational_point();
	test_G1 = BLS12_test_G1_SCM();
	test_G2 = BLS12_test_G2_SCM();
	test_G3 = BLS12_test_G3_EXP();
	test_pairing = BLS12_test_opt_ate_pairing();

printf("====================================================================================\n");
    printf("Lazy test\n");
	
	Lazy_Field = test_Field_Lazy(10000,10000,10000);
	Lazy_EFp = test_EFp_Lazy(1000,1000,100);
	Lazy_EFp2 = test_EFp2_Lazy(100,100,100);
	Lazy_EFp12 = test_EFp12_Lazy(10,10,10);
	Lazy_G1 = test_BLS12_G1_SCM_Lazy(100);
	Lazy_G2 = test_BLS12_G2_SCM_Lazy(10);
	Lazy_G3 = test_BLS12_G3_exp_Lazy(10);
	Lazy_pairing = test_BLS12_opt_ate_pairing_Lazy(10);

printf("====================================================================================\n");
    printf("Jacobian test\n");

	
	Jacobian_EFp = test_EFp_Jacobian(1000,1000,100);
	Jacobian_EFp2 = test_EFp2_Jacobian(1000,1000,100);
	Jacobian_G1 = test_BLS12_G1_SCM_Jacobian(100);


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
