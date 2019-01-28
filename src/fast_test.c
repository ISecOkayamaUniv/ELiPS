#include <ELiPS/fast_test.h>

int Fast_test_BLS12_G1_SCM(int scm){
    int i;
    float scm_time=0;
    struct timeval tv_A,tv_B;
    EFp12 A_EFp12,test;
    EFp12_init(&A_EFp12);
    EFp12_init(&test);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

    printf("====================================================================================\n");
    printf("BLS12_G1_SCM test\n");

for(i=0;i<scm;i++){
    
    mpz_urandomm(scalar,state,order_z);
    BLS12_EFp12_generate_G1(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BLS12_G1_SCM(&test,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

}
    printf("BLS12 FAST G1 SCM.    : %.4f[ms]\n",scm_time/scm);

	return 0;
}

int Fast_test_BLS12_G2_SCM(int scm){
    int i;
    float scm_time=0;
    struct timeval tv_A,tv_B;
    EFp12 A_EFp12,test;
    EFp12_init(&A_EFp12);
    EFp12_init(&test);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

    printf("====================================================================================\n");
    printf("BLS12_G2_SCM test\n");

for(i=0;i<scm;i++){

    mpz_urandomm(scalar,state,order_z);
    EFp12_generate_G2(&A_EFp12);

    gettimeofday(&tv_A,NULL);
    BLS12_G2_SCM(&test,&A_EFp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

}
    printf("BLS12 G2 SCM.      : %.4f[ms]\n",scm_time/scm);

	return 0;
}

int Fast_test_BLS12_G3_EXP(int exp){
    int i;
    float exp_time=0,exp_lazy_time=0;
    struct timeval tv_A,tv_B;
    Fp12 A_Fp12,B_Fp12,test;
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12_init(&A_Fp12);
    Fp12_init(&test);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


    printf("====================================================================================\n");
    printf("BLS12_G3_exp test\n");

for(i=0;i<exp;i++){
    mpz_urandomm(scalar,state,order_z);
    BLS12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);
    BLS12_PAIRING(&A_Fp12,&P,&Q);

    gettimeofday(&tv_A,NULL);
    BLS12_G3_EXP(&test,&A_Fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_time+=timedifference_msec(tv_A,tv_B);

}
    printf("BLS12 G3 exp.      : %.4f[ms]\n",exp_time/exp);

	return 0;

}

int Fast_test_BLS12_pairing(int pairing){
    int i;
    float opt_time=0;
    struct timeval tv_A,tv_B;

    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);

    Fp12 Z,test;
    Fp12_init(&Z);
    Fp12_init(&test);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

    printf("====================================================================================\n");
    printf("BLS12_Opt-ate pairing\n\n");

for(i=0;i<pairing;i++){

    BLS12_EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);

    gettimeofday(&tv_A,NULL);
    BLS12_PAIRING(&test,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_time+=timedifference_msec(tv_A,tv_B);
}

    printf("BLS12 optimal ate. : %.4f[ms]\n",opt_time/pairing);

	return 0;
}
