#include <ELiPS/bn12_test.h>
int bn12_test_rational_point(){
    printf("====================================================================================\n");
    efp12_t test_G1,test_G2;
    efp12_init(&test_G1);
    efp12_init(&test_G2);
    
    bn12_generate_g1(&test_G1);
    efp12_printf("G1\n",&test_G1);
    printf("\n");
    efp12_scm(&test_G1,&test_G1,order_z);
    efp12_printf("G1 test\n",&test_G1);
    printf("\n");
    
    bn12_generate_g2(&test_G2);
    efp12_printf("G2\n",&test_G2);
    printf("\n");
    efp12_scm(&test_G2,&test_G2,order_z);
    efp12_printf("G2 test\n",&test_G2);
    printf("\n");

    if(test_G1.infinity!=1 || test_G2.infinity!=1) return 1;
    else return 0;
}

void bn12_test_ate_pairing(){
    printf("====================================================================================\n");
    printf("Plain-ate pairing\n\n");
    efp12_t P,Q,s1P,s2P,s1Q,s2Q;
    efp12_init(&P);
    efp12_init(&Q);
    efp12_init(&s1P);
    efp12_init(&s2P);
    efp12_init(&s1Q);
    efp12_init(&s2Q);
    fp12_t Z,test1,test2,test3;
    fp12_init(&Z);
    fp12_init(&test1);
    fp12_init(&test2);
    fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order_z);
	mpz_urandomm(s2,state,order_z);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order_z);
    
    bn12_generate_g1(&P);
    efp12_printf("P\n",&P);
    printf("\n\n");
    bn12_generate_g2(&Q);
    efp12_printf("Q\n",&Q);
    printf("\n\n");
    efp12_scm(&s1P,&P,s1);
    efp12_scm(&s2P,&P,s2);
    efp12_scm(&s1Q,&Q,s1);
    efp12_scm(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain_ate(P,Q)^s1*s2\n");
    bn12_ate_pairing(&Z,&P,&Q);
    fp12_pow(&test1,&Z,s12);
    fp12_printf("",&test1);
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("plain_ate([s1]P,[s2Q])\n");
    bn12_ate_pairing(&test2,&s1P,&s2Q);
    fp12_printf("",&test2);
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("plain_ate([s2]P,[s1]Q)\n");
    bn12_ate_pairing(&test3,&s2P,&s1Q);
    fp12_printf("",&test3);
    printf("\n\n");
    
    printf("bilinear test\n");
    if(fp12_cmp(&test1,&test2)==0 && fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
}

int bn12_test_opt_ate_pairing(int pairing){
    int i,n=0;
    float opt_time=0,opt_compress_time=0,opt_compress_lazy_time=0,opt_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("bn12_Opt-ate pairing\n\n");

    efp12_t P,Q,s1P,s2P,s1Q,s2Q;
    efp12_init(&P);
    efp12_init(&Q);
    efp12_init(&s1P);
    efp12_init(&s2P);
    efp12_init(&s1Q);
    efp12_init(&s2Q);

    fp12_t Z,testA,testB,testC,test1,test2,test3,test4;
    fp12_init(&Z);
    fp12_init(&testA);
    fp12_init(&testB);
    fp12_init(&testC);
    fp12_init(&test1);
    fp12_init(&test2);
    fp12_init(&test3);
    fp12_init(&test4);

    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

    mpz_urandomm(s1,state,order_z);
    mpz_urandomm(s2,state,order_z);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order_z);

    bn12_generate_g1(&P);
    bn12_generate_g2(&Q);

    efp12_scm(&s1P,&P,s1);
    efp12_scm(&s2P,&P,s2);
    efp12_scm(&s1Q,&Q,s1);
    efp12_scm(&s2Q,&Q,s2);

    
    bn12_optate_pairing(&Z,&P,&Q);
    fp12_pow(&testA,&Z,s12);
    bn12_optate_pairing(&testB,&s1P,&s2Q);    
    bn12_optate_pairing(&testC,&s2P,&s1Q);
    
    printf("bilinear test\n");
    if(fp12_cmp(&testA,&testB)!=0 && fp12_cmp(&testA,&testC)!=0){
        printf("bilinear failed!!\n\n");
	return 1;
    }


for(i=0;i<pairing;i++){

    bn12_generate_g1(&P);
    bn12_generate_g2(&Q);

    gettimeofday(&tv_A,NULL);
    bn12_optate_pairing(&test1,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bn12_optate_pairing_lazy(&test2,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    
    if(fp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    fp12_printf("",&test1);
    	fp12_printf("\n",&test2);
    	printf("\n\n");
    	return 1;
    }
}

    printf("bn12 opt ate.                : %.4f[ms]\n",opt_time/pairing);
    printf("bn12 opt ate lazy.           : %.4f[ms]\n",opt_lazy_time/pairing);

    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);

    return 0;
}

void bn12_test_x_ate_pairing(){
    printf("====================================================================================\n");
    printf("X-ate pairing\n\n");
    efp12_t P,Q,s1P,s2P,s1Q,s2Q;
    efp12_init(&P);
    efp12_init(&Q);
    efp12_init(&s1P);
    efp12_init(&s2P);
    efp12_init(&s1Q);
    efp12_init(&s2Q);
    fp12_t Z,test1,test2,test3;
    fp12_init(&Z);
    fp12_init(&test1);
    fp12_init(&test2);
    fp12_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order_z);
	mpz_urandomm(s2,state,order_z);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order_z);
    
    bn12_generate_g1(&P);
    efp12_printf("P\n",&P);
    printf("\n\n");
    bn12_generate_g2(&Q);
    efp12_printf("Q\n",&Q);
    printf("\n\n");
    efp12_scm(&s1P,&P,s1);
    efp12_scm(&s2P,&P,s2);
    efp12_scm(&s1Q,&Q,s1);
    efp12_scm(&s2Q,&Q,s2);
    /*
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("x_ate(P,Q)^s1*s2\n");
    bn12_X_ate_pairing(&Z,&P,&Q);
    fp12_pow(&test1,&Z,s12);
    fp12_printf("",&test1);
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("x_ate([s1]P,[s2]Q)\n");
    bn12_X_ate_pairing(&test2,&s1P,&s2Q);
    fp12_printf("",&test2);
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("x_ate([s2]P,[s1]Q)\n");
    bn12_X_ate_pairing(&test3,&s2P,&s1Q);
    fp12_printf("",&test3);
    printf("\n\n");
    */
    printf("bilinear test\n");
    if(fp12_cmp(&test1,&test2)==0 && fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
}

/*----------------------------------------------------------------------------*/
//bn12_scm
int bn12_test_g1_scm(int scm){
    int i,n=0;
    float scm_time=0,scm_2split_time=0,scm_2split_JSF_time=0,scm_lazy_time=0,scm_jacobian_time=0,scm_jacobian_table_time=0;
    struct timeval tv_A,tv_B;
    printf("================================================================================\n");
    printf("G1 scm test\n\n");
    efp12_t A_efp12,test1,test2,test3,test4,test5,test6;
    efp12_init(&A_efp12);
    efp12_init(&test1);
    efp12_init(&test2);
    efp12_init(&test3);
    efp12_init(&test4);
    efp12_init(&test5);
    efp12_init(&test6);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<scm;i++){
    
    mpz_urandomm(scalar,state,order_z);
    bn12_generate_g1(&A_efp12);

    gettimeofday(&tv_A,NULL);
    bn12_g1_scm_plain(&test1,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g1_scm_2split(&test2,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g1_scm_2split_JSF(&test3,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_JSF_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g1_scm_2split_JSF_lazy(&test4,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(efp12_cmp(&test1,&test2)!=0 || efp12_cmp(&test1,&test3)!=0 || efp12_cmp(&test1,&test4)!=0){
        printf("failed!\n\n");
	efp12_printf("test1=",&test1);
	efp12_printf("\ntest2=",&test2);
	efp12_printf("\ntest3=",&test3);
	efp12_printf("\ntest4=",&test4);
	printf("\n\n");
	return 1;
    }
}
    printf("bn12 G1 scm.                           : %.4f[ms]\n",scm_time/scm);
    printf("bn12 G1 scm 2split.                    : %.4f[ms]\n",scm_2split_time/scm);
    printf("bn12 G1 scm 2split JSF.                : %.4f[ms]\n",scm_2split_JSF_time/scm);
    printf("bn12 G1 scm 2split JSF lazy.           : %.4f[ms]\n",scm_lazy_time/scm);

    mpz_clear(scalar);

    return 0;
}

int bn12_test_g2_scm(int scm){
    int i,n=0;
    float scm_time=0,scm_2split_time=0,scm_2split_JSF_time=0,scm_4split_time=0,scm_lazy_time=0,scm_jacobian_time=0,scm_jacobian_table_time=0;
    struct timeval tv_A,tv_B;
    printf("================================================================================\n");
    printf("G2 scm\n\n");
    efp12_t A_efp12,test1,test2,test3,test4,test5,test6,test7;
    efp12_init(&A_efp12);
    efp12_init(&test1);
    efp12_init(&test2);
    efp12_init(&test3);
    efp12_init(&test4);
    efp12_init(&test5);
    efp12_init(&test6);
    efp12_init(&test7);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<scm;i++){

    mpz_urandomm(scalar,state,order_z);
    //mpz_tdiv_q_2exp(scalar,scalar,30);//relic
    bn12_generate_g2(&A_efp12);

    gettimeofday(&tv_A,NULL);
    bn12_g2_scm_plain(&test1,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g2_scm_2split(&test2,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g2_scm_2split_JSF(&test3,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_JSF_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g2_scm_4split(&test4,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_4split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g2_scm_4split_lazy(&test5,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(efp12_cmp(&test1,&test2)!=0 || efp12_cmp(&test1,&test3)!=0 || efp12_cmp(&test1,&test4)!=0 || efp12_cmp(&test1,&test5)!=0){
        printf("failed!\n\n");
	efp12_printf("test1=",&test1);
	efp12_printf("\ntest2=",&test2);
	efp12_printf("\ntest3=",&test3);
	efp12_printf("\ntest4=",&test4);
	efp12_printf("\ntest5=",&test5);
	printf("\n\n");
	return 1;
    }
}
    printf("bn12 G2 scm.                           : %.4f[ms]\n",scm_time/scm);
    printf("bn12 G2 scm 2split.                    : %.4f[ms]\n",scm_2split_time/scm);
    printf("bn12 G2 scm 2split JSF.                : %.4f[ms]\n",scm_2split_JSF_time/scm);
    printf("bn12 G2 scm 4split.                    : %.4f[ms]\n",scm_4split_time/scm);
    printf("bn12 G2 scm 4split lazy.               : %.4f[ms]\n",scm_lazy_time/scm);
    mpz_clear(scalar);

    return 0;
}



int bn12_test_g3_exp(int exp){
    int i,n=0;
    float exp_time=0,exp_2split_time=0,exp_2split_JSF_time=0,exp_4split_time=0,exp_lazy_time=0,exp_gs_time=0,exp_gs_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("================================================================================\n");
    printf("G3 Exp.\n\n");
    efp12_t P,Q;
    fp12_t A_fp12,test1,test2,test3,test4,test5,test6,test7;
    efp12_init(&P);
    efp12_init(&Q);
    fp12_init(&A_fp12);
    fp12_init(&test1);
    fp12_init(&test2);
    fp12_init(&test3);
    fp12_init(&test4);
    fp12_init(&test5);
    fp12_init(&test6);
    fp12_init(&test7);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    

for(i=0;i<exp;i++){

    mpz_urandomm(scalar,state,order_z);
    bn12_generate_g1(&P);
    bn12_generate_g2(&Q);
    bn12_optate_pairing(&A_fp12,&P,&Q);

    
    gettimeofday(&tv_A,NULL);
    bn12_g3_exp_plain(&test1,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g3_exp_2split(&test2,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_2split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g3_exp_2split_JSF(&test3,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_2split_JSF_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bn12_g3_exp_4split(&test4,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_4split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bn12_g3_exp_4split_lazy(&test5,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    if(fp12_cmp(&test1,&test2)!=0 || fp12_cmp(&test1,&test3)!=0 || fp12_cmp(&test1,&test4)!=0 || fp12_cmp(&test1,&test5)!=0){
        printf("failed!\n\n");
	fp12_printf("test1=",&test1);
	fp12_printf("\ntest2=",&test2);
	fp12_printf("\ntest3=",&test3);
	fp12_printf("\ntest4=",&test4);
	fp12_printf("\ntest5=",&test5);
	printf("\n\n");
	return 1;
    }
}
    printf("bn12 G3 exp.                           : %.4f[ms]\n",exp_time/exp);
    printf("bn12 G3 exp 2split.                    : %.4f[ms]\n",exp_2split_time/exp);
    printf("bn12 G3 exp 2split JSF.                : %.4f[ms]\n",exp_2split_JSF_time/exp);
    printf("bn12 G3 exp 4split.                    : %.4f[ms]\n",exp_4split_time/exp);
    printf("bn12 G3 exp 4split lazy.               : %.4f[ms]\n",exp_lazy_time/exp);

    mpz_clear(scalar);

    return 0;
}
