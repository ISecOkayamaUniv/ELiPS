#include <ELiPS/bls12_test.h>
int bls12_test_rational_point(){
    printf("====================================================================================\n");
    efp12_t test_g1,test_g2;
    efp12_init(&test_g1);
    efp12_init(&test_g2);
    
    bls12_generate_g1(&test_g1);
    efp12_printf("g1\n",&test_g1);
    printf("\n");
    efp12_scm(&test_g1,&test_g1,order_z);
    efp12_printf("g1 test\n",&test_g1);
    printf("\n");
    
    bls12_generate_g2(&test_g2);
    efp12_printf("g2\n",&test_g2);
    printf("\n");
    efp12_scm(&test_g2,&test_g2,order_z);
    efp12_printf("g2 test\n",&test_g2);
    printf("\n");

    if(test_g1.infinity!=1 || test_g2.infinity!=1) return 1;
    else return 0;
}

void bls12_test_plain_ate_pairing(){
    printf("====================================================================================\n");
    printf("ate pairing\n\n");
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
    
    bls12_generate_g1(&P);
    efp12_printf("P\n",&P);
    printf("\n\n");
    bls12_generate_g2(&Q);
    efp12_printf("Q\n",&Q);
    printf("\n\n");
    efp12_scm(&s1P,&P,s1);
    efp12_scm(&s2P,&P,s2);
    efp12_scm(&s1Q,&Q,s1);
    efp12_scm(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain_ate(P,Q)^s1*s2\n");
    bls12_ate_pairing(&Z,&P,&Q);
    fp12_pow(&test1,&Z,s12);
    fp12_printf("",&test1);
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("plain_ate([s1]P,[s2Q])\n");
    bls12_ate_pairing(&test2,&s1P,&s2Q);
    fp12_printf("",&test2);
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("plain_ate([s2]P,[s1]Q)\n");
    bls12_ate_pairing(&test3,&s2P,&s1Q);
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

int bls12_test_opt_ate_pairing(int pairing){
    int i,n=0;
    float opt_time=0,opt_compress_time=0,opt_compress_lazy_time=0,opt_lazy_time=0,opt_compress_lazy_montgomery_time=0;
    cost tmp,opt_cost,opt_compress_cost,opt_compress_lazy_cost,opt_lazy_cost,opt_compress_lazy_montgomery_cost;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("bls12_Opt-ate pairing\n\n");

    efp12_t P,Q,s1P,s2P,s1Q,s2Q;
    efp12_init(&P);
    efp12_init(&Q);
    efp12_init(&s1P);
    efp12_init(&s2P);
    efp12_init(&s1Q);
    efp12_init(&s2Q);

    fp12_t Z,testA,testB,testC,test1,test2,test3,test4,test5;
    fp12_init(&Z);
    fp12_init(&testA);
    fp12_init(&testB);
    fp12_init(&testC);
    fp12_init(&test1);
    fp12_init(&test2);
    fp12_init(&test3);
    fp12_init(&test4);
    fp12_init(&test5);
    
    cost_init(&tmp);
    cost_init(&opt_cost);
    cost_init(&opt_compress_cost);
    cost_init(&opt_lazy_cost);
    cost_init(&opt_compress_lazy_cost);
    cost_init(&opt_compress_lazy_montgomery_cost);
    
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

    bls12_generate_g1(&P);
    bls12_generate_g2(&Q);

    efp12_scm(&s1P,&P,s1);
    efp12_scm(&s2P,&P,s2);
    efp12_scm(&s1Q,&Q,s1);
    efp12_scm(&s2Q,&Q,s2);

    
    bls12_optate_pairing(&Z,&P,&Q);
    fp12_pow(&testA,&Z,s12);
    bls12_optate_pairing(&testB,&s1P,&s2Q);    
    bls12_optate_pairing(&testC,&s2P,&s1Q);
    
    printf("bilinear test\n");
    if(fp12_cmp(&testA,&testB)!=0 || fp12_cmp(&testA,&testC)!=0){
        printf("bilinear failed!!\n\n");
	return 1;
    }


MILLER_OPT=0;
FINALEXP_OPT=0;
MILLER_OPT_MONTGOMERY=0;
FINALEXP_OPT_MONTGOMERY=0;
cost_init(&MILLER_OPT_MONTGOMERY_COST);
cost_init(&FINALEXP_OPT_MONTGOMERY_COST);
for(i=0;i<pairing;i++){

    bls12_generate_g1(&P);
    bls12_generate_g2(&Q);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_optate_pairing_basic(&test1,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&opt_cost,&tmp);
    
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_optate_pairing_lazy(&test2,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_lazy_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&opt_lazy_cost,&tmp);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_optate_pairing_compress(&test3,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_compress_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&opt_compress_cost,&tmp);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_optate_pairing_compress_lazy(&test4,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_compress_lazy_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&opt_compress_lazy_cost,&tmp);
    
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_optate_pairing_compress_lazy_montgomery(&test5,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_compress_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&opt_compress_lazy_montgomery_cost,&tmp);
    
    if(fp12_cmp(&test1,&test2)!=0 || fp12_cmp(&test1,&test3)!=0 || fp12_cmp(&test1,&test4)!=0){
        printf("failed!\n\n");
	    fp12_printf("",&test1);
    	fp12_printf("\n",&test2);
    	printf("\n\n");
    	return 1;
    }
}
    cost_substruction(&FINALEXP_OPT_MONTGOMERY_COST,&opt_compress_lazy_montgomery_cost,&MILLER_OPT_MONTGOMERY_COST);

    printf("bls12 opt ate.                                    : %.4f[ms]\n",opt_time/pairing);
    printf("bls12 opt ate lazy.                               : %.4f[ms]\n",opt_lazy_time/pairing);
    printf("bls12 opt ate compress.                           : %.4f[ms]\n",opt_compress_time/pairing);
    printf("bls12 opt ate compress lazy.                      : %.4f[ms]\n",opt_compress_lazy_time/pairing);
    printf("bls12 opt ate compress lazy(MILLER_OPTATE).       : %.4f[ms]\n",MILLER_OPT/pairing);
    printf("bls12 opt ate compress lazy(FINALEXP_OPT).        : %.4f[ms]\n",FINALEXP_OPT/pairing);
    printf("bls12 opt ate compress lazy montgomery.           : %.4f[ms]\n",opt_compress_lazy_montgomery_time/pairing);
    printf("bls12 opt ate compress lazy(MILLER_OPTATE_MONT).  : %.4f[ms]\n",MILLER_OPT_MONTGOMERY/pairing);
    printf("bls12 opt ate compress lazy(FINALEXP_OPT_MONT).   : %.4f[ms]\n",FINALEXP_OPT_MONTGOMERY/pairing);

    #ifdef DEBUG_COST_A
    printf("*********bls12 opt ate fp COST.********         \n");
    cost_printf("bls12 opt ate",&opt_cost,pairing);
    cost_printf("bls12 opt ate compress",&opt_compress_cost,pairing);
    cost_printf("bls12 opt ate lazy",&opt_lazy_cost,pairing);
    cost_printf("bls12 opt ate compress lazy",&opt_compress_lazy_cost,pairing);
    cost_printf("bls12 opt ate compress lazy montgomery",&opt_compress_lazy_montgomery_cost,pairing);
    cost_printf("bls12 opt ate compress lazy(MILLER_OPTATE_MONT)",&MILLER_OPT_MONTGOMERY_COST,pairing);
    cost_printf("bls12 opt ate compress lazy(FINALEXP_OPT_MONT)",&FINALEXP_OPT_MONTGOMERY_COST,pairing);
    printf("***************************************         \n");
    #endif
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);

    return 0;
}

/*----------------------------------------------------------------------------*/
//bls12_scm
int bls12_test_g1_scm(int scm){
    int i,n=0;
    float scm_time=0,scm_2split_time=0,scm_2split_jsf_time=0,scm_lazy_time=0,scm_jacobian_time=0,scm_mixture_time=0,scm_jacobian_table_time=0,scm_2split_2naf_time=0,scm_2split_3naf_M_time=0,scm_2split_3naf_I_time=0,scm_2split_5naf_I_time=0,scm_2split_5naf_I_mixture_time=0,scm_2split_5naf_I_mixture_lazy_time=0,scm_2split_5naf_I_mixture_lazy_montgomery_time=0,scm_2split_7naf_I_mixture_time=0;
    cost tmp,scm_cost,scm_2split_cost,scm_2split_jsf_cost,scm_lazy_cost,scm_jacobian_cost,scm_mixture_cost,scm_jacobian_table_cost,scm_2split_2naf_cost,scm_2split_3naf_M_cost,scm_2split_3naf_I_cost,scm_2split_5naf_I_cost,scm_2split_5naf_I_mixture_cost,scm_2split_5naf_I_mixture_lazy_cost,scm_2split_5naf_I_mixture_lazy_montgomery_cost,scm_2split_7naf_I_mixture_cost;
    struct timeval tv_A,tv_B;
    printf("================================================================================\n");
    printf("g1 scm test\n\n");
    efp12_t A_efp12,test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11,test12,test13,test14,test15;
    efp12_init(&A_efp12);
    efp12_init(&test1);
    efp12_init(&test2);
    efp12_init(&test3);
    efp12_init(&test4);
    efp12_init(&test5);
    efp12_init(&test6);
    efp12_init(&test7);
    efp12_init(&test8);
    efp12_init(&test9);
    efp12_init(&test10);
    efp12_init(&test11);
    efp12_init(&test12);
    efp12_init(&test13);
    efp12_init(&test14);
    efp12_init(&test15);

    mpz_t scalar;
    mpz_init(scalar);
    
    cost_init(&tmp);
    cost_init(&scm_cost);
    cost_init(&scm_2split_5naf_I_mixture_lazy_montgomery_cost);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    //mpz_set_ui(scalar,1234567);
    bls12_generate_g1(&A_efp12);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_basic(&test1,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&scm_cost,&tmp);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split(&test2,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_jsf(&test3,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_jsf_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_jsf_lazy(&test4,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_jsf_jacobian_lazy(&test5,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_jacobian_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_jsf_mixture_lazy(&test7,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_mixture_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_2naf(&test8,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_2naf_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_3naf_shamia(&test9,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_3naf_M_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_3naf_interleaving(&test10,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_3naf_I_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_5naf_interleaving(&test11,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_5naf_I_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_5naf_interleaving_mixture(&test12,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_5naf_I_mixture_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_5naf_interleaving_mixture_lazy(&test13,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_5naf_I_mixture_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_5naf_interleaving_mixture_lazy_montgomery(&test14,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_5naf_I_mixture_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&scm_2split_5naf_I_mixture_lazy_montgomery_cost,&tmp);
    
    gettimeofday(&tv_A,NULL);
    bls12_g1_scm_2split_7naf_interleaving_mixture_lazy(&test15,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_7naf_I_mixture_time+=timedifference_msec(tv_A,tv_B);
    
    if(efp12_cmp(&test1,&test2)!=0 || efp12_cmp(&test1,&test3)!=0 || efp12_cmp(&test1,&test4)!=0 || efp12_cmp(&test1,&test5)!=0  || efp12_cmp(&test1,&test7)!=0 || efp12_cmp(&test1,&test8)!=0 || efp12_cmp(&test1,&test9)!=0 || efp12_cmp(&test1,&test10)!=0 || efp12_cmp(&test1,&test11)!=0 || efp12_cmp(&test1,&test12)!=0 || efp12_cmp(&test1,&test13)!=0 || efp12_cmp(&test1,&test14)!=0 || efp12_cmp(&test1,&test15)!=0){
        printf("failed!\n\n");
	efp12_printf("test1=",&test1);
	efp12_printf("\ntest2=",&test2);
	efp12_printf("\ntest3=",&test3);
	efp12_printf("\ntest4=",&test4);
	efp12_printf("\ntest5=",&test5);
	efp12_printf("\ntest7=",&test7);
	efp12_printf("\ntest8=",&test8);
	efp12_printf("\ntest9=",&test9);
	efp12_printf("\ntest10=",&test10);
	efp12_printf("\ntest11=",&test11);
	efp12_printf("\ntest12=",&test12);
	efp12_printf("\ntest13=",&test13);
	efp12_printf("\ntest14=",&test14);
	efp12_printf("\ntest15=",&test15);
	printf("\n\n");
	return 1;
    }
}
    printf("bls12 g1 scm.                                                : %.4f[ms]\n",scm_time/scm);
    printf("bls12 g1 scm 2split.                                         : %.4f[ms]\n",scm_2split_time/scm);
    printf("bls12 g1 scm 2split jsf.                                     : %.4f[ms]\n",scm_2split_jsf_time/scm);
    printf("bls12 g1 scm 2split jsf lazy.                                : %.4f[ms]\n",scm_lazy_time/scm);
    printf("bls12 g1 scm 2split jsf jacobian lazy.                       : %.4f[ms]\n",scm_jacobian_time/scm);
    printf("bls12 g1 scm 2split jsf mixture lazy.                        : %.4f[ms]\n",scm_mixture_time/scm);
    printf("bls12 g1 scm 2split 2naf.                                    : %.4f[ms]\n",scm_2split_2naf_time/scm);
    printf("bls12 g1 scm 2split 3naf shamir.                              : %.4f[ms]\n",scm_2split_3naf_M_time/scm);
    printf("bls12 g1 scm 2split 3naf interleave.                         : %.4f[ms]\n",scm_2split_3naf_I_time/scm);
    printf("bls12 g1 scm 2split 5naf interleave.                         : %.4f[ms]\n",scm_2split_5naf_I_time/scm);
    printf("bls12 g1 scm 2split 5naf interleave mixture.                 : %.4f[ms]\n",scm_2split_5naf_I_mixture_time/scm);
    printf("bls12 g1 scm 2split 5naf interleave mixture lazy.            : %.4f[ms]\n",scm_2split_5naf_I_mixture_lazy_time/scm);
    printf("bls12 g1 scm 2split 5naf interleave mixture lazy montgomery. : %.4f[ms]\n",scm_2split_5naf_I_mixture_lazy_montgomery_time/scm);
    printf("bls12 g1 scm 2split 7naf interleave mixture.                 : %.4f[ms]\n",scm_2split_7naf_I_mixture_time/scm);
    
    #ifdef DEBUG_COST_A
    printf("*********bls12 g2 scm fp COST.********         \n");
    cost_printf("bls12 g1 scm",&scm_cost,scm);
    cost_printf("bls12 g1 scm 2split 5naf interleave mixture lazy montgomery",&scm_2split_5naf_I_mixture_lazy_montgomery_cost,scm);
    printf("***************************************         \n");
    #endif
    
    mpz_clear(scalar);

    return 0;
}

int bls12_test_g2_scm(int scm){
    int i,n=0;
    float scm_time=0,scm_2split_time=0,scm_2split_jsf_time=0,scm_4split_time=0,scm_lazy_time=0,scm_jacobian_time=0,scm_jacobian_table_time=0,scm_mixture_time=0,scm_2split_3naf_I_mixture_time=0,scm_2split_5naf_I_mixture_time=0,scm_2split_5naf_I_mixture_montgomery_time=0;
    cost tmp,scm_cost,scm_2split_cost,scm_2split_jsf_cost,scm_4split_cost,scm_lazy_cost,scm_jacobian_cost,scm_jacobian_table_cost,scm_mixture_cost,scm_2split_3naf_I_mixture_cost,scm_2split_5naf_I_mixture_cost,scm_2split_5naf_I_mixture_montgomery_cost;
    struct timeval tv_A,tv_B;
    printf("================================================================================\n");
    printf("g2 scm\n\n");
    efp12_t A_efp12,test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11;
    efp12_init(&A_efp12);
    efp12_init(&test1);
    efp12_init(&test2);
    efp12_init(&test3);
    efp12_init(&test4);
    efp12_init(&test5);
    efp12_init(&test6);
    efp12_init(&test7);
    efp12_init(&test8);
    efp12_init(&test9);
    efp12_init(&test10);
    efp12_init(&test11);

    mpz_t scalar;
    mpz_init(scalar);
    
    cost_init(&tmp);
    cost_init(&scm_cost);
    cost_init(&scm_2split_5naf_I_mixture_montgomery_cost);
    
    gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    gmp_randseed_ui(state,1);

for(i=0;i<scm;i++){

    mpz_urandomm(scalar,state,order_z);
    //mpz_set_ui(scalar,1234567);
    //mpz_tdiv_q_2exp(scalar,scalar,30);//relic
    bls12_generate_g2(&A_efp12);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_basic(&test1,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&scm_cost,&tmp);

    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_2split(&test2,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_2split_jsf(&test3,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_jsf_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split(&test4,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_4split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split_lazy(&test5,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split_jacobian_lazy(&test6,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_jacobian_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split_mixture_lazy(&test8,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_mixture_time+=timedifference_msec(tv_A,tv_B);
        
    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split_3naf_interleaving_mixture_lazy(&test9,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_3naf_I_mixture_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split_5naf_interleaving_mixture_lazy(&test10,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_5naf_I_mixture_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split_5naf_interleaving_mixture_lazy_montgomery(&test11,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_5naf_I_mixture_montgomery_time+=timedifference_msec(tv_A,tv_B);
    
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g2_scm_4split_5naf_interleaving_mixture_lazy_montgomery(&test11,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_2split_5naf_I_mixture_montgomery_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&scm_2split_5naf_I_mixture_montgomery_cost,&tmp);

    if(efp12_cmp(&test1,&test2)!=0 || efp12_cmp(&test1,&test3)!=0 || efp12_cmp(&test1,&test4)!=0 || efp12_cmp(&test1,&test5)!=0 || efp12_cmp(&test1,&test6)!=0 || efp12_cmp(&test1,&test8)!=0 || efp12_cmp(&test1,&test9)!=0 || efp12_cmp(&test1,&test10)!=0 || efp12_cmp(&test1,&test11)!=0){
        printf("failed!\n\n");
	efp12_printf("test1=",&test1);
	efp12_printf("\ntest2=",&test2);
	efp12_printf("\ntest3=",&test3);
	efp12_printf("\ntest4=",&test4);
	efp12_printf("\ntest5=",&test5);
	efp12_printf("\ntest6=",&test6);
	efp12_printf("\ntest8=",&test8);
	efp12_printf("\ntest9=",&test9);
	efp12_printf("\ntest10=",&test10);
	efp12_printf("\ntest11=",&test11);
	printf("\n\n");
	return 1;
    }
}
    printf("bls12 g2 scm.                                                : %.4f[ms]\n",scm_time/scm);
    printf("bls12 g2 scm 2split.                                         : %.4f[ms]\n",scm_2split_time/scm);
    printf("bls12 g2 scm 2split jsf.                                     : %.4f[ms]\n",scm_2split_jsf_time/scm);
    printf("bls12 g2 scm 4split.                                         : %.4f[ms]\n",scm_4split_time/scm);
    printf("bls12 g2 scm 4split lazy.                                    : %.4f[ms]\n",scm_lazy_time/scm);
    printf("bls12 g2 scm 4split jacobian lazy.                           : %.4f[ms]\n",scm_jacobian_time/scm);
    printf("bls12 g2 scm 4split mixture lazy.                            : %.4f[ms]\n",scm_mixture_time/scm);
    printf("bls12 g2 scm 4split 3naf interleave mixture.                 : %.4f[ms]\n",scm_2split_3naf_I_mixture_time/scm);
    printf("bls12 g2 scm 4split 5naf interleave mixture.                 : %.4f[ms]\n",scm_2split_5naf_I_mixture_time/scm);
    printf("bls12 g2 scm 4split 5naf interleave mixture montgomery.      : %.4f[ms]\n",scm_2split_5naf_I_mixture_montgomery_time/scm);

    #ifdef DEBUG_COST_A
    printf("*********bls12 g2 scm fp COST.********         \n");
    cost_printf("bls12 g2 scm",&scm_cost,scm);
    cost_printf("bls12 g2 scm 4split 5naf interleave mixture montgomery",&scm_2split_5naf_I_mixture_montgomery_cost,scm);
    printf("***************************************         \n");
    #endif
    
    mpz_clear(scalar);

    return 0;
}



int bls12_test_g3_exp(int exp){
    int i,n=0;
    float exp_time=0,exp_2split_time=0,exp_2split_jsf_time=0,exp_4split_time=0,exp_lazy_time=0,exp_gs_time=0,exp_gs_lazy_time=0,exp_5naf_gs_lazy_time=0,exp_5naf_gs_lazy_montgomery_time=0;
    cost tmp,exp_cost,exp_2split_cost,exp_2split_jsf_cost,exp_4split_cost,exp_lazy_cost,exp_gs_cost,exp_gs_lazy_cost,exp_5naf_gs_lazy_cost,exp_5naf_gs_lazy_montgomery_cost;
    struct timeval tv_A,tv_B;
    printf("================================================================================\n");
    printf("G3 Exp.\n\n");
    efp12_t P,Q;
    fp12_t A_fp12,test1,test2,test3,test4,test5,test6,test7,test8,test9;
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
    fp12_init(&test8);
    fp12_init(&test9);
    mpz_t scalar;
    mpz_init(scalar);
    
    cost_init(&tmp);
    cost_init(&exp_cost);
    cost_init(&exp_5naf_gs_lazy_montgomery_cost);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
for(i=0;i<exp;i++){

    mpz_urandomm(scalar,state,order_z);
    bls12_generate_g1(&P);
    bls12_generate_g2(&Q);
    bls12_optate_pairing_basic(&A_fp12,&P,&Q);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_basic(&test1,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&exp_cost,&tmp);

    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_2split(&test2,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_2split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_2split_jsf(&test3,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_2split_jsf_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_4split(&test4,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_4split_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_4split_lazy(&test5,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_4split_GS(&test6,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_gs_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_4split_GS_lazy(&test7,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_gs_lazy_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_4split_5naf_interleaving_GS_lazy(&test8,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_5naf_gs_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g3_exp_4split_5naf_interleaving_GS_lazy_montgomery(&test9,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_5naf_gs_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&exp_5naf_gs_lazy_montgomery_cost,&tmp);
    
    if(fp12_cmp(&test1,&test2)!=0 || fp12_cmp(&test1,&test3)!=0 || fp12_cmp(&test1,&test4)!=0 || fp12_cmp(&test1,&test5)!=0 || fp12_cmp(&test1,&test6)!=0 || fp12_cmp(&test1,&test7)!=0 || fp12_cmp(&test1,&test8)!=0 || fp12_cmp(&test1,&test9)!=0){
        printf("failed!\n\n");
	fp12_printf("test1=",&test1);
	fp12_printf("\ntest2=",&test2);
	fp12_printf("\ntest3=",&test3);
	fp12_printf("\ntest4=",&test4);
	fp12_printf("\ntest5=",&test5);
	fp12_printf("\ntest6=",&test6);
	fp12_printf("\ntest7=",&test7);
	fp12_printf("\ntest8=",&test8);
	fp12_printf("\ntest9=",&test9);
	printf("\n\n");
	return 1;
    }
}
    printf("bls12 G3 exp.                           : %.4f[ms]\n",exp_time/exp);
    printf("bls12 G3 exp 2split.                    : %.4f[ms]\n",exp_2split_time/exp);
    printf("bls12 G3 exp 2split jsf.                : %.4f[ms]\n",exp_2split_jsf_time/exp);
    printf("bls12 G3 exp 4split.                    : %.4f[ms]\n",exp_4split_time/exp);
    printf("bls12 G3 exp 4split lazy.               : %.4f[ms]\n",exp_lazy_time/exp);
    printf("bls12 G3 exp GS.                        : %.4f[ms]\n",exp_gs_time/exp);
    printf("bls12 G3 exp GS lazy.                   : %.4f[ms]\n",exp_gs_lazy_time/exp);
    printf("bls12 G3 exp 5naf GS lazy.              : %.4f[ms]\n",exp_5naf_gs_lazy_time/exp);
    printf("bls12 G3 exp 5naf GS lazy montgomery.   : %.4f[ms]\n",exp_5naf_gs_lazy_montgomery_time/exp);
    
    #ifdef DEBUG_COST_A
    printf("*********bls12 G3 exp fp COST.********         \n");
    cost_printf("bls12 G3 exp",&exp_cost,exp);
    cost_printf("bls12 G3 exp 5naf GS lazy montgomery",&exp_5naf_gs_lazy_montgomery_cost,exp);
    printf("***************************************         \n");
    #endif

    mpz_clear(scalar);

    return 0;
}
