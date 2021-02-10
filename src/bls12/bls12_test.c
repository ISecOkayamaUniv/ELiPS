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


int bls12_test_opt_ate_pairing(int pairing){
    int i,n=0;
    float opt_time=0,opt_affine_time=0;
    cost tmp,opt_cost,opt_affine_cost;
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

    fp12_t Z,testA,testB,testC,test1,test2,test3;
    fp12_init(&Z);
    fp12_init(&testA);
    fp12_init(&testB);
    fp12_init(&testC);
    fp12_init(&test1);
    fp12_init(&test2);
    fp12_init(&test3);

    cost_init(&tmp);
    cost_init(&opt_cost);
    cost_init(&opt_affine_cost);

    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);

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

    bls12_optate_pairing_basic(&Z,&P,&Q);
    fp12_pow(&testA,&Z,s12);
    bls12_optate_pairing_basic(&testB,&s1P,&s2Q);
    bls12_optate_pairing_basic(&testC,&s2P,&s1Q);

    printf("bilinear test\n");
    fp12_printf("testA=",&testA);
    if(fp12_cmp(&testA,&testB)!=0 || fp12_cmp(&testA,&testC)!=0){
        printf("bilinear failed!!\n\n");
	return 1;
    }
    printf("bilinear succeced!!\n\n");

    MILLER_OPT_AFFINE = 0;
    FINALEXP_OPT_AFFINE = 0;
    MILLER_OPT_PROJECTIVE = 0;
    FINALEXP_OPT_PROJECTIVE = 0;
    cost_init(&MILLER_OPT_AFFINE_COST);
    cost_init(&FINALEXP_OPT_AFFINE_COST);
    cost_init(&MILLER_OPT_PROJECTIVE_COST);
    cost_init(&FINALEXP_OPT_PROJECTIVE_COST);

    bls12_generate_g2(&Q);

for(i=0;i<pairing;i++){

    bls12_generate_g1(&P);

    bls12_optate_pairing_basic(&test1,&P,&Q);

    // cost_zero();
    // gettimeofday(&tv_A,NULL);
    // bls12_optate_pairing_affine(&test2,&P,&Q);
    // gettimeofday(&tv_B,NULL);
    // opt_affine_time+=timedifference_msec(tv_A,tv_B);
    // cost_check(&tmp);
    // cost_addition(&opt_affine_cost,&tmp);


    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_optate_pairing(&test3,&P,&Q);
    gettimeofday(&tv_B,NULL);
    opt_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&opt_cost,&tmp);

    // if(fp12_cmp(&test1,&test2) != 0){
    //     printf("pairing affine failed!\n\n");
    //     fp12_printf("test2=",&test2);
	//     printf("\n\n");
	//     return 1;
    // }
    if(fp12_cmp(&test1,&test3) != 0){
        printf("pairing projective failed!\n\n");
	    printf("\n\n");
	    return 1;
    }
}
    //cost_substruction(&FINALEXP_OPT_AFFINE_COST, &opt_affine_cost, &MILLER_OPT_AFFINE_COST);
    cost_substruction(&FINALEXP_OPT_PROJECTIVE_COST, &opt_cost, &MILLER_OPT_PROJECTIVE_COST);

    // printf("bls12 opt ate affine.            : %.4f[ms]\n",opt_affine_time/pairing);
    // printf("bls12 opt ate affine (MILLER).   : %.4f[ms]\n",MILLER_OPT_AFFINE/pairing);
    // printf("bls12 opt ate affine (FINALEXP). : %.4f[ms]\n",FINALEXP_OPT_AFFINE/pairing);

    printf("bls12 opt ate.                   : %.4f[ms]\n",opt_time/pairing);
    printf("bls12 opt ate (MILLER).          : %.4f[ms]\n",MILLER_OPT_PROJECTIVE/pairing);
    printf("bls12 opt ate (FINALEXP).        : %.4f[ms]\n",FINALEXP_OPT_PROJECTIVE/pairing);

    #ifdef DEBUG_COST_A
    printf("*********bls12 opt ate fp COST.********         \n");
    // cost_printf("bls12 opt ate affine", &opt_affine_cost, pairing);
    // cost_printf("bls12 opt ate affine (MILLER)", &MILLER_OPT_AFFINE_COST, pairing);
    // cost_printf("bls12 opt ate affine (FINALEXP)", &FINALEXP_OPT_AFFINE_COST, pairing);


    cost_printf("bls12 opt ate ", &opt_cost, pairing);
    cost_printf("bls12 opt ate (MILLER)", &MILLER_OPT_PROJECTIVE_COST, pairing);
    cost_printf("bls12 opt ate (FINALEXP)", &FINALEXP_OPT_PROJECTIVE_COST, pairing);
    printf("***************************************         \n");
    #endif

    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);

    return 0;
}
int bls12_test_symmmetric_opt_ate_pairing(){
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("bls12_symmetric_Opt-ate pairing\n\n");

    sym_t A,B,s1A,s2A,s1B,s2B;
    sym_init(&A);
    sym_init(&B);
    sym_init(&s1A);
    sym_init(&s2A);
    sym_init(&s1B);
    sym_init(&s2B);

    fp12_t test1,test2,test3;
    fp12_init(&test1);
    fp12_init(&test2);
    fp12_init(&test3);
    
    mpz_t a,s12,s1,s2;
    mpz_init(a);
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    mpz_urandomm(a,state,order_z);
    mpz_urandomm(s1,state,order_z);
    mpz_urandomm(s2,state,order_z);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order_z);

    bls12_generate_symmetric_point(&A,a);
    bls12_generate_symmetric_point(&B,a);

    bls12_sym_scm(&s1A,&A,s1);
    bls12_sym_scm(&s2A,&A,s2);
    bls12_sym_scm(&s1B,&B,s1);
    bls12_sym_scm(&s2B,&B,s2);
        
    bls12_symmetric_optate_pairing(&test1,&A,&B);
    fp12_pow(&test1,&test1,s12);
    bls12_symmetric_optate_pairing(&test2,&s1A,&s2B);    
    bls12_symmetric_optate_pairing(&test3,&s2A,&s1B);
    
    printf("symmmetric pairing bilinear test\n");
    if(fp12_cmp(&test1,&test2)!=0 || fp12_cmp(&test1,&test3)!=0){
        printf("bilinear failed!!\n\n");
	return 1;
    }
    printf("bilinear succeced!!\n\n");
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);

    return 0;
}
/*----------------------------------------------------------------------------*/
//bls12_scm
int bls12_test_g1_scm(int scm){
    int i,n=0;
    float scm_time=0;
    cost tmp,scm_cost;
    struct timeval tv_A,tv_B;

    efp12_t A_efp12,test1,test2;
    efp12_init(&A_efp12);
    efp12_init(&test1);
    efp12_init(&test2);

    mpz_t scalar;
    mpz_init(scalar);

    cost_init(&tmp);
    cost_init(&scm_cost);

    printf("================================================================================\n");
    printf("g1 scm test\n\n");

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    bls12_generate_g1(&A_efp12);

    //basic type
    bls12_g1_scm_basic(&test1,&A_efp12,scalar);

    //faster type
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g1_scm(&test2,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&scm_cost,&tmp);

    if(efp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    efp12_printf("test1=",&test1);
	    efp12_printf("\ntest2=",&test2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("bls12 g1 scm.     : %.4f[ms]\n",scm_time/scm);

    #ifdef DEBUG_COST_A
    printf("*********bls12 g2 scm fp COST.********         \n");
    cost_printf("bls12 g1 scm",&scm_cost,scm);
    printf("***************************************         \n");
    #endif

    mpz_clear(scalar);

    return 0;
}

int bls12_test_g2_scm(int scm){
    int i,n=0;
    float scm_time=0;
    cost tmp,scm_cost;
    struct timeval tv_A,tv_B;

    efp12_t A_efp12,test1,test2;
    efp12_init(&A_efp12);
    efp12_init(&test1);
    efp12_init(&test2);

    mpz_t scalar;
    mpz_init(scalar);

    cost_init(&tmp);
    cost_init(&scm_cost);

    printf("================================================================================\n");
    printf("g2 scm\n\n");

    gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    gmp_randseed_ui(state,1);
    mpz_urandomm(scalar,state,order_z);
    bls12_generate_g2(&A_efp12);

for(i=0;i<scm;i++){
    

    //basic type
    bls12_g2_scm_basic(&test1,&A_efp12,scalar);

    //faster type
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g2_scm(&test2,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&scm_cost,&tmp);

    if(efp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    efp12_printf("test1=",&test1);
	    efp12_printf("\ntest2=",&test2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("bls12 g2 scm.     : %.4f[ms]\n",scm_time/scm);

    #ifdef DEBUG_COST_A
    printf("*********bls12 g2 scm fp COST.********         \n");
    cost_printf("bls12 g2 scm",&scm_cost,scm);
    printf("***************************************         \n");
    #endif

    mpz_clear(scalar);

    return 0;
}



int bls12_test_g3_exp(int exp){
    int i,n=0;
    float exp_time=0;
    cost tmp,exp_cost;
    struct timeval tv_A,tv_B;

    efp12_t P,Q;
    fp12_t A_fp12,test1,test2;
    efp12_init(&P);
    efp12_init(&Q);
    fp12_init(&A_fp12);
    fp12_init(&test1);
    fp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);

    cost_init(&tmp);
    cost_init(&exp_cost);

    printf("================================================================================\n");
    printf("G3 Exp.\n\n");

    bls12_generate_g2(&Q);

for(i=0;i<exp;i++){

    mpz_urandomm(scalar,state,order_z);
    bls12_generate_g1(&P);
    bls12_optate_pairing(&A_fp12,&P,&Q);

    //basic type
    bls12_g3_exp_basic(&test1,&A_fp12,scalar);

    //faster type
    cost_zero();
    gettimeofday(&tv_A,NULL);
    bls12_g3_exp(&test2,&A_fp12,scalar);
    gettimeofday(&tv_B,NULL);
    exp_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&exp_cost,&tmp);

    if(fp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    fp12_printf("test1=",&test1);
	    fp12_printf("\ntest2=",&test2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("bls12 G3 exp.     : %.4f[ms]\n",exp_time/exp);

    #ifdef DEBUG_COST_A
    printf("*********bls12 G3 exp fp COST.********         \n");
    cost_printf("bls12 G3 exp",&exp_cost,exp);
    printf("***************************************         \n");
    #endif

    mpz_clear(scalar);

    return 0;
}
