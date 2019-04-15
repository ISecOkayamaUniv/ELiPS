#include <ELiPS/bn12_test.h>
void BN12_test_rational_point(){
    printf("====================================================================================\n");
    EFp12 test_G1,test_G2;
    EFp12_init(&test_G1);
    EFp12_init(&test_G2);
    
    BN12_EFp12_generate_G1(&test_G1);
    EFp12_printf(&test_G1,"G1\n");
    printf("\n");
    EFp12_SCM(&test_G1,&test_G1,order_z);
    EFp12_printf(&test_G1,"G1 test\n");
    printf("\n");
    
    EFp12_generate_G2(&test_G2);
    EFp12_printf(&test_G2,"G2\n");
    printf("\n");
    EFp12_SCM(&test_G2,&test_G2,order_z);
    EFp12_printf(&test_G2,"G2 test\n");
    printf("\n");
}
/*----------------------------------------------------------------------------*/
//bn12_pairing
void BN12_test_plain_ate_pairing(){
    printf("====================================================================================\n");
    printf("Plain-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
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
    
    BN12_EFp12_generate_G1(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM(&s1P,&P,s1);
    EFp12_SCM(&s2P,&P,s2);
    EFp12_SCM(&s1Q,&Q,s1);
    EFp12_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain_ate(Q,P)^s1*s2\n");
    BN12_Plain_ate_pairing(&Z,&Q,&P);
    Fp12_pow(&test1,&Z,s12);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_PLAINATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("plain_ate([s2]Q,[s1]P)\n");
    BN12_Plain_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_PLAINATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("plain_ate([s1]Q,[s2]P)\n");
    BN12_Plain_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_PLAINATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
}

void BN12_test_opt_ate_pairing(){
    printf("====================================================================================\n");
    printf("Opt-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
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
    
    BN12_EFp12_generate_G1(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM(&s1P,&P,s1);
    EFp12_SCM(&s2P,&P,s2);
    EFp12_SCM(&s1Q,&Q,s1);
    EFp12_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("opt_ate(Q,P)^s1*s2\n");
    BN12_Opt_ate_pairing(&Z,&Q,&P);
    Fp12_pow(&test1,&Z,s12);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_OPTATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("opt_ate([s2]Q,[s1]P)\n");
    BN12_Opt_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_OPTATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("opt_ate([s1]Q,[s2]P)\n");
    BN12_Opt_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_OPTATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
}

void BN12_test_x_ate_pairing(){
    printf("====================================================================================\n");
    printf("X-ate pairing\n\n");
    EFp12 P,Q,s1P,s2P,s1Q,s2Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1P);
    EFp12_init(&s2P);
    EFp12_init(&s1Q);
    EFp12_init(&s2Q);
    Fp12 Z,test1,test2,test3;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
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
    
    BN12_EFp12_generate_G1(&P);
    EFp12_printf(&P,"P\n");
    printf("\n\n");
    EFp12_generate_G2(&Q);
    EFp12_printf(&Q,"Q\n");
    printf("\n\n");
    EFp12_SCM(&s1P,&P,s1);
    EFp12_SCM(&s2P,&P,s2);
    EFp12_SCM(&s1Q,&Q,s1);
    EFp12_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("x_ate(Q,P)^s1*s2\n");
    BN12_X_ate_pairing(&Z,&Q,&P);
    Fp12_pow(&test1,&Z,s12);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_XATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("x_ate([s2]Q,[s1]P)\n");
    BN12_X_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_XATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("x_ate([s1]Q,[s2]P)\n");
    BN12_X_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (X-ate) : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Final Exp. (total) : %.2f[ms]\n",MILLER_XATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
}
/*----------------------------------------------------------------------------*/
//BN12_SCM
void BN12_test_G1_SCM(){
    printf("====================================================================================\n");
    printf("G1 SCM\n\n");
    EFp12 P,test1,test2,test3;
    EFp12_init(&P);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp12_init(&test3);
    mpz_t scalar;
    mpz_init(scalar);
    
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    //mpz_urandomm(scalar,state,order);
    
    mpz_set_str(scalar,"6433987872172996767370742739789192855206084819751810317146156999801714004530205817441825065255166612237926117694268872766909697908042480267",10);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    BN12_EFp12_generate_G1(&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp12_G1_SCM_plain(&test1,&P,scalar);
    printf("G1 SCM (plain) : %.2f[ms]\n",G1SCM_PLAIN);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    BN12_EFp12_G1_SCM_2split(&test2,&P,scalar);
    printf("G1 SCM (2split) : %.2f[ms]\n",G1SCM_2SPLIT);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    BN12_EFp12_G1_SCM_2split_JSF(&test3,&P,scalar);
    printf("G1 SCM (2split-JSF) : %.2f[ms]\n",G1SCM_2SPLIT_JSF);
    EFp12_printf(&test3,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0
    && Fp12_cmp(&test1.x,&test3.x)==0 && Fp12_cmp(&test1.y,&test3.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
}

void BN12_test_G2_SCM(){
    printf("====================================================================================\n");
    printf("G2 SCM\n\n");
    EFp12 Q,test1,test2,test3,test4;
    EFp12_init(&Q);
    EFp12_init(&test1);
    EFp12_init(&test2);
    EFp12_init(&test3);
    EFp12_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    //mpz_urandomm(scalar,state,order);
    
    mpz_set_str(scalar,"6433987872172996767370742739789192855206084819751810317146156999801714004530205817441825065255166612237926117694268872766909697908042480267",10);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2(&Q);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp12_G2_SCM_plain(&test1,&Q,scalar);
    printf("G2 SCM (plain) : %.2f[ms]\n",G2SCM_PLAIN);
    EFp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    EFp12_G2_SCM_2split(&test2,&Q,scalar);
    printf("G2 SCM (2split) : %.2f[ms]\n",G2SCM_2SPLIT);
    EFp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    EFp12_G2_SCM_2split_JSF(&test3,&Q,scalar);
    printf("G2 SCM (2split-JSF) : %.2f[ms]\n",G2SCM_2SPLIT_JSF);
    EFp12_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("4split\n");
    BN12_EFp12_G2_SCM_4split(&test4,&Q,scalar);
    printf("G2 SCM (4split) : %.2f[ms]\n",G2SCM_4SPLIT);
    EFp12_printf(&test4,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0
    && Fp12_cmp(&test1.x,&test3.x)==0 && Fp12_cmp(&test1.y,&test3.y)==0
    && Fp12_cmp(&test1.x,&test4.x)==0 && Fp12_cmp(&test1.y,&test4.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
}
void BN12_test_G3_EXP(){
    printf("====================================================================================\n");
    printf("G3 Exp.\n\n");
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12 Z,test1,test2,test3,test4;
    Fp12_init(&Z);
    Fp12_init(&test1);
    Fp12_init(&test2);
    Fp12_init(&test3);
    Fp12_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    //mpz_urandomm(scalar,state,order);
    
    mpz_set_str(scalar,"6433987872172996767370742739789192855206084819751810317146156999801714004530205817441825065255166612237926117694268872766909697908042480267",10);
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    BN12_EFp12_generate_G1(&P);
    printf("generating rational point Q in G2\n\n");
    EFp12_generate_G2(&Q);
    printf("x-ate(Q,P)\n");
    BN12_Opt_ate_pairing(&Z,&Q,&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    Fp12_G3_EXP_plain(&test1,&Z,scalar);
    printf("G3 exp (plain) : %.2f[ms]\n",G3EXP_PLAIN);
    Fp12_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    Fp12_G3_EXP_2split(&test2,&Z,scalar);
    printf("G3 exp (2split) : %.2f[ms]\n",G3EXP_2SPLIT);
    Fp12_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    Fp12_G3_EXP_2split_JSF(&test3,&Z,scalar);
    printf("G3 exp (2split-JSF) : %.2f[ms]\n",G3EXP_2SPLIT_JSF);
    Fp12_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("4split\n");
    BN12_Fp12_G3_EXP_4split(&test4,&Z,scalar);
    printf("G3 exp (4split) : %.2f[ms]\n",G3EXP_4SPLIT);
    Fp12_printf(&test4,"");
    printf("\n\n");
    
    if(Fp12_cmp(&test1,&test2)==0 && Fp12_cmp(&test1,&test3)==0 && Fp12_cmp(&test1,&test4)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
}

/*----------------------------------------------------------------------------*/
//BN12_compare
void BN12_compare_pairings(){
    printf("====================================================================================\n");
    printf("Ate-based pairing\n\n");
    EFp12 P,Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    Fp12 Z;
    Fp12_init(&Z);
    EFp2 twisted_Q;
    EFp2_init(&twisted_Q);
    
    printf("generating rational point P in G1...\n\n");
    BN12_EFp12_generate_G1(&P);
    //EFp12_printf(&P,"P\n");
    //printf("\n\n");
    printf("generating rational point Q in G2...\n\n");
    EFp12_generate_G2(&Q);
    //EFp12_printf(&Q,"Q\n");
    //printf("\n\n");
    
    /*printf("------------------------------------------------------------------------------------\n");
    printf("Plain-ate pairing\n\n");
    Plain_ate_pairing(&Z,&Q,&P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Pairing (total) : %.2f[ms]\n",MILLER_PLAINATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    printf("\n");*/
    
    printf("------------------------------------------------------------------------------------\n");
    printf("Opt-ate pairing\n\n");
    BN12_Opt_ate_pairing(&Z,&Q,&P);
    printf("Miller's Algo. (Opt-ate)   : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Pairing (total) : %.2f[ms]\n",MILLER_OPTATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    printf("\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("X-ate pairing\n\n");
    BN12_X_ate_pairing(&Z,&Q,&P);
    printf("Miller's Algo. (X-ate)     : %.2f[ms]\n",MILLER_XATE);
    printf("Final Exp. (easy)  : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. (hard)  : %.2f[ms]\n",FINALEXP_OPT_HARD);
    printf("Pairing (total) : %.2f[ms]\n",MILLER_XATE+FINALEXP_OPT_EASY+FINALEXP_OPT_HARD);
    printf("\n");
    
    /*printf("------------------------------------------------------------------------------------\n");
    printf("Sextic twist\n\n");
    gettimeofday(&tv_start,NULL);
    EFp12_to_EFp2(&twisted_Q,&Q);
    gettimeofday(&tv_end,NULL);
    printf("EFp12 to EFp2     : %.2f[us]\n",timedifference_usec(tv_start,tv_end));
    
    gettimeofday(&tv_start,NULL);
    EFp2_to_EFp12(&Q,&twisted_Q);
    gettimeofday(&tv_end,NULL);
    printf("EFp2 to EFp12     : %.2f[us]\n",timedifference_usec(tv_start,tv_end));
    */
    
}
