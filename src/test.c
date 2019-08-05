#include <ELiPS/test.h>
/*----------------------------------------------------------------------------*/
//test
int test_Field(int fp,int fp2,int fp6,int fp12,int sqr){
    int i,j,n=0;
    float add_time=0,add_lazy_time=0;
    float shift2_time=0,shift4_time=0,shift8_time=0;
    float mul_time=0,mul_lazy_time=0;
    float inv_time=0,inv_lazy_time=0;
    float sqr_time=0,sqr_lazy_time=0;
    float com_time=0,rec_time=0,gs_time=0;
    float com_lazy_time=0,rec_lazy_time=0,gs_lazy_time=0;
    EFp12 P,Q;
    Fp12 tmp,t0,t1,f;
    struct timeval tv_A,tv_B;
    
    printf("====================================================================================\n");
    printf("Field test\n");
    mp_limb_t test1_mpn[FPLIMB2],test2_mpn[FPLIMB2];
    Fp A_Fp,B_Fp,test1_Fp,test2_Fp;
    Fp2 A_Fp2,B_Fp2,test1_Fp2,test2_Fp2;
    Fp6 A_Fp6,B_Fp6,test1_Fp6,test2_Fp6,test3_Fp6;
    Fp12 A_Fp12,B_Fp12,test1_Fp12,test2_Fp12,test3_Fp12,test4_Fp12;
    Fp12 test0_Fp12,test1l_Fp12,test2l_Fp12,test3l_Fp12;

    Fp_init(&A_Fp);
    Fp_init(&B_Fp);
    Fp_init(&test1_Fp);
    Fp_init(&test2_Fp);
    
    Fp2_init(&A_Fp2);
    Fp2_init(&B_Fp2);
    Fp2_init(&test1_Fp2);
    Fp2_init(&test2_Fp2);

    
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
    printf("Fp_add test\n");
add_time=0,add_lazy_time=0;
n=1000;
for(i=0;i<fp;i++){
    Fp_set_random(&A_Fp,state);
    Fp_set_random(&B_Fp,state);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp_add(&test1_Fp,&A_Fp,&B_Fp);
    gettimeofday(&tv_B,NULL);
    add_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp_add_lazy(test2_mpn,FPLIMB,A_Fp.x0,FPLIMB,B_Fp.x0,FPLIMB);
    gettimeofday(&tv_B,NULL);
    add_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    Fp_mod(&test2_Fp,test2_mpn,FPLIMB);

    if(Fp_cmp(&test1_Fp,&test2_Fp)!=0){
        printf("failed!\n\n");
    	Fp_printf("",&test1_Fp);
	    Fp_printf("\n",&test2_Fp);
	    printf("\n\n");
	    return 1;
    }
}
    printf("Fp add.      : %.6f[ms]\n",add_time/fp);
    printf("Fp add lazy. : %.6f[ms]\n",add_lazy_time/fp);
    
printf("------------------------------------------------------------------------------------\n");
    printf("Fp_mul test\n");
mul_time=0,mul_lazy_time=0;
n=1000;
for(i=0;i<fp;i++){
    Fp_set_random(&A_Fp,state);
    Fp_set_random(&B_Fp,state);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp_mul(&test1_Fp,&A_Fp,&B_Fp);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp_mul_lazy(test2_mpn,A_Fp.x0,B_Fp.x0);
    gettimeofday(&tv_B,NULL);
    mul_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    Fp_mod(&test2_Fp,test2_mpn,FPLIMB2);

    if(Fp_cmp(&test1_Fp,&test2_Fp)!=0){
        printf("failed!\n\n");
    	Fp_printf("",&test1_Fp);
	    Fp_printf("\n",&test2_Fp);
	    printf("\n\n");
	    return 1;
    }
}
    printf("Fp mul.      : %.6f[ms]\n",mul_time/fp);
    printf("Fp mul lazy. : %.6f[ms]\n",mul_lazy_time/fp);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("Fp_sqr test\n");
sqr_time=0,sqr_lazy_time=0;
n=1000;
for(i=0;i<fp;i++){
    Fp_set_random(&A_Fp,state);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp_sqr(&test1_Fp,&A_Fp);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp_sqr_lazy(test2_mpn,A_Fp.x0);
    gettimeofday(&tv_B,NULL);
    sqr_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    Fp_mod(&test2_Fp,test2_mpn,FPLIMB2);

    if(Fp_cmp(&test1_Fp,&test2_Fp)!=0){
        printf("failed!\n\n");
    	Fp_printf("",&test1_Fp);
	    Fp_printf("\n",&test2_Fp);
	    printf("\n\n");
	    return 1;
    }
}
    printf("Fp sqr.      : %.6f[ms]\n",sqr_time/fp);
    printf("Fp sqr lazy. : %.6f[ms]\n",sqr_lazy_time/fp);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("Fp_inv test\n");
inv_time=0,inv_lazy_time=0;
n=1000;
for(i=0;i<fp;i++){
    Fp_set_random(&A_Fp,state);
    Fp_set_random(&B_Fp,state);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp_inv(&test1_Fp,&A_Fp);
    gettimeofday(&tv_B,NULL);
    inv_time+=timedifference_msec(tv_A,tv_B)/n;

}
    printf("Fp inv.      : %.6f[ms]\n",inv_time/fp);
    
printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_add test\n");
add_time=0,add_lazy_time=0;
for(i=0;i<fp2;i++){
    Fp2_set_random(&A_Fp2,state);
    Fp2_set_random(&B_Fp2,state);
    
    gettimeofday(&tv_A,NULL);
    Fp2_add(&test1_Fp2,&A_Fp2,&B_Fp2);
    gettimeofday(&tv_B,NULL);
    add_time+=timedifference_msec(tv_A,tv_B);

}
    printf("Fp2 add.      : %.6f[ms]\n",add_time/fp2);
    
printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_shift test\n");
    shift2_time=0,shift4_time=0,shift8_time=0;
for(i=0;i<fp2;i++){
    Fp2_set_random(&A_Fp2,state);
    Fp2_set_random(&B_Fp2,state);

    gettimeofday(&tv_A,NULL);
    Fp2_lshift2(&test1_Fp2,&A_Fp2);
    gettimeofday(&tv_B,NULL);
    shift2_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    Fp2_lshift(&test1_Fp2,&A_Fp2,2);
    gettimeofday(&tv_B,NULL);
    shift4_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    Fp2_lshift(&test1_Fp2,&A_Fp2,3);
    gettimeofday(&tv_B,NULL);
    shift8_time+=timedifference_msec(tv_A,tv_B);
    
    }

    printf("Fp2 shift2.      : %.6f[ms]\n",shift2_time/fp2);
    printf("Fp2 shift4.      : %.6f[ms]\n",shift4_time/fp2);
    printf("Fp2 shift8.      : %.6f[ms]\n",shift8_time/fp2);
    
    
printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_mul test\n");
mul_time=0,mul_lazy_time=0;
n=100;
for(i=0;i<fp2;i++){
    Fp2_set_random(&A_Fp2,state);
    Fp2_set_random(&B_Fp2,state);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp2_mul(&test1_Fp2,&A_Fp2,&B_Fp2);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp2_mul_lazy(&test2_Fp2,&A_Fp2,&B_Fp2);
    gettimeofday(&tv_B,NULL);
    mul_lazy_time+=timedifference_msec(tv_A,tv_B)/n;


    if(Fp2_cmp(&test1_Fp2,&test2_Fp2)!=0){
        printf("failed!\n\n");
	    Fp2_printf("",&test1_Fp2);
	    Fp2_printf("\n",&test2_Fp2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("Fp2 mul.      : %.6f[ms]\n",mul_time/fp2);
    printf("Fp2 mul lazy. : %.6f[ms]\n",mul_lazy_time/fp2);


printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_sqr test\n");
sqr_time=0,sqr_lazy_time=0;
n=100;
for(i=0;i<fp2;i++){
    Fp2_set_random(&A_Fp2,state);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp2_sqr(&test1_Fp2,&A_Fp2);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp2_sqr_lazy(&test2_Fp2,&A_Fp2);
    gettimeofday(&tv_B,NULL);
    sqr_lazy_time+=timedifference_msec(tv_A,tv_B)/n;


    if(Fp2_cmp(&test1_Fp2,&test2_Fp2)!=0){
        printf("failed!\n\n");
	    Fp2_printf("",&test1_Fp2);
	    Fp2_printf("\n",&test2_Fp2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("Fp2 sqr.      : %.6f[ms]\n",sqr_time/fp2);
    printf("Fp2 sqr lazy. : %.6f[ms]\n",sqr_lazy_time/fp2);
    
printf("------------------------------------------------------------------------------------\n");
    printf("Fp2_inv test\n");
inv_time=0,inv_lazy_time=0;
n=1000;
for(i=0;i<fp2;i++){
    Fp2_set_random(&A_Fp2,state);
    Fp2_set_random(&B_Fp2,state);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp2_inv(&test1_Fp2,&A_Fp2);
    gettimeofday(&tv_B,NULL);
    inv_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp2_inv_lazy(&test1_Fp2,&A_Fp2);
    gettimeofday(&tv_B,NULL);
    inv_lazy_time+=timedifference_msec(tv_A,tv_B)/n;

}
    printf("Fp2 inv.      : %.6f[ms]\n",inv_time/fp2);
    printf("Fp2 inv_lazy. : %.6f[ms]\n",inv_lazy_time/fp2);

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
	    Fp6_printf("",&test1_Fp6);
	    Fp6_printf("\n",&test2_Fp6);
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
	    Fp6_printf("",&test1_Fp6);
	    Fp6_printf("\n",&test2_Fp6);
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
	    Fp12_printf("",&test3_Fp12);
	    Fp12_printf("\n",&test4_Fp12);
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
	    Fp12_printf("",&test1_Fp12);
	    Fp12_printf("\n",&test2_Fp12);
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
	    Fp12_printf("",&test1_Fp12);
	    Fp12_printf("\n",&test2_Fp12);
	    printf("\n\n");
	    return 1;
    }
}
    printf("Fp12 sqr_cyclotomic.      : %.4f[ms]\n",sqr_time/fp12);
    printf("Fp12 sqr_cyclotomic lazy. : %.4f[ms]\n",sqr_lazy_time/fp12);

printf("------------------------------------------------------------------------------------\n");
    printf("Fp12_sqr_compressed test\n");
n=10;
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
    for(j=0;j<n;j++)Fp12_sqr_compressed(&test1_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    com_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp12_sqr_compressed_lazy(&test1l_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    com_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp12_sqr_recover_g1(&test2_Fp12,&test1_Fp12);
    for(j=0;j<n;j++)Fp12_sqr_recover_g0(&test2_Fp12,&test2_Fp12);
    gettimeofday(&tv_B,NULL);
    rec_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp12_sqr_recover_g1_lazy(&test2l_Fp12,&test1l_Fp12);
    for(j=0;j<n;j++)Fp12_sqr_recover_g0_lazy(&test2l_Fp12,&test2l_Fp12);
    gettimeofday(&tv_B,NULL);
    rec_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp12_sqr_GS(&test3_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    gs_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)Fp12_sqr_GS_lazy(&test3l_Fp12,&tmp);
    gettimeofday(&tv_B,NULL);
    gs_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    
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

void test_Fp(int fp){
    int i,j,n=10000;
    float add_time=0,dbl_time=0,mul_time=0,sqr_time=0,mod_time=0,mod_mont_time=0;
    struct timespec tv_a,tv_b;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("Field Fp test\n");
    mp_limb_t A[FPLIMB],B[FPLIMB],C[FPLIMB],C1[FPLIMB2];
    mp_limb_t At[FPLIMB2],Bt[FPLIMB2],Ct[FPLIMB2],test_mul[FPLIMB2];
    Fp test1,test2;

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<fp;i++){
    mpn_random(A,FPLIMB);
    mpn_random(B,FPLIMB);

    mpn_mod(A,A,FPLIMB);
    mpn_mod(B,B,FPLIMB);

    mpn_mul_n(At,A,A,FPLIMB);
    mpn_mul_n(Bt,B,B,FPLIMB);

    clock_gettime(CLOCK_REALTIME, &tv_a);
    for(j=0;j<n;j++)    mpn_add_n(test_mul,At,At,FPLIMB2);
    clock_gettime(CLOCK_REALTIME, &tv_b);
    add_time+=timedifference_nsec(tv_a,tv_b);

    clock_gettime(CLOCK_REALTIME, &tv_a);
    for(j=0;j<n;j++)    mpn_lshift(test_mul,At,FPLIMB2,1);
    clock_gettime(CLOCK_REALTIME, &tv_b);
    dbl_time+=timedifference_nsec(tv_a,tv_b);


    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    mpn_sqr(test_mul,A,FPLIMB);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    mpn_mul_n(test_mul,A,A,FPLIMB);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    Fp_mod(&test1,test_mul,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    mod_time+=timedifference_msec(tv_A,tv_B);
    }
    add_time=add_time/(1000*1000);
    dbl_time=dbl_time/(1000*1000);
    printf("mpn add.            : %.6f[ms]\n",add_time/(fp*n));
    printf("mpn dbl.            : %.6f[ms]\n",dbl_time/(fp*n));
    printf("mpn sqr.            : %.6f[ms]\n",sqr_time/(fp*n));
    printf("mpn mul.            : %.6f[ms]\n",mul_time/(fp*n));
    printf("mpn mod.            : %.6f[ms]\n",mod_time/(fp*n));
    printf("mpn add + mul * 2.  : %.6f[ms]\n",(add_time+(mul_time*2))/(fp*n));
}

int test_mod(int mod){
    int i,j,n=0;
    float add_time=0,dbl_time=0,mul_time=0,sqr_time=0,mod_time=0,mod_mont_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("Field Mod test\n");
    mp_limb_t A[FPLIMB],B[FPLIMB],C[FPLIMB],C1[FPLIMB2];
    mp_limb_t At[FPLIMB2],Bt[FPLIMB2],Ct[FPLIMB2],test_mul[FPLIMB2];
    Fp Af,Bf,Cf;
    Fp test1,test2;

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
mod_time=0;
mod_mont_time=0;
for(i=0;i<mod;i++){
    mpn_random(A,FPLIMB);
    mpn_random(B,FPLIMB);

    mpn_mod(Af.x0,A,FPLIMB);
    mpn_mod(Bf.x0,B,FPLIMB);

    gettimeofday(&tv_A,NULL);
    mpn_mul_n(Ct,A,B,FPLIMB);
    Fp_mod(&test1,Ct,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    mod_time+=timedifference_msec(tv_A,tv_B);
    
    Fp_to_montgomery(&Af,&Af);
    Fp_to_montgomery(&Bf,&Bf);
    
    gettimeofday(&tv_A,NULL);
    Fp_mulmod_montgomery(&Cf,&Af,&Bf);
    gettimeofday(&tv_B,NULL);
    mod_mont_time+=timedifference_msec(tv_A,tv_B);

    Fp_mod_montgomery(&test2,&Cf);
    
    if(Fp_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    Fp_printf("",&test1);
	    Fp_printf("\n",&test2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("Fp mod.            : %.6f[ms]\n",mod_time/mod);
    printf("Fp mod montgomery. : %.6f[ms]\n",mod_mont_time/mod);

    return 0;
    
}

int test_EFp(int ecd,int eca,int scm){
    int i,j,n=0;
    float ecd_time=0,ecd_lazy_time=0,ecd_Projective_lazy_time=0,ecd_Jacobian_time=0,ecd_Jacobian_lazy_time=0,ecd_Jacobian_lazy_montgomery_time=0;
    float eca_time=0,eca_lazy_time=0,eca_Projective_lazy_time=0,eca_Jacobian_time=0,eca_Jacobian_lazy_time=0,eca_Jacobian_lazy_montgomery_time=0,eca_Mixture_time=0,eca_Mixture_lazy_time=0,eca_Mixture_lazy_montgomery_time=0;
    float scm_time=0,scm_lazy_time=0,scm_Jacobian_time=0,scm_Jacobian_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp A_EFp,B_EFp,test0,test1,test2,test3,test4,test5,test6,test7,test8;
    EFp A_EFpm,B_EFpm;
    EFpP A_EFpP,B_EFpP,testP1,testP2,testP3;
    EFpJ A_EFpJ,B_EFpJ,testJ1,testJ2,testJ3;
    EFpJ A_EFpJm,B_EFpJm;
    EFpJ A_EFpM,B_EFpM,testM1,testM2,testM3;
    EFpJ A_EFpMm,B_EFpMm,testMm;
    EFp_init(&A_EFp);
    EFp_init(&B_EFp);
    EFp_init(&test0);
    EFp_init(&test1);
    EFp_init(&test2);
    EFp_init(&test3);
    EFp_init(&test4);
    EFp_init(&test5);
    EFp_init(&test6);
    EFp_init(&test7);
    EFp_init(&test8);
    
    EFpP_init(&A_EFpP);
    EFpP_init(&B_EFpP);
    EFpP_init(&testP1);
    EFpP_init(&testP2);
    EFpP_init(&testP3);
    
    EFpJ_init(&A_EFpJ);
    EFpJ_init(&B_EFpJ);
    EFpJ_init(&B_EFpJm);
    EFpJ_init(&testJ1);
    EFpJ_init(&testJ2);
    EFpJ_init(&testJ3);
    EFpJ_init(&A_EFpM);
    EFpJ_init(&B_EFpM);
    EFpJ_init(&testM1);
    EFpJ_init(&testM2);
    EFpJ_init(&testM3);
    
    EFp_init(&A_EFpm);
    EFp_init(&B_EFpm);
    EFpJ_init(&A_EFpMm);
    EFpJ_init(&B_EFpMm);
    EFpJ_init(&A_EFpJm);
    EFpJ_init(&B_EFpJm);
    EFpJ_init(&testMm);
    
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    gmp_randseed_ui(state,1);


printf("------------------------------------------------------------------------------------\n");
    printf("EFp_ECD test\n");
n=100;
for(i=0;i<ecd;i++){
    EFp_rational_point(&B_EFp);
    EFp_to_EFpP(&B_EFpP,&B_EFp);
    EFp_to_EFpJ(&B_EFpJ,&B_EFp);
	EFp_to_montgomery(&B_EFpm,&B_EFp);
    EFp_to_EFpJ_montgomery(&B_EFpJm,&B_EFpm);
	//Fp_to_montgomery(B_EFpJm.z.x0,B_EFpJm.z.x0);
    EFp_ECD(&B_EFp,&B_EFp);
    EFp_ECD_Projective_lazy(&B_EFpP,&B_EFpP);
    EFp_ECD_Jacobian(&B_EFpJ,&B_EFpJ);
    EFp_ECD_Jacobian_lazy_montgomery(&B_EFpJm,&B_EFpJm);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECD(&test0,&B_EFp);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECD_lazy(&test1,&B_EFp);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECD_Projective_lazy(&testP1,&B_EFpP);
    gettimeofday(&tv_B,NULL);
    ecd_Projective_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Projective(&test2,&testP1);
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECD_Jacobian(&testJ1,&B_EFpJ);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Jacobian(&test3,&testJ1);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECD_Jacobian_lazy(&testJ2,&B_EFpJ);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Jacobian(&test4,&testJ2);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECD_Jacobian_lazy_montgomery(&testJ3,&B_EFpJm);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Jacobian_montgomery(&test5,&testJ3);
    EFp_mod_montgomery(&test5,&test5);
    
    if(EFp_cmp(&test0,&test1)!=0 || EFp_cmp(&test1,&test2)!=0 || EFp_cmp(&test1,&test3)!=0 || EFp_cmp(&test1,&test4)!=0 || EFp_cmp(&test1,&test5)!=0){
        printf("failed!\n\n");
	    EFp_println("test1=",&test1);
	    EFp_println("test2=",&test2);
	    EFp_println("test3=",&test3);
	    EFp_println("test4=",&test4);
	    EFp_println("test5=",&test5);
	    EFpP_printf("testP1=",&testP1);
	    EFpJ_printf("\ntestJ1=",&testJ1);
	    EFpJ_printf("\ntestJ2=",&testJ2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp ECD.                            : %.6f[ms]\n",ecd_time/ecd);
    printf("EFp ECD lazy.                       : %.6f[ms]\n",ecd_lazy_time/ecd);
    printf("EFp ECD Projective lazy.            : %.6f[ms]\n",ecd_Projective_lazy_time/ecd);
    printf("EFp ECD Jacobian.                   : %.6f[ms]\n",ecd_Jacobian_time/ecd);
    printf("EFp ECD Jacobian lazy.              : %.6f[ms]\n",ecd_Jacobian_lazy_time/ecd);
    printf("EFp ECD Jacobian lazy montgomery.   : %.6f[ms]\n",ecd_Jacobian_lazy_montgomery_time/ecd);


printf("------------------------------------------------------------------------------------\n");
    printf("EFp_ECA test\n");
    gmp_randseed_ui(state,1);
    n=100;
    EFp_rational_point(&A_EFp);
	EFp_to_montgomery(&A_EFpm,&A_EFp);
    EFp_to_EFpP(&A_EFpP,&A_EFp);
    EFp_to_EFpJ(&A_EFpJ,&A_EFp);
    EFp_to_EFpJ_montgomery(&A_EFpJm,&A_EFpm);
    EFp_to_EFpJ(&A_EFpM,&A_EFp);
    EFp_to_EFpJ_montgomery(&A_EFpMm,&A_EFpm);
	//Fp_to_montgomery(A_EFpMm.z.x0,A_EFpMm.z.x0);
	
sleep(1);
for(i=0;i<eca;i++){
    EFp_rational_point(&B_EFp);
	EFp_to_montgomery(&B_EFpm,&B_EFp);
    EFp_to_EFpP(&B_EFpP,&B_EFp);
    EFp_to_EFpJ(&B_EFpJ,&B_EFp);
    EFp_to_EFpJ_montgomery(&B_EFpJm,&B_EFpm);
	//Fp_to_montgomery(B_EFpJm.z.x0,B_EFpJm.z.x0);
    
    EFp_ECD(&B_EFp,&B_EFp);
    EFp_ECD_Projective_lazy(&B_EFpP,&B_EFpP);
    EFp_ECD_Jacobian(&B_EFpJ,&B_EFpJ);
    EFp_ECD_Jacobian_lazy_montgomery(&B_EFpJm,&B_EFpJm);
    
    
    //Fp_println("A_EFpJm=",&A_EFpJm.x);
    //Fp_println("A_EFpJm=",&A_EFpJm.y);
    //Fp_println("A_EFpJm=",&A_EFpJm.z);
    
    //Fp_println("A_EFpMm=",&A_EFpMm.x);
    //Fp_println("A_EFpMm=",&A_EFpMm.y);
    //Fp_println("A_EFpMm=",&A_EFpMm.z);
    
    //Fp_println("B_EFpJm=",&B_EFpJm.x);
    //Fp_println("B_EFpJm=",&B_EFpJm.y);
    //Fp_println("B_EFpJm=",&B_EFpJm.z);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA(&test0,&A_EFp,&B_EFp);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B)/n;
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_lazy(&test1,&A_EFp,&B_EFp);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_Projective_lazy(&testP1,&A_EFpP,&B_EFpP);
    gettimeofday(&tv_B,NULL);
    eca_Projective_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Projective(&test2,&testP1);
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_Jacobian(&testJ1,&A_EFpJ,&B_EFpJ);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Jacobian(&test3,&testJ1);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_Jacobian_lazy(&testJ2,&A_EFpJ,&B_EFpJ);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Jacobian(&test4,&testJ2);
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_Jacobian_lazy_montgomery(&testJ3,&B_EFpJm,&A_EFpJm);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B)/n;
    //Fp_println("testJ3.x=",&testJ3.x);
    //Fp_println("testJ3.y=",&testJ3.y);
    //Fp_println("testJ3.z=",&testJ3.z);
    EFp_Jacobian_montgomery(&test5,&testJ3);
    EFp_mod_montgomery(&test5,&test5);
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_Mixture(&testM1,&B_EFpJ,&A_EFpM);
    gettimeofday(&tv_B,NULL);
    eca_Mixture_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Jacobian(&test6,&testM1);
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_Mixture_lazy(&testM2,&B_EFpJ,&A_EFpM);
    gettimeofday(&tv_B,NULL);
    eca_Mixture_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    EFp_Jacobian(&test7,&testM2);
    
    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)EFp_ECA_Mixture_lazy_montgomery(&testMm,&B_EFpJm,&A_EFpMm);
    gettimeofday(&tv_B,NULL);
    eca_Mixture_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B)/n;
    //Fp_println("testMm.x=",&testMm.x);
    //Fp_println("testMm.y=",&testMm.y);
    //Fp_println("testMm.z=",&testMm.z);
    EFp_Jacobian_montgomery(&test8,&testMm);
    EFp_mod_montgomery(&test8,&test8);
    
    if(EFp_cmp(&test0,&test1)!=0 || EFp_cmp(&test1,&test2)!=0 || EFp_cmp(&test1,&test3)!=0 || EFp_cmp(&test1,&test4)!=0 || EFp_cmp(&test1,&test5)!=0 || EFp_cmp(&test1,&test6)!=0 || EFp_cmp(&test1,&test7)!=0 || EFp_cmp(&test1,&test8)!=0){
        printf("failed!\n\n");
        EFp_println("test0=",&test0);
	    if(EFp_cmp(&test0,&test1)!=0)EFp_println("test1=",&test1);
	    if(EFp_cmp(&test1,&test2)!=0)EFp_println("test2=",&test2);
	    if(EFp_cmp(&test1,&test3)!=0)EFp_println("test3=",&test3);
	    if(EFp_cmp(&test1,&test4)!=0)EFp_println("test4=",&test4);
	    if(EFp_cmp(&test1,&test5)!=0)EFp_println("test5=",&test5);
	    if(EFp_cmp(&test1,&test6)!=0)EFp_println("test6=",&test6);
	    if(EFp_cmp(&test1,&test7)!=0)EFp_println("test7=",&test7);
	    if(EFp_cmp(&test1,&test8)!=0)EFp_println("test8=",&test8);
	    //EFpP_printf("testP1=",&testP1);
	    //EFpJ_printf("\ntestJ1=",&testJ1);
	    //EFpJ_printf("\ntestJ2=",&testJ2);
	    //EFpJ_printf("\ntestM1=",&testM1);
	    //EFpJ_printf("\ntestMm=",&testMm);
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp ECA.                           : %.6f[ms]\n",eca_time/eca);
    printf("EFp ECA lazy.                      : %.6f[ms]\n",eca_lazy_time/eca);
    printf("EFp ECA Projective lazy.           : %.6f[ms]\n",eca_Projective_lazy_time/eca);
    printf("EFp ECA Jacobian.                  : %.6f[ms]\n",eca_Jacobian_time/eca);
    printf("EFp ECA Jacobian lazy.             : %.6f[ms]\n",eca_Jacobian_lazy_time/eca);
    printf("EFp ECA Jacobian lazy montgomery.  : %.6f[ms]\n",eca_Jacobian_lazy_montgomery_time/eca);
    printf("EFp ECA Mixture.                   : %.6f[ms]\n",eca_Mixture_time/eca);
    printf("EFp ECA Mixture lazy.              : %.6f[ms]\n",eca_Mixture_lazy_time/eca);
    printf("EFp ECA Mixture lazy montgomery.   : %.6f[ms]\n",eca_Mixture_lazy_montgomery_time/eca);

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
	    EFp_printf("",&test1);
	    EFp_printf("\n",&test2);
	    EFp_printf("\n",&test3);
	    printf("\n\n");
	    return 1;
    }
}

    printf("EFp SCM.               : %.6f[ms]\n",scm_time/scm);
    printf("EFp SCM lazy.          : %.6f[ms]\n",scm_lazy_time/scm);
    printf("EFp SCM Jacobian.      : %.6f[ms]\n",scm_Jacobian_time/scm);
    printf("EFp SCM Jacobian Lazy. : %.6f[ms]\n",scm_Jacobian_lazy_time/scm);

    return 0;
}
int test_EFp2(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0,ecd_Jacobian_time=0,ecd_Jacobian_lazy_time=0,ecd_Jacobian_lazy_montgomery_time=0;
    float eca_time=0,eca_lazy_time=0,eca_Jacobian_time=0,eca_Jacobian_lazy_time=0,eca_Jacobian_lazy_montgomery_time=0,eca_Mixture_lazy_montgomery_time=0;
    float scm_time=0,scm_lazy_time=0,scm_Jacobian_time=0,scm_Jacobian_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    EFp2 A_EFp2,B_EFp2,test0,test1,test2,test3,test4,test5;
    EFp2 A_EFp2m,B_EFp2m;
    EFpJ2 A_EFpJ2,B_EFpJ2,testZ1,testZ2,testZ3,testZ4,testZ5;
    EFpJ2 A_EFpJ2m,B_EFpJ2m;
    EFp2_init(&A_EFp2);
    EFp2_init(&B_EFp2);
    EFp2_init(&A_EFp2m);
    EFp2_init(&B_EFp2m);
    EFp2_init(&test0);
    EFp2_init(&test1);
    EFp2_init(&test2);
    EFp2_init(&test3);
    EFp2_init(&test4);
    EFpJ2_init(&A_EFpJ2);
    EFpJ2_init(&B_EFpJ2);
    EFpJ2_init(&testZ1);
    EFpJ2_init(&testZ2);
    EFpJ2_init(&testZ3);
    EFpJ2_init(&testZ4);
    EFpJ2_init(&testZ5);
    EFpJ2_init(&A_EFpJ2m);
    EFpJ2_init(&B_EFpJ2m);
    
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    gmp_randseed_ui(state,1);


printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_ECD test\n");
    
for(i=0;i<ecd;i++){

    EFp2_rational_point(&B_EFp2);
	EFp2_to_montgomery(&B_EFp2m,&B_EFp2);
	
    EFp2_to_EFpJ2(&B_EFpJ2,&B_EFp2);
    EFp2_to_EFpJ2_montgomery(&B_EFpJ2m,&B_EFp2m);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD(&test0,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD_lazy(&test1,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECD_Jacobian(&testZ2,&B_EFpJ2);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_time+=timedifference_msec(tv_A,tv_B);
    EFp2_Jacobian(&test2,&testZ2);

    gettimeofday(&tv_A,NULL);
    EFp2_ECD_Jacobian_lazy(&testZ3,&B_EFpJ2);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);
    EFp2_Jacobian(&test3,&testZ3);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECD_Jacobian_lazy_montgomery(&testZ4,&B_EFpJ2m);
    gettimeofday(&tv_B,NULL);
    ecd_Jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    EFp2_Jacobian_montgomery(&test4,&testZ4);
    EFp2_mod_montgomery(&test4,&test4);
        
    if(EFp2_cmp(&test0,&test1)!=0 || EFp2_cmp(&test1,&test2)!=0 || EFp2_cmp(&test1,&test3)!=0 || EFp2_cmp(&test1,&test4)!=0){
        printf("failed!\n\n");
	    EFp2_printf("",&test1);
	    EFp2_printf("\n",&test2);
	    EFp2_printf("\n",&test3);
	    EFp2_printf("\n",&test4);
	    EFpJ2_printf("\ntestZ2=",&testZ2);
	    EFpJ2_printf("\ntestZ3=",&testZ3);
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp2 ECD.                          : %.4f[ms]\n",ecd_time/ecd);
    printf("EFp2 ECD lazy.                     : %.4f[ms]\n",ecd_lazy_time/ecd);
    printf("EFp2 ECD Jacobian.                 : %.4f[ms]\n",ecd_Jacobian_time/ecd);
    printf("EFp2 ECD Jacobian lazy.            : %.4f[ms]\n",ecd_Jacobian_lazy_time/ecd);
    printf("EFp2 ECD Jacobian lazy montgomery. : %.4f[ms]\n",ecd_Jacobian_lazy_montgomery_time/ecd);


printf("------------------------------------------------------------------------------------\n");
    printf("EFp2_ECA test\n");
    
    EFp2_rational_point(&A_EFp2);
	EFp2_to_montgomery(&A_EFp2m,&A_EFp2);

    EFp2_to_EFpJ2(&A_EFpJ2,&A_EFp2);
    EFp2_to_EFpJ2_montgomery(&A_EFpJ2m,&A_EFp2m);
for(i=0;i<eca;i++){

    EFp2_rational_point(&B_EFp2);
	EFp2_to_montgomery(&B_EFp2m,&B_EFp2);
	
    EFp2_to_EFpJ2(&B_EFpJ2,&B_EFp2);
    EFp2_to_EFpJ2_montgomery(&B_EFpJ2m,&B_EFp2m);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA(&test0,&A_EFp2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECA_lazy(&test1,&A_EFp2,&B_EFp2);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Jacobian(&testZ2,&A_EFpJ2,&B_EFpJ2);
    EFp2_Jacobian(&test2,&testZ2);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Jacobian_lazy(&testZ3,&A_EFpJ2,&B_EFpJ2);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_lazy_time+=timedifference_msec(tv_A,tv_B);
    EFp2_Jacobian(&test3,&testZ3);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Jacobian_lazy_montgomery(&testZ4,&B_EFpJ2m,&A_EFpJ2m);
    gettimeofday(&tv_B,NULL);
    eca_Jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    EFp2_Jacobian_montgomery(&test4,&testZ4);
    EFp2_mod_montgomery(&test4,&test4);
    
    gettimeofday(&tv_A,NULL);
    EFp2_ECA_Mixture_lazy_montgomery(&testZ5,&B_EFpJ2m,&A_EFpJ2m);
    gettimeofday(&tv_B,NULL);
    eca_Mixture_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    EFp2_Jacobian_montgomery(&test5,&testZ5);
    EFp2_mod_montgomery(&test5,&test5);
    
    if(EFp2_cmp(&test0,&test1)!=0 || EFp2_cmp(&test1,&test2)!=0 || EFp2_cmp(&test1,&test3)!=0 || EFp2_cmp(&test1,&test4)!=0 || EFp2_cmp(&test1,&test5)!=0){
        printf("failed!\n\n");
	    EFp2_printf("test1=",&test1);
	    EFp2_printf("\ntest2=",&test2);
	    EFp2_printf("\ntest3=",&test3);
	    EFp2_printf("\ntest4=",&test4);
	    EFp2_printf("\ntest5=",&test5);
	    //EFpJ2_printf("\ntestZ2=",&testZ2);
	    //EFpJ2_printf("\ntestZ3=",&testZ3);
	    printf("\n\n");
	    return 1;
    }
}
    printf("EFp2 ECA.                          : %.4f[ms]\n",eca_time/eca);
    printf("EFp2 ECA lazy.                     : %.4f[ms]\n",eca_lazy_time/eca);
    printf("EFp2 ECA Jacobian.                 : %.4f[ms]\n",eca_Jacobian_time/eca);
    printf("EFp2 ECA Jacobian lazy.            : %.4f[ms]\n",eca_Jacobian_lazy_time/eca);
    printf("EFp2 ECA Jacobian lazy montgomery. : %.4f[ms]\n",eca_Jacobian_lazy_montgomery_time/eca);
    printf("EFp2 ECA Mixture lazy montgomery.  : %.4f[ms]\n",eca_Mixture_lazy_montgomery_time/eca);

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
	    EFp2_printf("",&test1);
	    EFp2_printf("\n",&test2);
	    EFp2_printf("\n",&test3);
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
	    EFp12_printf("",&test1);
	    EFp12_printf("\n",&test2);
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
	    EFp12_printf("",&test1);
	    EFp12_printf("\n",&test2);
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
	    EFp12_printf("",&test1);
	    EFp12_printf("\n",&test2);
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
    Fp12_printf("",&A_Fp12);
    printf("\n\n");
    
    printf("frobenius\n");
    Fp12_frobenius_map_p10(&test1,&A_Fp12);
    Fp12_printf("",&test1);
    printf("\n\n");
    
    mpz_pow_ui(exp,prime_z,10);
    Fp12_pow(&test2,&A_Fp12,exp);
    Fp12_printf("",&test2);
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
    EFp12_printf("",&test1);
    printf("\n\n");
    
    EFp2_skew_frobenius_map_p10(&twisted_Q,&twisted_Q);
    EFp2_to_EFp12(&test2,&twisted_Q);
    EFp12_printf("",&test2);
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
    EFp12_printf("Q\n",&Q);
    printf("\n\n");
    
    EFp12_to_EFp2(&twist_Q,&Q);
    EFp2_ECD(&twist_Q,&twist_Q);
    EFp2_to_EFp12(&test1,&twist_Q);
    EFp12_printf("",&test1);
    printf("\n\n");
    
    EFp12_ECD(&test2,&Q);
    EFp12_printf("",&test2);
    printf("\n\n");
    
    if(Fp12_cmp(&test1.x,&test2.x)==0 && Fp12_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
}
void test_All(){
	int point,Field,EFp,EFp2,EFp12,G1,G2,G3,pairing;
	
printf("====================================================================================\n");
    printf("Test time\n");
	
	point = BLS12_test_rational_point();
	Field = test_Field(100,100,100,10,10);
	EFp = test_EFp(100,100,10);
	EFp2 = test_EFp2(10,10,10);
	EFp12 = test_EFp12(10,10,10);
	G1 = BLS12_test_G1_SCM(10);
	G2 = BLS12_test_G2_SCM(10);
	G3 = BLS12_test_G3_EXP(10);
	pairing = BLS12_test_opt_ate_pairing(10);
	
printf("====================================================================================\n");
    printf("Test Result\n\n");

	printf("Test Point      :");
	if(point==0) printf("Success\n");
	else printf("Failed\n");
	
	printf("Test Field      :");
	if(Field==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test EFp        :");
	if(EFp==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test EFp2       :");
	if(EFp2==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test EFp12      :");
	if(EFp12==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G1         :");
	if(G1==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G2         :");
	if(G2==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test G3         :");
	if(G3==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test pairing    :");
	if(pairing==0) printf("Success\n");
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