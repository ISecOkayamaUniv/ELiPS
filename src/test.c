#include <ELiPS/test.h>
/*----------------------------------------------------------------------------*/
//test
/*
void test_fp(int fp_n){
    int i,j,n=10000;
    float add_time=0,dbl_time=0,mul_time=0,sqr_time=0,mod_time=0,mod_mont_time=0;
    //struct timespec tv_a,tv_b;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("Field fp test\n");
    mp_limb_t A[FPLIMB],B[FPLIMB],C[FPLIMB],C1[FPLIMB2];
    mp_limb_t At[FPLIMB2],Bt[FPLIMB2],Ct[FPLIMB2],test_mul[FPLIMB2];
    fp_t test1,test2;

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

for(i=0;i<fp_n;i++){
    mpn_random(A,FPLIMB);
    mpn_random(B,FPLIMB);

    mpn_mod(A,A,FPLIMB);
    mpn_mod(B,B,FPLIMB);

    mpn_mul_n(At,A,A,FPLIMB);
    mpn_mul_n(Bt,B,B,FPLIMB);

	//TODO:Raspberrypi error

    //clock_gettime(CLOCK_REALTIME, &tv_a);
    //for(j=0;j<n;j++)    mpn_add_n(test_mul,At,At,FPLIMB2);
    //clock_gettime(CLOCK_REALTIME, &tv_b);
    //add_time+=timedifference_nsec(tv_a,tv_b);

    //clock_gettime(CLOCK_REALTIME, &tv_a);
    //for(j=0;j<n;j++)    mpn_lshift(test_mul,At,FPLIMB2,1);
    //clock_gettime(CLOCK_REALTIME, &tv_b);
    //dbl_time+=timedifference_nsec(tv_a,tv_b);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    mpn_add_n(test_mul,At,At,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    add_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    mpn_lshift(test_mul,At,FPLIMB2,1);
    gettimeofday(&tv_B,NULL);
    dbl_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    mpn_sqr(test_mul,A,FPLIMB);
    gettimeofday(&tv_B,NULL);
    sqr_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    mpn_mul_n(test_mul,A,A,FPLIMB);
    gettimeofday(&tv_B,NULL);
    mul_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)    fp_mod(&test1,test_mul,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    mod_time+=timedifference_msec(tv_A,tv_B);
    }
    add_time=add_time/(1000*1000);
    dbl_time=dbl_time/(1000*1000);
    printf("mpn add.            : %.6f[ms]\n",add_time/(fp_n*n));
    printf("mpn dbl.            : %.6f[ms]\n",dbl_time/(fp_n*n));
    printf("mpn sqr.            : %.6f[ms]\n",sqr_time/(fp_n*n));
    printf("mpn mul.            : %.6f[ms]\n",mul_time/(fp_n*n));
    printf("mpn mod.            : %.6f[ms]\n",mod_time/(fp_n*n));
    printf("mpn add + mul * 2.  : %.6f[ms]\n",(add_time+(mul_time*2))/(fp_n*n));
}
*/
int test_mod(int mod){
    int i,j,n=0;
    float add_time=0,dbl_time=0,mul_time=0,sqr_time=0,mod_time=0,mod_mont_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    printf("Field Mod test\n");
    mp_limb_t A[FPLIMB],B[FPLIMB],C[FPLIMB],C1[FPLIMB2];
    mp_limb_t At[FPLIMB2],Bt[FPLIMB2],Ct[FPLIMB2],test_mul[FPLIMB2];
    fp_t Af,Bf,Cf;
    fp_t test1,test2;

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));

mod_time=0;
mod_mont_time=0;
for(i=0;i<mod;i++){
    mpn_random(A,FPLIMB);
    mpn_random(B,FPLIMB);

    mpn_mod(Af.x0,A,FPLIMB);
    mpn_mod(Bf.x0,B,FPLIMB);


    mpn_mul_n(Ct,Af.x0,Bf.x0,FPLIMB);

    gettimeofday(&tv_A,NULL);
    fp_mod(&test1,Ct,FPLIMB2);
    gettimeofday(&tv_B,NULL);
    mod_time+=timedifference_msec(tv_A,tv_B);



    fp_to_montgomery(&Af,&Af);
    fp_to_montgomery(&Bf,&Bf);

    mpn_mul_n(Ct,Af.x0,Bf.x0,FPLIMB);

    gettimeofday(&tv_A,NULL);
    mpn_mod_montgomery(Cf.x0,FPLIMB,Ct,FPLIMB2);
//    fp_mulmod_montgomery(&Cf,&Af,&Bf);
    gettimeofday(&tv_B,NULL);
    mod_mont_time+=timedifference_msec(tv_A,tv_B);

    fp_mod_montgomery(&test2,&Cf);

    if(fp_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    fp_printf("",&test1);
	    fp_printf("\n",&test2);
	    printf("\n\n");
	    return 1;
    }
}
    printf("fp mod.            : %.6f[ms]\n",mod_time/mod);
    printf("fp mod montgomery. : %.6f[ms]\n",mod_mont_time/mod);

    return 0;

}
int test_efp(int ecd,int eca){
    int i,j,n=0;
    float ecd_time=0,ecd_lazy_time=0,ecd_projective_lazy_time=0,ecd_jacobian_time=0,ecd_jacobian_lazy_time=0,ecd_jacobian_lazy_montgomery_time=0;
    float eca_time=0,eca_lazy_time=0,eca_projective_lazy_time=0,eca_jacobian_time=0,eca_jacobian_lazy_time=0,eca_jacobian_lazy_montgomery_time=0,eca_mixture_time=0,eca_mixture_lazy_time=0,eca_mixture_lazy_montgomery_time=0;
    float scm_time=0,scm_lazy_time=0,scm_jacobian_time=0,scm_jacobian_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    efp_t A_efp,B_efp,test0,test1,test2,test3,test4,test5,test6,test7,test8;
    efp_t A_efpm,B_efpm;
    efp_projective_t A_efpP,B_efpP,testP1,testP2,testP3;
    efp_jacobian_t A_efpJ,B_efpJ,testJ1,testJ2,testJ3;
    efp_jacobian_t A_efpJm,B_efpJm;
    efp_jacobian_t A_efpM,B_efpM,testM1,testM2,testM3;
    efp_jacobian_t A_efpMm,B_efpMm,testMm;
    efp_init(&A_efp);
    efp_init(&B_efp);
    efp_init(&test0);
    efp_init(&test1);
    efp_init(&test2);
    efp_init(&test3);
    efp_init(&test4);
    efp_init(&test5);
    efp_init(&test6);
    efp_init(&test7);
    efp_init(&test8);

    efp_projective_init(&A_efpP);
    efp_projective_init(&B_efpP);
    efp_projective_init(&testP1);
    efp_projective_init(&testP2);
    efp_projective_init(&testP3);

    efp_jacobian_init(&A_efpJ);
    efp_jacobian_init(&B_efpJ);
    efp_jacobian_init(&B_efpJm);
    efp_jacobian_init(&testJ1);
    efp_jacobian_init(&testJ2);
    efp_jacobian_init(&testJ3);
    efp_jacobian_init(&A_efpM);
    efp_jacobian_init(&B_efpM);
    efp_jacobian_init(&testM1);
    efp_jacobian_init(&testM2);
    efp_jacobian_init(&testM3);

    efp_init(&A_efpm);
    efp_init(&B_efpm);
    efp_jacobian_init(&A_efpMm);
    efp_jacobian_init(&B_efpMm);
    efp_jacobian_init(&A_efpJm);
    efp_jacobian_init(&B_efpJm);
    efp_jacobian_init(&testMm);

    mpz_t scalar;
    mpz_init(scalar);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,1);


printf("------------------------------------------------------------------------------------\n");
    printf("efp_ecd test\n");
n=100;
for(i=0;i<ecd;i++){
    efp_set_random(&B_efp);
    efp_affine_to_projective(&B_efpP,&B_efp);
    efp_affine_to_jacobian(&B_efpJ,&B_efp);
	efp_to_montgomery(&B_efpm,&B_efp);
    efp_affine_to_jacobian_montgomery(&B_efpJm,&B_efpm);
	//fp_to_montgomery(B_efpJm.z.x0,B_efpJm.z.x0);
    efp_ecd(&B_efp,&B_efp);
    //efp_ecd_projective_lazy(&B_efpP,&B_efpP);
    efp_ecd_jacobian_lazy_montgomery(&B_efpJm,&B_efpJm);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)efp_ecd(&test1,&B_efp);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    //for(j=0;j<n;j++)efp_ecd_projective_lazy(&testP1,&B_efpP);
    gettimeofday(&tv_B,NULL);
    ecd_projective_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    efp_projective_to_affine(&test2,&testP1);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)efp_ecd_jacobian_lazy_montgomery(&testJ3,&B_efpJm);
    gettimeofday(&tv_B,NULL);
    ecd_jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B)/n;
    efp_jacobian_to_affine_montgomery(&test5,&testJ3);
    efp_mod_montgomery(&test5,&test5);

    if(efp_cmp(&test1,&test2)!=0 || efp_cmp(&test1,&test5)!=0){
        printf("failed!\n\n");
	    efp_println("test1=",&test1);
	    efp_println("test2=",&test2);
	    efp_println("test5=",&test5);
	    printf("\n\n");
	    return 1;
    }
}
    printf("efp ecd.                            : %.6f[ms]\n",ecd_time/ecd);
    printf("efp ecd projective lazy.            : %.6f[ms]\n",ecd_projective_lazy_time/ecd);
    printf("efp ecd jacobian lazy montgomery.   : %.6f[ms]\n",ecd_jacobian_lazy_montgomery_time/ecd);


printf("------------------------------------------------------------------------------------\n");
    printf("efp_eca test\n");
    gmp_randseed_ui(state,1);
    n=100;
    efp_set_random(&A_efp);
	efp_to_montgomery(&A_efpm,&A_efp);
    efp_affine_to_projective(&A_efpP,&A_efp);
    efp_affine_to_jacobian(&A_efpJ,&A_efp);
    efp_affine_to_jacobian_montgomery(&A_efpJm,&A_efpm);
    efp_affine_to_jacobian(&A_efpM,&A_efp);
    efp_affine_to_jacobian_montgomery(&A_efpMm,&A_efpm);
	//fp_to_montgomery(A_efpMm.z.x0,A_efpMm.z.x0);

sleep(1);
for(i=0;i<eca;i++){
    efp_set_random(&B_efp);
	efp_to_montgomery(&B_efpm,&B_efp);
    efp_affine_to_projective(&B_efpP,&B_efp);
    efp_affine_to_jacobian(&B_efpJ,&B_efp);
    efp_affine_to_jacobian_montgomery(&B_efpJm,&B_efpm);
	//fp_to_montgomery(B_efpJm.z.x0,B_efpJm.z.x0);

    efp_ecd(&B_efp,&B_efp);
    //efp_ecd_projective_lazy(&B_efpP,&B_efpP);\
    efp_ecd_jacobian_lazy_montgomery(&B_efpJm,&B_efpJm);


    //fp_println("A_efpJm=",&A_efpJm.x);
    //fp_println("A_efpJm=",&A_efpJm.y);
    //fp_println("A_efpJm=",&A_efpJm.z);

    //fp_println("A_efpMm=",&A_efpMm.x);
    //fp_println("A_efpMm=",&A_efpMm.y);
    //fp_println("A_efpMm=",&A_efpMm.z);

    //fp_println("B_efpJm=",&B_efpJm.x);
    //fp_println("B_efpJm=",&B_efpJm.y);
    //fp_println("B_efpJm=",&B_efpJm.z);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)efp_eca(&test1,&A_efp,&B_efp);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B)/n;

    gettimeofday(&tv_A,NULL);
    //for(j=0;j<n;j++)efp_eca_projective_lazy(&testP1,&A_efpP,&B_efpP);
    gettimeofday(&tv_B,NULL);
    eca_projective_lazy_time+=timedifference_msec(tv_A,tv_B)/n;
    efp_projective_to_affine(&test2,&testP1);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)efp_eca_jacobian_lazy_montgomery(&testJ3,&B_efpJm,&A_efpJm);
    gettimeofday(&tv_B,NULL);
    eca_jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B)/n;
    //fp_println("testJ3.x=",&testJ3.x);
    //fp_println("testJ3.y=",&testJ3.y);
    //fp_println("testJ3.z=",&testJ3.z);
    efp_jacobian_to_affine_montgomery(&test5,&testJ3);
    efp_mod_montgomery(&test5,&test5);

    gettimeofday(&tv_A,NULL);
    for(j=0;j<n;j++)efp_eca_mixture_lazy_montgomery(&testMm,&B_efpJm,&A_efpMm);
    gettimeofday(&tv_B,NULL);
    eca_mixture_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B)/n;
    //fp_println("testMm.x=",&testMm.x);
    //fp_println("testMm.y=",&testMm.y);
    //fp_println("testMm.z=",&testMm.z);
    efp_jacobian_to_affine_montgomery(&test8,&testMm);
    efp_mod_montgomery(&test8,&test8);

    if(efp_cmp(&test1,&test2)!=0 || efp_cmp(&test1,&test5)!=0 || efp_cmp(&test1,&test8)!=0){
        printf("failed!\n\n");
        efp_println("test0=",&test0);
	    if(efp_cmp(&test0,&test1)!=0)efp_println("test1=",&test1);
	    if(efp_cmp(&test1,&test2)!=0)efp_println("test2=",&test2);
	    if(efp_cmp(&test1,&test5)!=0)efp_println("test5=",&test5);
	    if(efp_cmp(&test1,&test8)!=0)efp_println("test8=",&test8);
	    //efp_projective_printf("testP1=",&testP1);
	    //efp_jacobian_printf("\ntestJ1=",&testJ1);
	    //efp_jacobian_printf("\ntestJ2=",&testJ2);
	    //efp_jacobian_printf("\ntestM1=",&testM1);
	    //efp_jacobian_printf("\ntestMm=",&testMm);
	    printf("\n\n");
	    return 1;
    }
}
    printf("efp eca.                           : %.6f[ms]\n",eca_time/eca);
    printf("efp eca projective lazy.           : %.6f[ms]\n",eca_projective_lazy_time/eca);
    printf("efp eca jacobian lazy montgomery.  : %.6f[ms]\n",eca_jacobian_lazy_montgomery_time/eca);
    printf("efp eca mixture lazy montgomery.   : %.6f[ms]\n",eca_mixture_lazy_montgomery_time/eca);


    return 0;
}
int test_efp2(int ecd,int eca){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0,ecd_projective_lazy_time=0,ecd_jacobian_time=0,ecd_jacobian_lazy_time=0,ecd_jacobian_lazy_montgomery_time=0;
    float eca_time=0,eca_lazy_time=0,eca_projective_lazy_time=0,eca_jacobian_time=0,eca_jacobian_lazy_time=0,eca_jacobian_lazy_montgomery_time=0,eca_mixture_lazy_montgomery_time=0;
    float scm_time=0,scm_lazy_time=0,scm_jacobian_time=0,scm_jacobian_lazy_time=0;
    cost tmp,eca_cost,eca_projective_lazy_cost;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    efp2_t A_efp2,B_efp2,test0,test1,test2,test3,test4,test5,test6;
    efp2_t A_efp2m,B_efp2m;
    efp2_projective_t A_efp2_projective,B_efp2_projective,testP1,testP2;
    efp2_jacobian_t A_efp2_jacobian,B_efp2_jacobian,testZ1,testZ2,testZ3,testZ4,testZ5,testZ6;
    efp2_jacobian_t A_efp2_jacobianm,B_efp2_jacobianm;
    efp2_init(&A_efp2);
    efp2_init(&B_efp2);
    efp2_init(&A_efp2m);
    efp2_init(&B_efp2m);
    efp2_init(&test0);
    efp2_init(&test1);
    efp2_init(&test2);
    efp2_init(&test3);
    efp2_init(&test4);
    efp2_projective_init(&A_efp2_projective);
    efp2_projective_init(&B_efp2_projective);
    efp2_projective_init(&testP1);
    efp2_projective_init(&testP2);
    efp2_jacobian_init(&A_efp2_jacobian);
    efp2_jacobian_init(&B_efp2_jacobian);
    efp2_jacobian_init(&testZ1);
    efp2_jacobian_init(&testZ2);
    efp2_jacobian_init(&testZ3);
    efp2_jacobian_init(&testZ4);
    efp2_jacobian_init(&testZ5);
    efp2_jacobian_init(&testZ6);
    efp2_jacobian_init(&A_efp2_jacobianm);
    efp2_jacobian_init(&B_efp2_jacobianm);

    mpz_t scalar;
    mpz_init(scalar);

    cost_init(&tmp);
    cost_init(&eca_cost);
    cost_init(&eca_projective_lazy_cost);

    gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    gmp_randseed_ui(state,1);


printf("------------------------------------------------------------------------------------\n");
    printf("efp2_ecd test\n");

for(i=0;i<ecd;i++){

    efp2_rational_point(&B_efp2);
	efp2_to_montgomery(&B_efp2m,&B_efp2);

    efp2_affine_to_projective(&B_efp2_projective,&B_efp2);
    efp2_affine_to_jacobian(&B_efp2_jacobian,&B_efp2);
    efp2_affine_to_jacobian_montgomery(&B_efp2_jacobianm,&B_efp2m);

    gettimeofday(&tv_A,NULL);
    efp2_ecd(&test1,&B_efp2);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    //efp2_ecd_projective_lazy(&testP2,&B_efp2_projective);
    gettimeofday(&tv_B,NULL);
    ecd_projective_lazy_time+=timedifference_msec(tv_A,tv_B);
    efp2_projective_to_affine(&test2,&testP2);

    gettimeofday(&tv_A,NULL);
    efp2_ecd_jacobian_lazy_montgomery(&testZ5,&B_efp2_jacobianm);
    gettimeofday(&tv_B,NULL);
    ecd_jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    efp2_jacobian_to_affine_montgomery(&test5,&testZ5);
    efp2_mod_montgomery(&test5,&test5);

    if(efp2_cmp(&test1,&test2)!=0 || efp2_cmp(&test1,&test5)!=0){
        printf("failed!\n\n");
	    efp2_printf("",&test1);
	    efp2_printf("\n",&test2);
	    efp2_printf("\n",&test5);
	    printf("\n\n");
	    return 1;
    }
}
    printf("efp2 ecd.                          : %.4f[ms]\n",ecd_time/ecd);
    printf("efp2 ecd projective lazy.          : %.4f[ms]\n",ecd_projective_lazy_time/ecd);
    printf("efp2 ecd jacobian lazy montgomery. : %.4f[ms]\n",ecd_jacobian_lazy_montgomery_time/ecd);


printf("------------------------------------------------------------------------------------\n");
    printf("efp2_eca test\n");

    efp2_rational_point(&A_efp2);
	efp2_to_montgomery(&A_efp2m,&A_efp2);

    efp2_affine_to_projective(&A_efp2_projective,&A_efp2);
    efp2_affine_to_jacobian(&A_efp2_jacobian,&A_efp2);
    efp2_affine_to_jacobian_montgomery(&A_efp2_jacobianm,&A_efp2m);
for(i=0;i<eca;i++){

    efp2_rational_point(&B_efp2);
	efp2_to_montgomery(&B_efp2m,&B_efp2);

    efp2_affine_to_projective(&B_efp2_projective,&B_efp2);
    efp2_affine_to_jacobian(&B_efp2_jacobian,&B_efp2);
    efp2_affine_to_jacobian_montgomery(&B_efp2_jacobianm,&B_efp2m);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    efp2_eca(&test1,&A_efp2,&B_efp2);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&eca_cost,&tmp);

    cost_zero();
    gettimeofday(&tv_A,NULL);
    //efp2_eca_projective_lazy(&testP2,&A_efp2_projective,&B_efp2_projective);
    gettimeofday(&tv_B,NULL);
    eca_projective_lazy_time+=timedifference_msec(tv_A,tv_B);
    cost_check(&tmp);
    cost_addition(&eca_projective_lazy_cost,&tmp);
    efp2_projective_to_affine(&test2,&testP2);

    gettimeofday(&tv_A,NULL);
    efp2_eca_jacobian_lazy_montgomery(&testZ5,&B_efp2_jacobianm,&A_efp2_jacobianm);
    gettimeofday(&tv_B,NULL);
    eca_jacobian_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    efp2_jacobian_to_affine_montgomery(&test5,&testZ5);
    efp2_mod_montgomery(&test5,&test5);

    gettimeofday(&tv_A,NULL);
    efp2_eca_mixture_lazy_montgomery(&testZ6,&B_efp2_jacobianm,&A_efp2_jacobianm);
    gettimeofday(&tv_B,NULL);
    eca_mixture_lazy_montgomery_time+=timedifference_msec(tv_A,tv_B);
    efp2_jacobian_to_affine_montgomery(&test6,&testZ6);
    efp2_mod_montgomery(&test6,&test6);

    if(efp2_cmp(&test1,&test2)!=0 || efp2_cmp(&test1,&test5)!=0 || efp2_cmp(&test1,&test6)!=0){
        printf("failed!\n\n");
	    efp2_printf("test1=",&test1);
	    efp2_printf("\ntest2=",&test2);
	    efp2_printf("\ntest5=",&test5);
	    efp2_printf("\ntest6=",&test6);
	    //efp2_jacobian_printf("\ntestZ2=",&testZ2);
	    //efp2_jacobian_printf("\ntestZ3=",&testZ3);
	    printf("\n\n");
	    return 1;
    }
}
    printf("efp2 eca.                          : %.4f[ms]\n",eca_time/eca);
    printf("efp2 eca projective lazy.          : %.4f[ms]\n",eca_projective_lazy_time/eca);
    printf("efp2 eca jacobian lazy montgomery. : %.4f[ms]\n",eca_jacobian_lazy_montgomery_time/eca);
    printf("efp2 eca mixture lazy montgomery.  : %.4f[ms]\n",eca_mixture_lazy_montgomery_time/eca);

    #ifdef DEBUG_COST_A
    printf("*********bls12 fp2 mul fp COST.********         \n");
    cost_printf("efp2 eca lazy cost",&eca_cost,eca);
    cost_printf("efp2 eca projective lazy cost",&eca_projective_lazy_cost,eca);
    printf("***************************************         \n");
    #endif


    return 0;
}
int test_efp12(int ecd,int eca,int scm){
    int i,n=0;
    float ecd_time=0,ecd_lazy_time=0;
    float eca_time=0,eca_lazy_time=0;
    float scm_time=0,scm_lazy_time=0;
    struct timeval tv_A,tv_B;
    printf("====================================================================================\n");
    efp12_t A_efp12,B_efp12,test1,test2;
    efp12_init(&A_efp12);
    efp12_init(&B_efp12);
    efp12_init(&test1);
    efp12_init(&test2);
    mpz_t scalar;
    mpz_init(scalar);

    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));


printf("------------------------------------------------------------------------------------\n");
    printf("efp12_ecd test\n");

for(i=0;i<ecd;i++){

    efp12_rational_point(&B_efp12);

    gettimeofday(&tv_A,NULL);
    efp12_ecd(&test1,&B_efp12);
    gettimeofday(&tv_B,NULL);
    ecd_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    efp12_ecd_lazy(&test2,&B_efp12);
    gettimeofday(&tv_B,NULL);
    ecd_lazy_time+=timedifference_msec(tv_A,tv_B);

    if(efp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    efp12_printf("",&test1);
	    efp12_printf("\n",&test2);
	printf("\n\n");
	return 1;
    }
}
    printf("efp12 ecd.      : %.4f[ms]\n",ecd_time/ecd);
    printf("efp12 ecd lazy. : %.4f[ms]\n",ecd_lazy_time/ecd);

printf("------------------------------------------------------------------------------------\n");
    printf("efp12_eca test\n");

    efp12_rational_point(&A_efp12);
for(i=0;i<eca;i++){

    efp12_rational_point(&B_efp12);

    gettimeofday(&tv_A,NULL);
    efp12_eca(&test1,&A_efp12,&B_efp12);
    gettimeofday(&tv_B,NULL);
    eca_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    efp12_eca_lazy(&test2,&A_efp12,&B_efp12);
    gettimeofday(&tv_B,NULL);
    eca_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(efp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    efp12_printf("",&test1);
	    efp12_printf("\n",&test2);
	printf("\n\n");
	return 1;
    }
}
    printf("efp12 eca.      : %.4f[ms]\n",eca_time/eca);
    printf("efp12 eca lazy. : %.4f[ms]\n",eca_lazy_time/eca);

printf("------------------------------------------------------------------------------------\n");
    printf("efp12_scm test\n");

for(i=0;i<scm;i++){
    mpz_urandomm(scalar,state,order_z);
    efp12_rational_point(&A_efp12);

    gettimeofday(&tv_A,NULL);
    efp12_scm(&test1,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_time+=timedifference_msec(tv_A,tv_B);

    gettimeofday(&tv_A,NULL);
    efp12_scm_lazy(&test2,&A_efp12,scalar);
    gettimeofday(&tv_B,NULL);
    scm_lazy_time+=timedifference_msec(tv_A,tv_B);


    if(efp12_cmp(&test1,&test2)!=0){
        printf("failed!\n\n");
	    efp12_printf("",&test1);
	    efp12_printf("\n",&test2);
	printf("\n\n");
	return 1;
    }
}
    printf("efp12 scm.      : %.4f[ms]\n",scm_time/scm);
    printf("efp12 scm lazy. : %.4f[ms]\n",scm_lazy_time/scm);

	return 0;
}
void test_Frobenius_map(){
    printf("====================================================================================\n");
    fp12_t A_fp12,test1,test2;
    fp12_init(&A_fp12);
    fp12_init(&test1);
    fp12_init(&test2);
    mpz_t exp;
    mpz_init(exp);

    gmp_randinit_default (state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    fp12_set_random(&A_fp12,state);
    fp12_printf("",&A_fp12);
    printf("\n\n");

    printf("frobenius\n");
    fp12_frobenius_map_p10(&test1,&A_fp12);
    fp12_printf("",&test1);
    printf("\n\n");

    mpz_pow_ui(exp,prime_z,10);
    fp12_pow(&test2,&A_fp12,exp);
    fp12_printf("",&test2);
    printf("\n");

    if(fp12_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }

    mpz_clear(exp);
}

void test_skew_frobenius_map(){
    printf("====================================================================================\n");
    efp12_t Q,test1,test2;
    efp12_init(&Q);
    efp12_init(&test1);
    efp12_init(&test2);
    efp2_t twisted_Q;
    efp2_init(&twisted_Q);

    bls12_generate_g2(&Q);
    efp12_to_efp2(&twisted_Q,&Q);

    fp12_frobenius_map_p10(&test1.x,&Q.x);
    fp12_frobenius_map_p10(&test1.y,&Q.y);
    efp12_printf("",&test1);
    printf("\n\n");

    efp2_skew_frobenius_map_p10(&twisted_Q,&twisted_Q);
    efp2_to_efp12(&test2,&twisted_Q);
    efp12_printf("",&test2);
    printf("\n\n");

    if(fp12_cmp(&test1.x,&test2.x)==0 && fp12_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
}


void test_twist(){
    printf("====================================================================================\n");
    efp12_t Q,test1,test2;
    efp12_init(&Q);
    efp12_init(&test1);
    efp12_init(&test2);
    efp2_t twist_Q;
    efp2_init(&twist_Q);


    bls12_generate_g2(&Q);
    efp12_printf("Q\n",&Q);
    printf("\n\n");

    efp12_to_efp2(&twist_Q,&Q);
    efp2_ecd(&twist_Q,&twist_Q);
    efp2_to_efp12(&test1,&twist_Q);
    efp12_printf("",&test1);
    printf("\n\n");

    efp12_ecd(&test2,&Q);
    efp12_printf("",&test2);
    printf("\n\n");

    if(fp12_cmp(&test1.x,&test2.x)==0 && fp12_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }

}
void test_All(){
	int point,Field,efp,efp2,efp12,G1,G2,G3,pairing;

printf("====================================================================================\n");
    printf("Test time\n");

	point = bls12_test_rational_point();
	//Field = test_field(100,100,100,10,10);
	//efp = test_efp(100,100,10);
	//efp2 = test_efp2(10,10,10);
	//efp12 = test_efp12(10,10,10);
	G1 = bls12_test_g1_scm(100);
	G2 = bls12_test_g2_scm(10);
	G3 = bls12_test_g3_exp(100);
	pairing = bls12_test_opt_ate_pairing(100);

printf("====================================================================================\n");
    printf("Test Result\n\n");

	printf("Test Point      :");
	if(point==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test Field      :");
	if(Field==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test efp        :");
	if(efp==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test efp2       :");
	if(efp2==0) printf("Success\n");
	else printf("Failed\n");

	printf("Test efp12      :");
	if(efp12==0) printf("Success\n");
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

void bench_fp(int scm){

}
