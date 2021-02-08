#include "ELiPS/bls12.h"
#include "test_g1_scm.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    //Declarate variable
    scalar_t s;
    efp12_t P;
    efp12_t ans0,ans1,ans2,ans3;
    efp12_t ans_wnaf_1,ans_wnaf_3,ans_wnaf_5,ans_wnaf_7,ans_jsf;
    float test0,test1,test2;
    float test_wnaf_3,test_wnaf_5,test_wnaf_7;
    float test_jsf;

    struct timeval tv_A,tv_B;

    //Initialize variable
    scalar_init(s);
    efp12_init(&P);
    efp12_init(&ans0);
    efp12_init(&ans1);
    efp12_init(&ans2);
    efp12_init(&ans3);
    efp12_init(&ans_wnaf_1);
    efp12_init(&ans_wnaf_3);
    efp12_init(&ans_wnaf_5);
    efp12_init(&ans_wnaf_7);
    efp12_init(&ans_jsf);


    //Initialize and print parameter
    bls12_init();
    bls12_print_parameters();

    int cnt = 100;

    //Generate rationalpoit P on G1 and Q on G2
    bls12_generate_g1(&P);

    cost cost_g1_basic,cost_g1,cost_naf_3,cost_naf_5,cost_naf_7,cost_jsf;
    cost tmp;
    //cost_addition(&MILLER_OPT_AFFINE_COST,&tmp);


    for(int i=0;i<cnt;i++){
        //Generate random scalar
        scalar_random_order(s);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g1_scm_basic(&ans0,&P,s);
        gettimeofday(&tv_B,NULL);
        test0+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_g1_basic,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g1_scm(&ans1,&P,s);
        gettimeofday(&tv_B,NULL);
        test1+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_g1,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g1_scm_w_naf(&ans_wnaf_3,&P,s,3);
        gettimeofday(&tv_B,NULL);
        test_wnaf_3+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_naf_3,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g1_scm_w_naf(&ans_wnaf_5,&P,s,5);
        gettimeofday(&tv_B,NULL);
        test_wnaf_5+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_naf_5,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g1_scm_w_naf(&ans_wnaf_7,&P,s,7);
        gettimeofday(&tv_B,NULL);
        test_wnaf_7+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_naf_7,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
		test_g1_scm_jsf(&ans_jsf,&P,s);
		gettimeofday(&tv_B,NULL);
        test_jsf+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_jsf,&tmp);


        if(efp12_cmp(&ans0,&ans1) != 0){
            printf("test1 failed!\n\n");
            return 1;
        }

        if(efp12_cmp(&ans0,&ans_wnaf_3) != 0){
            printf("test wnaf 3 failed!\n\n");
            return 1;
        }
        if(efp12_cmp(&ans0,&ans_wnaf_5) != 0){
            printf("test wnaf 5 failed!\n\n");
            return 1;
        }
        if(efp12_cmp(&ans0,&ans_wnaf_7) != 0){
            printf("test wnaf 7 failed!\n\n");
            return 1;
        }
		if(efp12_cmp(&ans0,&ans_jsf) != 0){
            printf("test jsf failed!\n\n");
            return 1;
        }
    }

    //Compare answer
    printf("**********result***********\n");
    printf("test g1 scm basic      : %.4f[ms]\n",test0/cnt);
    printf("test g1 scm            : %.4f[ms]\n",test1/cnt);
    printf("test g1 scm naf 3      : %.4f[ms]\n",test_wnaf_3/cnt);
    printf("test g1 scm naf 5      : %.4f[ms]\n",test_wnaf_5/cnt);
    printf("test g1 scm naf 7      : %.4f[ms]\n\n",test_wnaf_7/cnt);
    printf("test g1 scm jsf        : %.4f[ms]\n",test_jsf/cnt);

    #ifdef DEBUG_COST_A
    printf("*********COST********         \n");
    cost_printf("test g1 scm basic", &cost_g1_basic, cnt);
    cost_printf("test g1 scm      ", &cost_g1, cnt);
    cost_printf("test g1 scm naf 3", &cost_naf_3, cnt);
    cost_printf("test g1 scm naf 5", &cost_naf_5, cnt);
    cost_printf("test g1 scm naf 7", &cost_naf_7, cnt);
    cost_printf("test g1 scm jsf  ", &cost_jsf, cnt);
    printf("***************************************         \n");
    #endif

    //Clear variable
    scalar_clear(s);

    return 0;
}
