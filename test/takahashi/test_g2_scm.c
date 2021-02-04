#include "ELiPS/bls12.h"
#include "test_g2_scm.h"

/*============================================================================*/
/* main                                                                       */
/*==========================================s==================================*/
int main(void){
    //Declarate variable
    scalar_t s;
    efp12_t Q;
    efp12_t ans0,ans1,ans2,ans3;
    efp12_t ans_wnaf_1,ans_wnaf_3,ans_wnaf_5,ans_wnaf_7,ans_jsf;
    float test0,test1,test2;
    float test_wnaf_3,test_wnaf_5,test_wnaf_7;
	float test_jsf;
    struct timeval tv_A,tv_B;

    //Initialize variable
    scalar_init(s);
    efp12_init(&Q);
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

    int cnt = 10000;
     //Generate rationalpoit P on g2 and Q on G2
    bls12_generate_g2(&Q);

    cost cost_g2_basic,cost_g2,cost_naf_3,cost_naf_5,cost_naf_7,cost_jsf;
    cost tmp;

    for(int i=0;i<cnt;i++){
        //Generate random scalar
        scalar_random_order(s);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g2_scm_basic(&ans0,&Q,s);
        gettimeofday(&tv_B,NULL);
        test0+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_g2_basic,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g2_scm(&ans1,&Q,s);
        gettimeofday(&tv_B,NULL);
        test1+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_g2,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g2_scm_w_naf(&ans_wnaf_3,&Q,s,3);
        gettimeofday(&tv_B,NULL);
        test_wnaf_3+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_naf_3,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g2_scm_w_naf(&ans_wnaf_5,&Q,s,5);
        gettimeofday(&tv_B,NULL);
        test_wnaf_5+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_naf_5,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        test_g2_scm_w_naf(&ans_wnaf_7,&Q,s,7);
        gettimeofday(&tv_B,NULL);
        test_wnaf_7+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&cost_naf_7,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
		test_g2_scm_jsf(&ans_jsf,&Q,s);
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
            printf("test jsf    failed!\n\n");
            return 1;
        }
    }

    //Compare answer
    printf("**********result***********\n");
    printf("test g2 scm basic      : %.4f[ms]\n",test0/cnt);
    printf("test g2 scm            : %.4f[ms]\n",test1/cnt);
    printf("test g2 scm naf 3      : %.4f[ms]\n",test_wnaf_3/cnt);
    printf("test g2 scm naf 5      : %.4f[ms]\n",test_wnaf_5/cnt);
    printf("test g2 scm naf 7      : %.4f[ms]\n",test_wnaf_7/cnt);
	printf("test g2 scm jsf        : %.4f[ms]\n",test_jsf/cnt);

    #ifdef DEBUG_COST_A
    printf("*********COST********         \n");
    cost_printf("test g2 scm basic", &cost_g2_basic, cnt);
    cost_printf("test g2 scm      ", &cost_g2, cnt);
    cost_printf("test g2 scm naf 3", &cost_naf_3, cnt);
    cost_printf("test g2 scm naf 5", &cost_naf_5, cnt);
    cost_printf("test g2 scm naf 7", &cost_naf_7, cnt);
    cost_printf("test g2 scm jsf  ", &cost_jsf, cnt);
    printf("***************************************         \n");
    #endif
    //Clear variable
    scalar_clear(s);

    return 0;
}
