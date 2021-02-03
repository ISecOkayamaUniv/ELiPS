#include "ELiPS/bls12.h"
#include "test_g3_exp.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    //Declarate variable
    scalar_t s;

    efp12_t P,Q;
    fp12_t E;
    fp12_t ans0,ans1,ans2,ans3;
    fp12_t ans_wnaf_1,ans_wnaf_3,ans_wnaf_5,ans_wnaf_7;
    float test0,test1,test2;
    float test_wnaf_3,test_wnaf_5,test_wnaf_7;
    struct timeval tv_A,tv_B;
    
    //Initialize variable
    scalar_init(s);
    efp12_init(&Q);
    fp12_init(&ans0);
    fp12_init(&ans1);
    fp12_init(&ans2);
    fp12_init(&ans3);
    fp12_init(&ans_wnaf_1);
    fp12_init(&ans_wnaf_3);
    fp12_init(&ans_wnaf_5);
    fp12_init(&ans_wnaf_7);

    //Initialize and print parameter
    bls12_init();
    bls12_print_parameters();
    
    int cnt = 10;
    
    for(int i=0;i<cnt;i++){
        //Generate random scalar
        scalar_random_order(s);
        
        //Generate rationalpoit P on G1 and Q on G2
        bls12_generate_g1(&P);
        bls12_generate_g2(&Q);
        bls12_optate_pairing(&E,&P,&Q);
        
        gettimeofday(&tv_A,NULL);
        test_g3_exp_basic(&ans0,&E,s);
        gettimeofday(&tv_B,NULL);
        test0+=timedifference_msec(tv_A,tv_B);
        
        gettimeofday(&tv_A,NULL);
        test_g3_exp(&ans1,&E,s);
        gettimeofday(&tv_B,NULL);
        test1+=timedifference_msec(tv_A,tv_B);
        
        gettimeofday(&tv_A,NULL);
        test_g3_exp_w_naf(&ans_wnaf_3,&E,s,3);
        gettimeofday(&tv_B,NULL);
        test_wnaf_3+=timedifference_msec(tv_A,tv_B);
        
        gettimeofday(&tv_A,NULL);
        test_g3_exp_w_naf(&ans_wnaf_5,&E,s,5);
        gettimeofday(&tv_B,NULL);
        test_wnaf_5+=timedifference_msec(tv_A,tv_B);

        gettimeofday(&tv_A,NULL);
        test_g3_exp_w_naf(&ans_wnaf_7,&E,s,7);
        gettimeofday(&tv_B,NULL);
        test_wnaf_7+=timedifference_msec(tv_A,tv_B);
        
        
        if(fp12_cmp(&ans0,&ans1) != 0){
            printf("test1 failed!\n\n");
            return 1;
        }
        
        if(fp12_cmp(&ans0,&ans_wnaf_3) != 0){
            printf("test wnaf 3 failed!\n\n");
            return 1;
        }
        if(fp12_cmp(&ans0,&ans_wnaf_5) != 0){
            printf("test wnaf 5 failed!\n\n");
            return 1;
        }
        if(fp12_cmp(&ans0,&ans_wnaf_7) != 0){
            printf("test wnaf 7 failed!\n\n");
            return 1;
        }
    }

    //Compare answer
    printf("**********result***********\n");
    printf("test g3 exp basic      : %.4f[ms]\n",test0/cnt);
    printf("test g3 exp            : %.4f[ms]\n",test1/cnt);
    printf("test g3 exp naf 3      : %.4f[ms]\n",test_wnaf_3/cnt);
    printf("test g3 exp naf 5      : %.4f[ms]\n",test_wnaf_5/cnt);
    printf("test g3 exp naf 7      : %.4f[ms]\n",test_wnaf_7/cnt);
    
    
    //Clear variable
    scalar_clear(s);
    
    return 0;
}
