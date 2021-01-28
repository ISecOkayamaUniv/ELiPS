#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    bls12_init();
    bls12_print_parameters();
    fr_order_init();
    fr_t s_128,s_256,s;
    g1_t P_random,P_ans;
    int cnt=100000;
    float time_128,time_256;
    fr_init(&s);
    fr_init(&s_128);
    fr_init(&s_256);
    g1_init(&P_random);
    g1_init(&P_ans);

    g1_set_random(&P_random,state);
    fr_set_random(&s,state);
    for(int i=0;i<4;i++){
        s_256.x0[i]=s.x0[i];
    }
    for(int i=0;i<2;i++){
        s_128.x0[i]=s_256.x0[i];
    }
    fr_println("s_128=",&s_128);
    fr_println("s_256=",&s_256);
    gettimeofday(&tv_start,NULL);
    for(int i=0;i<cnt;i++){
        g1_scm(&P_ans,&P_random,&s_128);
    }
    gettimeofday(&tv_end,NULL);
    time_128=timedifference_msec(tv_start,tv_end);
    printf("s_128[P] time: %.4f[ms]\n",time_128/cnt);

    gettimeofday(&tv_start,NULL);
    for(int i=0;i<cnt;i++){
        g1_scm(&P_ans,&P_random,&s_256);
    }
    gettimeofday(&tv_end,NULL);
    time_256=timedifference_msec(tv_start,tv_end);
    printf("s_256[P] time: %.4f[ms]\n",time_256/cnt);
    return 0;
}
