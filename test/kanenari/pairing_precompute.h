#include "ELiPS/bls12.h"
#include "ELiPS/fr.h"

fp2_t miller_table[3][10000];
float MILLER_OPT_PROJECTIVE_PRECOMPUTE,FINALEXP_OPT_PROJECTIVE_PRECOMPUTE;
cost MILLER_OPT_PROJECTIVE_PRECOMPUTE_COST,FINALEXP_OPT_PROJECTIVE_PRECOMPUTE_COST;


void ff_ltt_precompute(efp2_projective_t *T,int i){
    static fp2_t tmp0_fp2,tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;
    static fpd2_t u0, u1;

    fp2_sqr_lazy_montgomery(&tmp0_fp2, &T->z);
    fp2_sqr_lazy_montgomery(&tmp1_fp2, &T->y);
    fp2_add_nonmod_single(&tmp5_fp2, &tmp0_fp2, &tmp1_fp2);
    fp2_mul_3_twist_b(&tmp2_fp2,&tmp0_fp2);

    fp2_sqr_lazy_montgomery(&miller_table[0][i], &T->x);

    fp2_mul_lazy_montgomery(&tmp4_fp2, &T->x, &T->y);
    fp_r1shift(&tmp4_fp2.x0, &tmp4_fp2.x0);
    fp_r1shift(&tmp4_fp2.x1, &tmp4_fp2.x1);
    fp2_add_nonmod_single(&tmp3_fp2, &tmp2_fp2, &tmp2_fp2);
    fp2_add_nonmod_single(&tmp3_fp2, &tmp3_fp2, &tmp2_fp2);
    fp2_sub_nonmod_single(&T->x, &tmp1_fp2, &tmp3_fp2);
    fp2_mul_lazy_montgomery(&T->x, &T->x, &tmp4_fp2);

    fp2_add_nonmod_single(&tmp3_fp2, &tmp1_fp2, &tmp3_fp2);
    fp_r1shift(&tmp3_fp2.x0, &tmp3_fp2.x0);
    fp_r1shift(&tmp3_fp2.x1, &tmp3_fp2.x1);

    //Lazy
    fp2_sqr_nonmod_montgomery(&u0, &tmp2_fp2);
    fp2_add_nonmod_double(&u1, &u0, &u0);
    fp2_add_nonmod_double(&u1, &u1, &u0);
    fp2_sqr_nonmod_montgomery(&u0, &tmp3_fp2);
    fp2_sub_nonmod_double(&u0, &u0, &u1);

    fp2_add_nonmod_single(&tmp3_fp2, &T->y, &T->z);
    fp2_sqr_lazy_montgomery(&tmp3_fp2, &tmp3_fp2);
    fp2_sub_nonmod_single(&miller_table[1][i], &tmp3_fp2, &tmp5_fp2);

    fp2_mod_montgomery_double(&T->y, &u0);

    fp2_mul_lazy_montgomery(&T->z, &tmp1_fp2, &miller_table[1][i]);
#ifdef TWIST_PHI_INV
    fp2_sub_nonmod_single(&miller_table[2][i], &tmp2_fp2, &tmp1_fp2);
#endif
#ifdef TWIST_PHI
    fp2_sub_nonmod_single(&miller_table[2][i], &tmp2_fp2, &tmp1_fp2);
#endif
}
void f_ltq_precompute(efp2_projective_t *T,efp2_projective_t *Q,int i){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12;

    fp2_mul_lazy_montgomery(&tmp1_fp2, &T->z, &Q->x);
    fp2_sub_nonmod_single(&miller_table[2][i], &T->x, &tmp1_fp2);
    fp2_mul_lazy_montgomery(&tmp2_fp2, &T->z, &Q->y);
    fp2_sub_nonmod_single(&miller_table[0][i], &T->y, &tmp2_fp2);

    fp2_sqr_lazy_montgomery(&tmp3_fp2, &miller_table[2][i]);
    fp2_mul_lazy_montgomery(&T->x, &T->x, &tmp3_fp2);
    fp2_mul_lazy_montgomery(&tmp3_fp2, &tmp3_fp2, &miller_table[2][i]);
    fp2_sqr_lazy_montgomery(&tmp4_fp2, &miller_table[0][i]);
    fp2_mul_lazy_montgomery(&tmp4_fp2, &tmp4_fp2, &T->z);
    fp2_add_nonmod_single(&tmp4_fp2, &tmp3_fp2, &tmp4_fp2);


    fp2_mul_lazy_montgomery(&tmp5_fp2, &Q->x, &miller_table[0][i]);

    fp2_sub_nonmod_single(&tmp4_fp2, &tmp4_fp2, &T->x);
    fp2_sub_nonmod_single(&tmp4_fp2, &tmp4_fp2, &T->x);
    fp2_sub_nonmod_single(&T->x, &T->x, &tmp4_fp2);
    fp2_mul_lazy_montgomery(&tmp2_fp2, &miller_table[0][i], &T->x);
    fp2_mul_lazy_montgomery(&T->y, &tmp3_fp2, &T->y);
    fp2_sub_nonmod_single(&T->y, &tmp2_fp2, &T->y);
    fp2_mul_lazy_montgomery(&T->x, &miller_table[2][i], &tmp4_fp2);
    fp2_mul_lazy_montgomery(&T->z, &T->z, &tmp3_fp2);

    fp2_mul_lazy_montgomery(&tmp3_fp2, &Q->y, &miller_table[2][i]);
    #ifdef TWIST_PHI_INV
    fp2_sub_nonmod_single(&miller_table[1][i], &tmp5_fp2, &tmp3_fp2);
    #endif
    #ifdef TWIST_PHI
    fp2_sub_nonmod_single(&miller_table[1][i], &tmp5_fp2, &tmp3_fp2);
    #endif
}

void ff_ltt_using_precompute(fp12_t *f,efp_t *P,int i){
    //static fp2_t tmp0_fp2,tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2,buf_fp2;
    static fp2_t tmp0_fp2,tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12,tmp2_fp12;

    static fpd2_t u0, u1;

    fp12_sqr_lazy_montgomery(&tmp1_fp12, f);

#ifdef TWIST_PHI_INV
    fp2_set(&tmp2_fp12.x1.x1, &miller_table[2][i]);
    fp_mulmod_montgomery(&tmp2_fp12.x1.x0.x0, &miller_table[0][i].x0, &P->x);
    fp_mulmod_montgomery(&tmp2_fp12.x1.x0.x1, &miller_table[0][i].x1, &P->x);
    fp_mulmod_montgomery(&tmp2_fp12.x0.x0.x0, &miller_table[1][i].x0, &P->y);
    fp_mulmod_montgomery(&tmp2_fp12.x0.x0.x1, &miller_table[1][i].x1, &P->y);
#endif
#ifdef TWIST_PHI
    fp2_set(&tmp2_fp12.x0.x0, &miller_table[2][i]);
    fp_mulmod_montgomery(&tmp2_fp12.x0.x1.x0, &miller_table[0][i].x0, &P->x);
    fp_mulmod_montgomery(&tmp2_fp12.x0.x1.x1, &miller_table[0][i].x1, &P->x);
    fp_mulmod_montgomery(&tmp2_fp12.x1.x1.x0, &miller_table[1][i].x0, &P->y);
    fp_mulmod_montgomery(&tmp2_fp12.x1.x1.x1, &miller_table[1][i].x1, &P->y);
#endif
    fp12_6_sparse_mul_lazy_montgomery(f, &tmp1_fp12, &tmp2_fp12);
}
void f_ltq_using_precompute(fp12_t *f,efp_t *P,int i){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2,tmp5_fp2;
    static fp12_t tmp1_fp12;
    #ifdef TWIST_PHI_INV
    fp_mulmod_montgomery(&tmp1_fp12.x1.x0.x0, &miller_table[0][i].x0, &P->x);
    fp_mulmod_montgomery(&tmp1_fp12.x1.x0.x1, &miller_table[0][i].x1, &P->x);
    #endif
    #ifdef TWIST_PHI
    fp_mulmod_montgomery(&tmp1_fp12.x0.x1.x0, &miller_table[0][i].x0, &P->x);
    fp_mulmod_montgomery(&tmp1_fp12.x0.x1.x1, &miller_table[0][i].x1, &P->x);
    #endif
    #ifdef TWIST_PHI_INV
    fp2_set(&tmp1_fp12.x1.x1, &miller_table[1][i]);
    fp_mulmod_montgomery(&tmp1_fp12.x0.x0.x0, &miller_table[2][i].x0, &P->y);
    fp_mulmod_montgomery(&tmp1_fp12.x0.x0.x1, &miller_table[2][i].x1, &P->y);
    #endif
    #ifdef TWIST_PHI
    fp2_set(&tmp1_fp12.x0.x0, &miller_table[1][i]);
    fp_mulmod_montgomery(&tmp1_fp12.x1.x1.x0, &miller_table[2][i].x0, &P->y);
    fp_mulmod_montgomery(&tmp1_fp12.x1.x1.x1, &miller_table[2][i].x1, &P->y);
    #endif
    fp12_6_sparse_mul_lazy_montgomery(f,f,&tmp1_fp12);
}
void g1g2_to_g3_pairing_precompute_init(g2_t *Q){
    efp2_projective_t T;
    efp2_projective_t mapped_Q,mapped_Q_neg;
    int i;

    efp2_projective_init(&T);
    efp2_projective_init(&mapped_Q);
    efp2_projective_init(&mapped_Q_neg);
    efp2_affine_to_projective_montgomery(&mapped_Q, Q);
    efp2_projective_set(&mapped_Q_neg, &mapped_Q); //set mapped_Q_neg
    fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    efp2_projective_set(&T, &mapped_Q); //set T
    int cnt=0;
    #ifdef X_PLUS
    //miller
    for(i=bls12_X_length-1; i>=0; i--){
        switch(bls12_X_binary[i]){
            case 0:
                ff_ltt_precompute(&T,cnt);cnt++;
                break;
            case 1:
                ff_ltt_precompute(&T,cnt);cnt++;
                f_ltq_precompute(&T,&mapped_Q,cnt);cnt++;
                break;
            case -1:
                ff_ltt_precompute(&T,cnt);cnt++;
                f_ltq_precompute(&T,&mapped_Q_neg,cnt); cnt++;
                break;
            default:
                break;
        }
    }
    #endif
    #ifdef X_MINUS
    //miller
    for(i=bls12_X_length-1; i>=0; i--){
        switch(bls12_X_binary[i]){
            case 0:
                ff_ltt_precompute(&T,i);cnt++;
                break;
            case -1:
                ff_ltt_precompute(&T,i);cnt++;
                f_ltq_precompute(&T,&mapped_Q,i); cnt++;
                break;
            case 1:
                ff_ltt_precompute(&T,i);cnt++;
                f_ltq_precompute(&T,&mapped_Q_neg,i); cnt++;
                break;
            default:
                break;
        }
    }
    #endif
}
void g1g2_to_g3_miller_algo_using_precompute(g3_t *ANS,g1_t *P){
    efp_t mapped_P_ltt,mapped_P_ltq;
    fp12_t f;
    int i;
    efp_init(&mapped_P_ltt);
    efp_init(&mapped_P_ltq);
    fp12_init(&f);

    //set
    efp_set(&mapped_P_ltq,P);//set mapped_P
    fp_set_mpn(&f.x0.x0.x0,RmodP);

    //precompute P
    fp_add_nonmod_single(&mapped_P_ltt.x,&mapped_P_ltq.x,&mapped_P_ltq.x);
    fp_add_nonmod_single(&mapped_P_ltt.x,&mapped_P_ltt.x,&mapped_P_ltq.x);
    fp_set_neg(&mapped_P_ltt.y,&mapped_P_ltq.y);
    fp_set_neg(&mapped_P_ltq.x,&mapped_P_ltq.x);
    int cnt=0;

    #ifdef X_PLUS
    //miller
    for(i=bls12_X_length-1; i>=0; i--){
        switch(bls12_X_binary[i]){
            case 0:
                ff_ltt_using_precompute(&f,&mapped_P_ltt,cnt);cnt++;
                break;
            case 1:
                ff_ltt_using_precompute(&f,&mapped_P_ltt,cnt);cnt++;
                f_ltq_using_precompute(&f,&mapped_P_ltq,cnt);cnt++;
                break;
            case -1:
                ff_ltt_using_precompute(&f,&mapped_P_ltt,cnt);cnt++;
                f_ltq_using_precompute(&f,&mapped_P_ltq,cnt);cnt++;
                break;
            default:
                break;
        }
    }
    #endif
    #ifdef X_MINUS
    //miller
    for(i=bls12_X_length-1; i>=0; i--){
        switch(bls12_X_binary[i]){
            case 0:
                ff_ltt_using_precompute(&f,&mapped_P_ltt,i);cnt++;
                break;
            case -1:
                ff_ltt_using_precompute(&f,&mapped_P_ltt,i);cnt++;
                f_ltq_using_precompute(&f,&mapped_P_ltq,i);cnt++;
                break;
            case 1:
                ff_ltt_using_precompute(&f,&mapped_P_ltt,i);cnt++;
                f_ltq_using_precompute(&f,&mapped_P_ltq,i);cnt++;
                break;
            default:
                break;
        }
    }
    #endif
    fp12_set(ANS,&f);
}
void g1g2_to_g3_pairing_precompute(g3_t *ANS,g1_t *P){
    #ifdef DEBUG_COST_A
    cost tmp;
    #endif

    //Miller's Algo.
    gettimeofday(&tv_start,NULL);
    g1g2_to_g3_miller_algo_using_precompute(ANS,P);
    gettimeofday(&tv_end,NULL);
    MILLER_OPT_PROJECTIVE_PRECOMPUTE+=timedifference_msec(tv_start,tv_end);

    #ifdef DEBUG_COST_A
    cost_check(&tmp);
    cost_addition(&MILLER_OPT_PROJECTIVE_PRECOMPUTE_COST,&tmp);
    #endif

    //Final Exp.
    gettimeofday(&tv_start,NULL);
    g1g2_to_g3_final_exp(ANS,ANS);
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_PROJECTIVE_PRECOMPUTE+=timedifference_msec(tv_start,tv_end);
}
int bench_pairing(int pairing){
    int i,n=0;
    float opt_time=0,opt_affine_time=0,opt_precompute_time=0;
    cost tmp,opt_cost,opt_affine_cost;
    struct timeval tv_A,tv_B;

    g1_t P;
    g2_t Q;
    g1_init(&P);
    g2_init(&Q);
    g3_t test;


    printf("====================================================================================\n");
    printf("bls12_Opt-ate pairing\n\n");

    MILLER_OPT_PROJECTIVE = 0;
    FINALEXP_OPT_PROJECTIVE = 0;
    MILLER_OPT_AFFINE = 0;
    FINALEXP_OPT_AFFINE = 0;
    MILLER_OPT_PROJECTIVE_PRECOMPUTE = 0;
    FINALEXP_OPT_PROJECTIVE_PRECOMPUTE = 0;
    opt_time=0;
    opt_affine_time=0;
    cost_init(&MILLER_OPT_PROJECTIVE_COST);
    cost_init(&FINALEXP_OPT_PROJECTIVE_COST);
    cost_init(&MILLER_OPT_AFFINE_COST);
    cost_init(&FINALEXP_OPT_AFFINE_COST);
    cost_init(&MILLER_OPT_PROJECTIVE_PRECOMPUTE_COST);
    cost_init(&FINALEXP_OPT_PROJECTIVE_PRECOMPUTE_COST);

    g1_set_random(&P,state);
    g2_set_random(&Q,state);

    //init
    g1g2_to_g3_pairing_precompute_init(&Q);
    printf("1 ok\n");
    for(i=0;i<pairing;i++){
        g1_set_random_with_basepoint(&P,&P,state);
        //g2_set_random_with_basepoint(&Q,&Q,state);
        //g1g2_to_g3_pairing(&test1,&P,&Q);
        cost_zero();
        gettimeofday(&tv_A,NULL);
        g1g2_to_g3_pairing(&test,&P,&Q);
        gettimeofday(&tv_B,NULL);
        opt_time+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&opt_cost,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        g1g2_to_g3_pairing_affine(&test,&P,&Q);
        gettimeofday(&tv_B,NULL);
        opt_affine_time+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&opt_affine_cost,&tmp);

        cost_zero();
        gettimeofday(&tv_A,NULL);
        g1g2_to_g3_pairing_precompute(&test,&P);
        gettimeofday(&tv_B,NULL);
        opt_precompute_time+=timedifference_msec(tv_A,tv_B);
        cost_check(&tmp);
        cost_addition(&opt_affine_cost,&tmp);
    }
    cost_substruction(&FINALEXP_OPT_PROJECTIVE_COST, &opt_cost, &MILLER_OPT_PROJECTIVE_COST);
    cost_substruction(&FINALEXP_OPT_AFFINE_COST, &opt_affine_cost, &MILLER_OPT_AFFINE_COST);
    cost_substruction(&FINALEXP_OPT_PROJECTIVE_PRECOMPUTE_COST, &opt_cost, &MILLER_OPT_PROJECTIVE_PRECOMPUTE_COST);


    printf("bls12 opt ate                    : %.4f[ms]\n",opt_time/pairing);
    printf("bls12 opt ate (MILLER).          : %.4f[ms]\n",MILLER_OPT_PROJECTIVE/pairing);
    printf("bls12 opt ate (FINALEXP).        : %.4f[ms]\n\n",FINALEXP_OPT_PROJECTIVE/pairing);

    printf("bls12 opt ate affine.            : %.4f[ms]\n",opt_affine_time/pairing);
    printf("bls12 opt ate (MILLER).          : %.4f[ms]\n",MILLER_OPT_AFFINE/pairing);
    printf("bls12 opt ate (FINALEXP).        : %.4f[ms]\n\n",FINALEXP_OPT_AFFINE/pairing);

    printf("bls12 opt ate precompute.        : %.4f[ms]\n",opt_precompute_time/pairing);
    printf("bls12 opt ate (MILLER).          : %.4f[ms]\n",MILLER_OPT_PROJECTIVE_PRECOMPUTE/pairing);
    printf("bls12 opt ate (FINALEXP).        : %.4f[ms]\n",FINALEXP_OPT_PROJECTIVE_PRECOMPUTE/pairing);

    #ifdef DEBUG_COST_A
    printf("*********bls12 opt ate fp COST.********         \n");
    cost_printf("bls12 opt ate ", &opt_cost, pairing);
    cost_printf("bls12 opt ate (MILLER)", &MILLER_OPT_PROJECTIVE_COST, pairing);
    cost_printf("bls12 opt ate (FINALEXP)", &FINALEXP_OPT_PROJECTIVE_COST, pairing);
    printf("***************************************         \n");
    printf("*********bls12 opt ate affine fp COST.********         \n");
    cost_printf("bls12 opt ate ", &opt_affine_cost, pairing);
    cost_printf("bls12 opt ate (MILLER)", &MILLER_OPT_AFFINE_COST, pairing);
    cost_printf("bls12 opt ate (FINALEXP)", &FINALEXP_OPT_AFFINE_COST, pairing);
    printf("***************************************         \n");
    printf("*********bls12 opt ate precompute fp COST.********         \n");
    cost_printf("bls12 opt ate ", &opt_precompute_cost, pairing);
    cost_printf("bls12 opt ate (MILLER)", &MILLER_OPT_PROJECTIVE_PRECOMPUTE_COST, pairing);
    cost_printf("bls12 opt ate (FINALEXP)", &FINALEXP_OPT_PROJECTIVE_PRECOMPUTE_COST, pairing);
    printf("***************************************         \n");
    #endif

    return 0;
}
// int miller_test(){
//     g1_t P;
//     g2_t Q;
//     g3_t test1,test2;

//     g1_set_random(&P,state);
//     g2_set_random(&Q,state);
//     g1g2_to_g3_miller_algo(&test1,&P,&Q);
//     fp12_mod_montgomery(&test1,&test1);
//     g1g2_to_g3_pairing_precompute_init(&Q);
//     //getchar();
//     printf("\nprecompute start\n");
//     g1g2_to_g3_miller_algo_using_precompute(&test2,&P);
//     fp12_mod_montgomery(&test2,&test2);

//     if(fp12_cmp(&test1,&test2) != 0){
//         printf("pairing projective failed!\n\n");
//         printf("\n\n");
//         fp12_printf_montgomery("\ntest1\n",&test1);
//         fp12_printf_montgomery("\ntest2\n",&test2);
//         return 1;
//     }else{
//         printf("ok!");
//     }
// }
void billinear_test_precompute(){
    fr_t a,b,c;
    g1_t P,aP;
    g2_t Q,bQ;
    g3_t test1,test2;

    fr_set_random(&a,state);

    g1_set_random(&P,state);
    g2_set_random(&Q,state);
    g1_scm(&aP,&P,&a);

   g1g2_to_g3_pairing_precompute_init(&Q);

    //test1={Pairng(P,Q)^c
    g1g2_to_g3_pairing_precompute(&test1,&P);
    g3_exp(&test1,&test1,&a);

    //test2=Pairng(aP,bQ)
    g1g2_to_g3_pairing_precompute(&test2,&aP);

    if(fp12_cmp(&test1,&test2)==0) printf("ok!\n");
    else printf("ng\n");
}
// int main(){
//     bls12_init();
//     fr_order_init();
//     bench_pairing(100000);
//     //miller_test();
//     //billinear_test_precompute();
//     return 0;
// }
