#include <ELiPS/bls12_miller.h>

void bls12_miller_algo_for_optate_basic(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    efp12_t Buf;
    efp12_init(&Buf);
    efp2_t T;
    efp2_init(&T);
    efp2_t mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    efp2_init(&mapped_Q);
    efp2_init(&mapped_Q_neg);
    efp2_init(&mapped_Q1);
    efp2_init(&mapped_Q2_neg);
    efp_t mapped_P;
    efp_init(&mapped_P);
    fp12_t f;
    fp12_init(&f);
    fp_t L;
    fp_init(&L);
    int i;

    //set
    efp12_to_efp(&mapped_P,P);//set mapped_P
    efp12_to_efp2(&mapped_Q,Q);//set mapped_Q
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    efp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    efp2_set(&T,&mapped_Q);     //set T
    fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=bls12_X_length-1; i>=0; i--){
        switch(bls12_X_binary[i]){
            case 0:
                ff_ltt(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
    }

    fp12_set(ANS,&f);
}

void bls12_miller_algo_for_optate_affine(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    efp12_t Buf;
    efp12_init(&Buf);
    efp2_t T;
    efp2_init(&T);
    efp2_t mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    efp2_init(&mapped_Q);
    efp2_init(&mapped_Q_neg);
    efp2_init(&mapped_Q1);
    efp2_init(&mapped_Q2_neg);
    efp_t mapped_P;
    efp_init(&mapped_P);
    fp12_t f;
    fp12_init(&f);
    fp_t L;
    fp_init(&L);
    int i;

    //set
    efp12_to_efp(&mapped_P,P);//set mapped_P
    efp12_to_efp2(&mapped_Q,Q);//set mapped_Q
	efp_to_montgomery(&mapped_P,&mapped_P);
	efp2_to_montgomery(&mapped_Q,&mapped_Q);

    Pseudo_8_sparse_mapping_montgomery(&mapped_P,&mapped_Q,&L);
    efp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    efp2_set(&T,&mapped_Q);     //set T
    //TODO:1->RmodP?
    fp_set_ui(&f.x0.x0.x0,1);

    //miller
    for(i=bls12_X_length-1; i>=0; i--){
        switch(bls12_X_binary[i]){
            case 0:
                ff_ltt_lazy_montgomery(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt_lazy_montgomery(&f,&T,&mapped_P,&L);
                f_ltq_lazy_montgomery(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt_lazy_montgomery(&f,&T,&mapped_P,&L);
                f_ltq_lazy_montgomery(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
    }

    fp12_mod_montgomery(ANS,&f);
}

// void bls12_miller_algo_for_optate_projective(fp12_t *ANS,efp12_t *P,efp12_t *Q){
//     efp2_t tmp_Q;
//     efp2_init(&tmp_Q);
//     efp2_projective_t T;
//     efp2_t T2;
//     efp2_projective_init(&T);
//     efp2_projective_t mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
//     efp2_projective_init(&mapped_Q);
//     efp2_projective_init(&mapped_Q_neg);
//     efp2_projective_init(&mapped_Q1);
//     efp2_projective_init(&mapped_Q2_neg);
//     efp_t mapped_P;
//     efp_init(&mapped_P);
//     fp12_t f;
//     fp12_init(&f);
//     fp_t L;
//     fp_init(&L);
//     int i;

//     //set
//     efp12_to_efp(&mapped_P,P);//set mapped_P
//     efp12_to_efp2(&tmp_Q,Q);//set mapped_Q

// 	efp_to_montgomery(&mapped_P,&mapped_P);
// 	efp2_to_montgomery(&tmp_Q,&tmp_Q);

//     efp2_affine_to_projective_montgomery(&mapped_Q, &tmp_Q);
//     efp2_projective_set(&mapped_Q_neg, &mapped_Q); //set mapped_Q_neg
//     fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
//     efp2_projective_set(&T, &mapped_Q); //set T
//     //fp_set_ui(&f.x0.x0.x0,1);
//     fp_set_mpn(&f.x0.x0.x0,RmodP);
//     //miller
//     for(i=bls12_X_length-1; i>=0; i--){
//         switch(bls12_X_binary[i]){
//             case 0:
//                 ff_ltt_projective_lazy_montgomery(&f,&T,&mapped_P);
//                 break;
//             case 1:
//                 ff_ltt_projective_lazy_montgomery(&f,&T,&mapped_P);
//                 f_ltq_projective_lazy_montgomery(&f,&T,&mapped_Q,&mapped_P);
//                 break;
//             case -1:
//                 ff_ltt_projective_lazy_montgomery(&f,&T,&mapped_P);
//                 f_ltq_projective_lazy_montgomery(&f,&T,&mapped_Q_neg,&mapped_P);
//                 break;
//             default:
//                 break;
//         }
//     }

//     fp12_mod_montgomery(ANS,&f);
// }
void bls12_miller_algo_for_optate_projective(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    efp2_t tmp_Q;
    efp2_init(&tmp_Q);
    efp2_projective_t T;
    //efp2_t T2;
    efp2_projective_init(&T);
    efp2_projective_t mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    efp2_projective_init(&mapped_Q);
    efp2_projective_init(&mapped_Q_neg);
    efp2_projective_init(&mapped_Q1);
    efp2_projective_init(&mapped_Q2_neg);
    efp_t mapped_P_ltt,mapped_P_ltq;
    efp_init(&mapped_P_ltt);
    efp_init(&mapped_P_ltq);
    fp12_t f;
    fp12_init(&f);
    // fp_t L;
    // fp_init(&L);
    int i;

    //set
    efp12_to_efp(&mapped_P_ltq,P);//set mapped_P
    efp12_to_efp2(&tmp_Q,Q);//set mapped_Q

	efp_to_montgomery(&mapped_P_ltq,&mapped_P_ltq);
	efp2_to_montgomery(&tmp_Q,&tmp_Q);

    efp2_affine_to_projective_montgomery(&mapped_Q, &tmp_Q);
    efp2_projective_set(&mapped_Q_neg, &mapped_Q); //set mapped_Q_neg
    fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    efp2_projective_set(&T, &mapped_Q); //set T
    //fp_set_ui(&f.x0.x0.x0,1);
    fp_set_mpn(&f.x0.x0.x0,RmodP);

    //precompute P
    fp_add(&mapped_P_ltt.x,&mapped_P_ltq.x,&mapped_P_ltq.x);
    fp_add(&mapped_P_ltt.x,&mapped_P_ltt.x,&mapped_P_ltq.x);
    fp_set_neg(&mapped_P_ltt.y,&mapped_P_ltq.y);
    fp_set_neg(&mapped_P_ltq.x,&mapped_P_ltq.x);

    //miller
    for(i=bls12_X_length-1; i>=0; i--){
        switch(bls12_X_binary[i]){
            case 0:
                ff_ltt_projective_lazy_montgomery(&f,&T,&mapped_P_ltt);
                break;
            case 1:
                ff_ltt_projective_lazy_montgomery(&f,&T,&mapped_P_ltt);
                f_ltq_projective_lazy_montgomery(&f,&T,&mapped_Q,&mapped_P_ltq);
                break;
            case -1:
                ff_ltt_projective_lazy_montgomery(&f,&T,&mapped_P_ltt);
                f_ltq_projective_lazy_montgomery(&f,&T,&mapped_Q_neg,&mapped_P_ltq);
                break;
            default:
                break;
        }
    }

    fp12_mod_montgomery(ANS,&f);
}
