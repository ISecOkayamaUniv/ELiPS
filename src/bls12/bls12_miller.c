#include <ELiPS/bls12_miller.h>
//bls12
void bls12_miller_algo_for_plain_ate(fp12_t *ANS,efp12_t *P,efp12_t *Q){
    efp12_t Test;
	efp12_init(&Test);
	efp2_t T;
	efp2_init(&T);
	efp2_t mapped_Q;
	efp2_init(&mapped_Q);
	efp_t mapped_P;
	efp_init(&mapped_P);
	fp12_t f;
	fp12_init(&f);
	fp_t L;
	fp_init(&L);
	mpz_t loop;
	mpz_init(loop);
	mpz_sub_ui(loop,trace_z,1);
	int i,length;
	length=(int)mpz_sizeinbase(loop,2);
	char binary[length];
	mpz_get_str(binary,2,loop);
	
	efp12_to_efp(&mapped_P,P);
	efp12_to_efp2(&mapped_Q,Q);
	Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
	efp2_set(&T,&mapped_Q);     //set T
	fp_set_ui(&f.x0.x0.x0,1);   //set f
	
	//miller
    for(i=1; binary[i]!='\0'; i++){
		ff_ltt(&f,&T,&mapped_P,&L);
		if(binary[i]=='1'){
			f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
		}
	}
	
	fp12_set(ANS,&f);
	
	
	mpz_clear(loop);
}

void bls12_miller_algo_for_opt_ate(fp12_t *ANS,efp12_t *P,efp12_t *Q){
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
void bls12_miller_algo_for_opt_ate_lazy(fp12_t *ANS,efp12_t *P,efp12_t *Q){
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
                ff_ltt_lazy(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt_lazy(&f,&T,&mapped_P,&L);
                f_ltq_lazy(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt_lazy(&f,&T,&mapped_P,&L);
                f_ltq_lazy(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
    }
    
    fp12_set(ANS,&f);
}
void bls12_miller_algo_for_opt_ate_lazy_montgomery(fp12_t *ANS,efp12_t *P,efp12_t *Q){
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
