#include <ELiPS/bls12_miller.h>
//bls12
void BLS12_Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp12 Test;
	EFp12_init(&Test);
	EFp2 T;
	EFp2_init(&T);
	EFp2 mapped_Q;
	EFp2_init(&mapped_Q);
	EFp mapped_P;
	EFp_init(&mapped_P);
	Fp12 f;
	Fp12_init(&f);
	Fp L;
	Fp_init(&L);
	mpz_t loop;
	mpz_init(loop);
	mpz_sub_ui(loop,trace_z,1);
	int i,length;
	length=(int)mpz_sizeinbase(loop,2);
	char binary[length];
	mpz_get_str(binary,2,loop);
	
	EFp12_to_EFp(&mapped_P,P);
	EFp12_to_EFp2(&mapped_Q,Q);
	Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
	EFp2_set(&T,&mapped_Q);     //set T
	Fp_set_ui(&f.x0.x0.x0,1);   //set f
	
	//miller
    for(i=1; binary[i]!='\0'; i++){
		ff_ltt(&f,&T,&mapped_P,&L);
		if(binary[i]=='1'){
			f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
		}
	}
	
	Fp12_set(ANS,&f);
	
	
	mpz_clear(loop);
}

void BLS12_Miller_algo_for_opt_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2_neg);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    //set
    EFp12_to_EFp(&mapped_P,P);//set mapped_P
    EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);     //set T
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=BLS12_X_length-1; i>=0; i--){
        switch(BLS12_X_binary[i]){
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
    
    
    Fp12_set(ANS,&f);
}
void BLS12_Miller_algo_for_opt_ate_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2_neg);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    //set
    EFp12_to_EFp(&mapped_P,P);//set mapped_P
    EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);     //set T
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=BLS12_X_length-1; i>=0; i--){
        switch(BLS12_X_binary[i]){
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
    
    Fp12_set(ANS,&f);
}
void BLS12_Miller_algo_for_x_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    //set
    EFp12_to_EFp(&mapped_P,P);//set mapped_P
    EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);     //set T
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=BLS12_X_length-1; i>=0; i--){
        switch(BLS12_X_binary[i]){
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
    
    Fp12_frobenius_map_p3(&Buf.x,&f);       //f←f*f^(p^3)
    Fp12_mul(&f,&Buf.x,&f);
    EFp2_skew_frobenius_map_p3(&mapped_Q1,&T);//mapped_Q1←T^(p^3)
    EFp2_set(&mapped_Q2,&T);
    f_ltq(&f,&mapped_Q2,&mapped_Q1,&mapped_P,&L);
    Fp12_frobenius_map_p10(&Buf.x,&f);      //f←f*f^(p^10)
    Fp12_mul(&f,&Buf.x,&f);
    EFp2_skew_frobenius_map_p10(&T,&mapped_Q2);//T←Q2^(p^10)
    f_ltq(&f,&T,&mapped_Q2,&mapped_P,&L);
    
    Fp12_set(ANS,&f);
}
