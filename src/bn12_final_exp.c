#include <ELiPS/bn12_final_exp.h>
/*----------------------------------------------------------------------------*/
//bn12
void BN12_Final_exp_plain(Fp12 *ANS,Fp12 *A){
    Fp12 Tmp,Buf1,Buf2;
	Fp12_init(&Tmp);
	Fp12_set(&Tmp,A);
	Fp12_init(&Buf1);
	Fp12_init(&Buf2);
	mpz_t exp,buf;
	mpz_init(exp);
	mpz_init(buf);
	
	Fp12_frobenius_map_p6(&Buf1,&Tmp);
	Fp12_inv(&Buf2,&Tmp);
	Fp12_mul(&Tmp,&Buf1,&Buf2);
	
	Fp12_frobenius_map_p2(&Buf1,&Tmp);
	Fp12_mul(&Tmp,&Buf1,&Tmp);
	
	mpz_pow_ui(exp,prime_z,4);
	mpz_pow_ui(buf,prime_z,2);
	mpz_sub(exp,exp,buf);
	mpz_add_ui(exp,exp,1);
	mpz_tdiv_q(exp,exp,order_z);
	Fp12_pow(ANS,&Tmp,exp);
	
	mpz_clear(exp);
	mpz_clear(buf);
}


void BN12_Fp12_pow_X(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
	Fp12_init(&tmp);
	Fp12_init(&A_inv);
	Fp12_frobenius_map_p6(&A_inv,A);
	
	Fp12_set(&tmp,A);
	for(i=BN12_X_length-1; i>=0; i--){
		switch(BN12_X_binary[i]){
			case 0:
				Fp12_sqr_cyclotomic(&tmp,&tmp);
				break;
			case 1:
				Fp12_sqr_cyclotomic(&tmp,&tmp);
				Fp12_mul(&tmp,&tmp,A);
				break;
			case -1:
				Fp12_sqr_cyclotomic(&tmp,&tmp);
				Fp12_mul(&tmp,&tmp,&A_inv);
				break;
			default:
				break;
		}
	}
	Fp12_set(ANS,&tmp);
}
void BN12_Fp12_pow_X_lazy(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
	Fp12_init(&tmp);
	Fp12_init(&A_inv);
	Fp12_frobenius_map_p6(&A_inv,A);
	
	Fp12_set(&tmp,A);
	for(i=BN12_X_length-1; i>=0; i--){
		switch(BN12_X_binary[i]){
			case 0:
				Fp12_sqr_cyclotomic_lazy(&tmp,&tmp);
				break;
			case 1:
				Fp12_sqr_cyclotomic_lazy(&tmp,&tmp);
				Fp12_mul_lazy(&tmp,&tmp,A);
				break;
			case -1:
				Fp12_sqr_cyclotomic_lazy(&tmp,&tmp);
				Fp12_mul_lazy(&tmp,&tmp,&A_inv);
				break;
			default:
				break;
		}
	}
	Fp12_set(ANS,&tmp);
}
void BN12_Final_exp_optimal(Fp12 *ANS,Fp12 *A){
    Fp12 t0,t1,t2,t3,t4;
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t4);
    
    gettimeofday(&tv_start,NULL);
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv(&t1,A);//f^-1
    Fp12_mul(A,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,A);//f^(p^2)
    Fp12_mul(A,&t0,A);//f^(p^2)*f
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_EASY=timedifference_msec(tv_start,tv_end);
    
    gettimeofday(&tv_start,NULL);
    
    BN12_Fp12_pow_X(&t0,A);   //t0←f^(-u)
    Fp12_frobenius_map_p6(&t0,&t0);
    Fp12_sqr_cyclotomic(&t0,&t0);              //t0←t0^2
    Fp12_sqr_cyclotomic(&t1,&t0);              //t1←t0^2
    Fp12_mul(&t1,&t0,&t1);              //t1←t0*t1
    BN12_Fp12_pow_X(&t2,&t1);        //t2←t1^(-u)
    Fp12_frobenius_map_p6(&t2,&t2);
    Fp12_frobenius_map_p6(&t3,&t1);         //t3←t1^-1
    Fp12_mul(&t1,&t2,&t3);              //t1←t2*t3
    Fp12_sqr_cyclotomic(&t3,&t2);              //t3←t2^2
    BN12_Fp12_pow_X(&t4,&t3);        //t4←t3^(-u)
    Fp12_frobenius_map_p6(&t4,&t4);
    Fp12_frobenius_map_p6(&t4,&t4);         //t4←t4^(-1)
    Fp12_mul(&t4,&t4,&t1);              //t4←t4*t1
    Fp12_mul(&t3,&t4,&t0);              //t3←t4*t0
    Fp12_mul(&t0,&t2,&t4);              //t0←t2*t4
    Fp12_mul(&t0,&t0,A);             //t0←t0*f
    Fp12_frobenius_map_p1(&t2,&t3);         //t2←t3^p
    Fp12_mul(&t0,&t2,&t0);              //t0←t2*t0
    Fp12_frobenius_map_p2(&t2,&t4);         //t2←t4^(p^2)
    Fp12_mul(&t0,&t2,&t0);              //t0←t2*t0
    Fp12_frobenius_map_p6(&t2,A);            //t2←f^(-1)
    Fp12_mul(&t2,&t2,&t3);              //t2←t2*t3
    Fp12_frobenius_map_p3(&t2,&t2);         //t2←t2^(p^3)
    Fp12_mul(ANS,&t2,&t0);              //t0←t2*t0
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_HARD=timedifference_msec(tv_start,tv_end);
}
void BN12_Final_exp_optimal_lazy(Fp12 *ANS,Fp12 *A){
    Fp12 t0,t1,t2,t3,t4;
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t4);
    
    gettimeofday(&tv_start,NULL);
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv(&t1,A);//f^-1
    Fp12_mul_lazy(A,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,A);//f^(p^2)
    Fp12_mul(A,&t0,A);//f^(p^2)*f
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_EASY=timedifference_msec(tv_start,tv_end);
    
    gettimeofday(&tv_start,NULL);
    
    BN12_Fp12_pow_X_lazy(&t0,A);   //t0←f^(-u)
    Fp12_frobenius_map_p6(&t0,&t0);
    Fp12_sqr_cyclotomic_lazy(&t0,&t0);              //t0←t0^2
    Fp12_sqr_cyclotomic_lazy(&t1,&t0);              //t1←t0^2
    Fp12_mul_lazy(&t1,&t0,&t1);              //t1←t0*t1
    BN12_Fp12_pow_X_lazy(&t2,&t1);        //t2←t1^(-u)
    Fp12_frobenius_map_p6(&t2,&t2);
    Fp12_frobenius_map_p6(&t3,&t1);         //t3←t1^-1
    Fp12_mul_lazy(&t1,&t2,&t3);              //t1←t2*t3
    Fp12_sqr_cyclotomic_lazy(&t3,&t2);              //t3←t2^2
    BN12_Fp12_pow_X_lazy(&t4,&t3);        //t4←t3^(-u)
    Fp12_frobenius_map_p6(&t4,&t4);
    Fp12_frobenius_map_p6(&t4,&t4);         //t4←t4^(-1)
    Fp12_mul_lazy(&t4,&t4,&t1);              //t4←t4*t1
    Fp12_mul_lazy(&t3,&t4,&t0);              //t3←t4*t0
    Fp12_mul_lazy(&t0,&t2,&t4);              //t0←t2*t4
    Fp12_mul_lazy(&t0,&t0,A);             //t0←t0*f
    Fp12_frobenius_map_p1(&t2,&t3);         //t2←t3^p
    Fp12_mul_lazy(&t0,&t2,&t0);              //t0←t2*t0
    Fp12_frobenius_map_p2(&t2,&t4);         //t2←t4^(p^2)
    Fp12_mul_lazy(&t0,&t2,&t0);              //t0←t2*t0
    Fp12_frobenius_map_p6(&t2,A);            //t2←f^(-1)
    Fp12_mul_lazy(&t2,&t2,&t3);              //t2←t2*t3
    Fp12_frobenius_map_p3(&t2,&t2);         //t2←t2^(p^3)
    Fp12_mul_lazy(ANS,&t2,&t0);              //t0←t2*t0
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_HARD=timedifference_msec(tv_start,tv_end);
}
