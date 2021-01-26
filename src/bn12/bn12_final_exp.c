#include <ELiPS/bn12_final_exp.h>
/*----------------------------------------------------------------------------*/
//bn12
void bn12_final_exp_plain(fp12_t *ANS,fp12_t *A){
    fp12_t Tmp,Buf1,Buf2;
	fp12_init(&Tmp);
	fp12_set(&Tmp,A);
	fp12_init(&Buf1);
	fp12_init(&Buf2);
	mpz_t exp,buf;
	mpz_init(exp);
	mpz_init(buf);
	
	fp12_frobenius_map_p6(&Buf1,&Tmp);
	fp12_inv(&Buf2,&Tmp);
	fp12_mul(&Tmp,&Buf1,&Buf2);
	
	fp12_frobenius_map_p2(&Buf1,&Tmp);
	fp12_mul(&Tmp,&Buf1,&Tmp);
	
	mpz_pow_ui(exp,prime_z,4);
	mpz_pow_ui(buf,prime_z,2);
	mpz_sub(exp,exp,buf);
	mpz_add_ui(exp,exp,1);
	mpz_tdiv_q(exp,exp,order_z);
	fp12_pow(ANS,&Tmp,exp);
	
	mpz_clear(exp);
	mpz_clear(buf);
}


void bn12_fp12_pow_X(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
	fp12_init(&tmp);
	fp12_init(&A_inv);
	fp12_frobenius_map_p6(&A_inv,A);
	
	fp12_set(&tmp,A);
	for(i=bn12_X_length-1; i>=0; i--){
		switch(bn12_X_binary[i]){
			case 0:
				fp12_sqr_cyclotomic(&tmp,&tmp);
				break;
			case 1:
				fp12_sqr_cyclotomic(&tmp,&tmp);
				fp12_mul(&tmp,&tmp,A);
				break;
			case -1:
				fp12_sqr_cyclotomic(&tmp,&tmp);
				fp12_mul(&tmp,&tmp,&A_inv);
				break;
			default:
				break;
		}
	}
	fp12_set(ANS,&tmp);
}
void bn12_fp12_pow_X_lazy(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
	fp12_init(&tmp);
	fp12_init(&A_inv);
	fp12_frobenius_map_p6(&A_inv,A);
	
	fp12_set(&tmp,A);
	for(i=bn12_X_length-1; i>=0; i--){
		switch(bn12_X_binary[i]){
			case 0:
				fp12_sqr_cyclotomic_lazy(&tmp,&tmp);
				break;
			case 1:
				fp12_sqr_cyclotomic_lazy(&tmp,&tmp);
				fp12_mul_lazy(&tmp,&tmp,A);
				break;
			case -1:
				fp12_sqr_cyclotomic_lazy(&tmp,&tmp);
				fp12_mul_lazy(&tmp,&tmp,&A_inv);
				break;
			default:
				break;
		}
	}
	fp12_set(ANS,&tmp);
}
void bn12_final_exp_optimal(fp12_t *ANS,fp12_t *A){
    fp12_t t0,t1,t2,t3,t4;
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t4);
    
    gettimeofday(&tv_start,NULL);
    //f←f^(p^6)*f^-1
    fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    fp12_inv(&t1,A);//f^-1
    fp12_mul(A,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    fp12_frobenius_map_p2(&t0,A);//f^(p^2)
    fp12_mul(A,&t0,A);//f^(p^2)*f
    
    
    bn12_fp12_pow_X(&t0,A);   //t0←f^(-u)
    fp12_frobenius_map_p6(&t0,&t0);
    fp12_sqr_cyclotomic(&t0,&t0);              //t0←t0^2
    fp12_sqr_cyclotomic(&t1,&t0);              //t1←t0^2
    fp12_mul(&t1,&t0,&t1);              //t1←t0*t1
    bn12_fp12_pow_X(&t2,&t1);        //t2←t1^(-u)
    fp12_frobenius_map_p6(&t2,&t2);
    fp12_frobenius_map_p6(&t3,&t1);         //t3←t1^-1
    fp12_mul(&t1,&t2,&t3);              //t1←t2*t3
    fp12_sqr_cyclotomic(&t3,&t2);              //t3←t2^2
    bn12_fp12_pow_X(&t4,&t3);        //t4←t3^(-u)
    fp12_frobenius_map_p6(&t4,&t4);
    fp12_frobenius_map_p6(&t4,&t4);         //t4←t4^(-1)
    fp12_mul(&t4,&t4,&t1);              //t4←t4*t1
    fp12_mul(&t3,&t4,&t0);              //t3←t4*t0
    fp12_mul(&t0,&t2,&t4);              //t0←t2*t4
    fp12_mul(&t0,&t0,A);             //t0←t0*f
    fp12_frobenius_map_p1(&t2,&t3);         //t2←t3^p
    fp12_mul(&t0,&t2,&t0);              //t0←t2*t0
    fp12_frobenius_map_p2(&t2,&t4);         //t2←t4^(p^2)
    fp12_mul(&t0,&t2,&t0);              //t0←t2*t0
    fp12_frobenius_map_p6(&t2,A);            //t2←f^(-1)
    fp12_mul(&t2,&t2,&t3);              //t2←t2*t3
    fp12_frobenius_map_p3(&t2,&t2);         //t2←t2^(p^3)
    fp12_mul(ANS,&t2,&t0);              //t0←t2*t0
    
}
void bn12_final_exp_optimal_lazy(fp12_t *ANS,fp12_t *A){
    fp12_t t0,t1,t2,t3,t4;
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t4);
    
    gettimeofday(&tv_start,NULL);
    //f←f^(p^6)*f^-1
    fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    fp12_inv(&t1,A);//f^-1
    fp12_mul_lazy(A,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    fp12_frobenius_map_p2(&t0,A);//f^(p^2)
    fp12_mul(A,&t0,A);//f^(p^2)*f
    
    
    bn12_fp12_pow_X_lazy(&t0,A);   //t0←f^(-u)
    fp12_frobenius_map_p6(&t0,&t0);
    fp12_sqr_cyclotomic_lazy(&t0,&t0);              //t0←t0^2
    fp12_sqr_cyclotomic_lazy(&t1,&t0);              //t1←t0^2
    fp12_mul_lazy(&t1,&t0,&t1);              //t1←t0*t1
    bn12_fp12_pow_X_lazy(&t2,&t1);        //t2←t1^(-u)
    fp12_frobenius_map_p6(&t2,&t2);
    fp12_frobenius_map_p6(&t3,&t1);         //t3←t1^-1
    fp12_mul_lazy(&t1,&t2,&t3);              //t1←t2*t3
    fp12_sqr_cyclotomic_lazy(&t3,&t2);              //t3←t2^2
    bn12_fp12_pow_X_lazy(&t4,&t3);        //t4←t3^(-u)
    fp12_frobenius_map_p6(&t4,&t4);
    fp12_frobenius_map_p6(&t4,&t4);         //t4←t4^(-1)
    fp12_mul_lazy(&t4,&t4,&t1);              //t4←t4*t1
    fp12_mul_lazy(&t3,&t4,&t0);              //t3←t4*t0
    fp12_mul_lazy(&t0,&t2,&t4);              //t0←t2*t4
    fp12_mul_lazy(&t0,&t0,A);             //t0←t0*f
    fp12_frobenius_map_p1(&t2,&t3);         //t2←t3^p
    fp12_mul_lazy(&t0,&t2,&t0);              //t0←t2*t0
    fp12_frobenius_map_p2(&t2,&t4);         //t2←t4^(p^2)
    fp12_mul_lazy(&t0,&t2,&t0);              //t0←t2*t0
    fp12_frobenius_map_p6(&t2,A);            //t2←f^(-1)
    fp12_mul_lazy(&t2,&t2,&t3);              //t2←t2*t3
    fp12_frobenius_map_p3(&t2,&t2);         //t2←t2^(p^3)
    fp12_mul_lazy(ANS,&t2,&t0);              //t0←t2*t0
    
}
