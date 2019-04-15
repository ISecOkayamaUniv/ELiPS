#include <ELiPS/bls12_final_exp.h>

/*----------------------------------------------------------------------------*/
//bls12
void BLS12_Final_exp_plain(Fp12 *ANS,Fp12 *A){
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

void BLS12_Fp12_pow_X(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
	Fp12_init(&tmp);
	Fp12_init(&A_inv);
	Fp12_frobenius_map_p6(&A_inv,A);
	
	Fp12_set(&tmp,A);
	for(i=BLS12_X_length-1; i>=0; i--){
		switch(BLS12_X_binary[i]){
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
void BLS12_Fp12_pow_X_compress(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
    Fp12 tmp6,tmp9,tmp11,tmp77;
    
    Fp12_init(&tmp);
    Fp12_init(&A_inv);
    Fp12_frobenius_map_p6(&A_inv,A);
    Fp12_set(&tmp,A);
	
    for(i=0;i<77;i++){
        if(i<=10)	Fp12_sqr_GS(&tmp,&tmp);
	else		Fp12_sqr_compressed(&tmp,&tmp);
	if(i==5){
	Fp12_set(&tmp6,&tmp);
	}
	if(i==8){
	Fp12_set(&tmp9,&tmp);
	}
	if(i==10){
	Fp12_set(&tmp11,&tmp);
	}
    }
	Fp12_sqr_recover_g1(&tmp77,&tmp);
	Fp12_sqr_recover_g0(&tmp77,&tmp77);
    
    Fp12_mul(&tmp,&tmp77,&tmp11);
    Fp12_inv(&tmp9,&tmp9);
    Fp12_inv(&tmp6,&tmp6);
    Fp12_mul(&tmp,&tmp,&tmp9);
    Fp12_mul(&tmp,&tmp,&tmp6);
    Fp12_set(ANS,&tmp);
}
void BLS12_Fp12_pow_X_lazy(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
	Fp12_init(&tmp);
	Fp12_init(&A_inv);
	Fp12_frobenius_map_p6(&A_inv,A);
	
	Fp12_set(&tmp,A);
	for(i=BLS12_X_length-1; i>=0; i--){
		switch(BLS12_X_binary[i]){
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
void BLS12_Fp12_pow_X_compress_lazy(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
    Fp12 tmp6,tmp9,tmp11,tmp77;
    
    Fp12_init(&tmp);
    Fp12_init(&A_inv);
    Fp12_frobenius_map_p6(&A_inv,A);
    Fp12_set(&tmp,A);
	
    for(i=0;i<77;i++){
        if(i<=10)	Fp12_sqr_GS_lazy(&tmp,&tmp);
	else		Fp12_sqr_compressed_lazy(&tmp,&tmp);
	if(i==5){
	Fp12_set(&tmp6,&tmp);
	}
	if(i==8){
	Fp12_set(&tmp9,&tmp);
	}
	if(i==10){
	Fp12_set(&tmp11,&tmp);
	}
    }
	Fp12_sqr_recover_g1(&tmp77,&tmp);
	Fp12_sqr_recover_g0(&tmp77,&tmp77);
    
    Fp12_mul_lazy(&tmp,&tmp77,&tmp11);
    Fp12_inv_lazy(&tmp9,&tmp9);
    Fp12_inv_lazy(&tmp6,&tmp6);
    Fp12_mul_lazy(&tmp,&tmp,&tmp9);
    Fp12_mul_lazy(&tmp,&tmp,&tmp6);
    Fp12_set(ANS,&tmp);
}
void BLS12_Fp12_pow_X2(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
    Fp12_init(&tmp);
    Fp12_init(&A_inv);
    Fp12_frobenius_map_p6(&A_inv,A);
    
    Fp12_set(&tmp,A);
    for(i=BLS12_X2_length-1; i>=0; i--){
        switch(BLS12_X2_binary[i]){
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

void BLS12_Fp12_pow_X2_compress(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
    Fp12 tmp5,tmp8,tmp10,tmp76;
    
    Fp12_init(&tmp);
    Fp12_init(&A_inv);
    Fp12_frobenius_map_p6(&A_inv,A);
    Fp12_set(&tmp,A);
	
    for(i=0;i<76;i++){
        if(i<=10)	Fp12_sqr_GS(&tmp,&tmp);
	else		Fp12_sqr_compressed(&tmp,&tmp);
	if(i==4){
	Fp12_set(&tmp5,&tmp);
	}
	if(i==7){
	Fp12_set(&tmp8,&tmp);
	}
	if(i==9){
	Fp12_set(&tmp10,&tmp);
	}
    }
	Fp12_sqr_recover_g1(&tmp76,&tmp);
	Fp12_sqr_recover_g0(&tmp76,&tmp76);
    
    Fp12_mul(&tmp,&tmp76,&tmp10);
    Fp12_inv(&tmp8,&tmp8);
    Fp12_inv(&tmp5,&tmp5);
    Fp12_mul(&tmp,&tmp,&tmp8);
    Fp12_mul(&tmp,&tmp,&tmp5);
    Fp12_set(ANS,&tmp);
}
void BLS12_Fp12_pow_X2_lazy(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
    Fp12_init(&tmp);
    Fp12_init(&A_inv);
    Fp12_frobenius_map_p6(&A_inv,A);
    
    Fp12_set(&tmp,A);
    for(i=BLS12_X2_length-1; i>=0; i--){
        switch(BLS12_X2_binary[i]){
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
void BLS12_Fp12_pow_X2_compress_lazy(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
    Fp12 tmp5,tmp8,tmp10,tmp76;
    
    Fp12_init(&tmp);
    Fp12_init(&A_inv);
    Fp12_frobenius_map_p6(&A_inv,A);
    Fp12_set(&tmp,A);
	
    for(i=0;i<76;i++){
        if(i<=10)	Fp12_sqr_GS_lazy(&tmp,&tmp);
	else		Fp12_sqr_compressed_lazy(&tmp,&tmp);
	if(i==4){
	Fp12_set(&tmp5,&tmp);
	}
	if(i==7){
	Fp12_set(&tmp8,&tmp);
	}
	if(i==9){
	Fp12_set(&tmp10,&tmp);
	}
    }
	Fp12_sqr_recover_g1(&tmp76,&tmp);
	Fp12_sqr_recover_g0(&tmp76,&tmp76);
    
    Fp12_mul_lazy(&tmp,&tmp76,&tmp10);
    Fp12_inv_lazy(&tmp8,&tmp8);
    Fp12_inv_lazy(&tmp5,&tmp5);
    Fp12_mul_lazy(&tmp,&tmp,&tmp8);
    Fp12_mul_lazy(&tmp,&tmp,&tmp5);
    Fp12_set(ANS,&tmp);
}
void BLS12_Final_exp_optimal(Fp12 *ANS,Fp12 *A){
    Fp12 tmp,t0,t1,t2,t3,t4,t5, test;
    Fp12_init(&tmp);
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t5);
    Fp12_init(&t4);
    Fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv(&t1,A);//f^-1
    Fp12_mul(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    Fp12_mul(&tmp,&t0,&tmp);//f^(p^2)*f
    
    //HARD PART
    Fp12_sqr_cyclotomic(&t0, &tmp);
    BLS12_Fp12_pow_X(&t1, &t0);
    
    BLS12_Fp12_pow_X2(&t2,&t1);//t2:=t1^(u2);
    Fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    Fp12_mul(&t1,&t3,&t1);//t1:=t3*t1;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    BLS12_Fp12_pow_X(&t2,&t1);//t2:=t1^(u);
    BLS12_Fp12_pow_X(&t3,&t2);//t3:=t2^(u);
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    Fp12_mul(&t3,&t1,&t3);//t3:=t1*t3;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_frobenius_map_p3(&t1,&t1);//t1:=t1^(p^3);
    Fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    BLS12_Fp12_pow_X(&t2,&t3);//t2:=t3^(u);
    Fp12_mul(&t2,&t2,&t0);//t2:=t2*t0;
    Fp12_mul(&t2,&t2,&tmp);//t2:=t2*f;
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    Fp12_frobenius_map_p1(&t2,&t3);//t2:=t3^p;
    Fp12_mul(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
}

void BLS12_Final_exp_optimal_compress(Fp12 *ANS,Fp12 *A){
    Fp12 tmp,t0,t1,t2,t3,t4,t5, test;
    Fp12_init(&tmp);
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t5);
    Fp12_init(&t4);
    Fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv(&t1,A);//f^-1
    Fp12_mul(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    Fp12_mul(&tmp,&t0,&tmp);//f^(p^2)*f
    
    //HARD PART
    Fp12_sqr_cyclotomic(&t0, &tmp);
    BLS12_Fp12_pow_X_compress(&t1, &t0);
    
    BLS12_Fp12_pow_X2_compress(&t2,&t1);//t2:=t1^(u2);
    Fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    Fp12_mul(&t1,&t3,&t1);//t1:=t3*t1;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    BLS12_Fp12_pow_X_compress(&t2,&t1);//t2:=t1^(u);
    BLS12_Fp12_pow_X_compress(&t3,&t2);//t3:=t2^(u);
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    Fp12_mul(&t3,&t1,&t3);//t3:=t1*t3;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_frobenius_map_p3(&t1,&t1);//t1:=t1^(p^3);
    Fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    BLS12_Fp12_pow_X_compress(&t2,&t3);//t2:=t3^(u);
    Fp12_mul(&t2,&t2,&t0);//t2:=t2*t0;
    Fp12_mul(&t2,&t2,&tmp);//t2:=t2*f;
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    Fp12_frobenius_map_p1(&t2,&t3);//t2:=t3^p;
    Fp12_mul(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
}
void BLS12_Final_exp_optimal_lazy(Fp12 *ANS,Fp12 *A){
    Fp12 tmp,t0,t1,t2,t3,t4,t5, test;
    Fp12_init(&tmp);
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t5);
    Fp12_init(&t4);
    Fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv_lazy(&t1,A);//f^-1
    Fp12_mul_lazy(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    Fp12_mul_lazy(&tmp,&t0,&tmp);//f^(p^2)*f
    
    
    //HARDPART
    Fp12_sqr_cyclotomic_lazy(&t0, &tmp);
    BLS12_Fp12_pow_X_lazy(&t1, &t0);
    
    BLS12_Fp12_pow_X2_lazy(&t2,&t1);//t2:=t1^(u2);
    Fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    Fp12_mul_lazy(&t1,&t3,&t1);//t1:=t3*t1;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    BLS12_Fp12_pow_X_lazy(&t2,&t1);//t2:=t1^(u);
    BLS12_Fp12_pow_X_lazy(&t3,&t2);//t3:=t2^(u);
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    Fp12_mul_lazy(&t3,&t1,&t3);//t3:=t1*t3;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_frobenius_map_p3_lazy(&t1,&t1);//t1:=t1^(p^3);
    Fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    Fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    BLS12_Fp12_pow_X_lazy(&t2,&t3);//t2:=t3^(u);
    Fp12_mul_lazy(&t2,&t2,&t0);//t2:=t2*t0;
    Fp12_mul_lazy(&t2,&t2,&tmp);//t2:=t2*f;
    Fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    Fp12_frobenius_map_p1_lazy(&t2,&t3);//t2:=t3^p;
    Fp12_mul_lazy(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
    
}

void BLS12_Final_exp_optimal_compress_lazy(Fp12 *ANS,Fp12 *A){
    Fp12 tmp,t0,t1,t2,t3,t4,t5, test;
    Fp12_init(&tmp);
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t5);
    Fp12_init(&t4);
    Fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv_lazy(&t1,A);//f^-1
    Fp12_mul_lazy(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    Fp12_mul_lazy(&tmp,&t0,&tmp);//f^(p^2)*f
    
    
    //HARDPART
    Fp12_sqr_GS_lazy(&t0, &tmp);
    BLS12_Fp12_pow_X_compress_lazy(&t1, &t0);
    
    BLS12_Fp12_pow_X2_compress_lazy(&t2,&t1);//t2:=t1^(u2);
    Fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    Fp12_mul_lazy(&t1,&t3,&t1);//t1:=t3*t1;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    BLS12_Fp12_pow_X_compress_lazy(&t2,&t1);//t2:=t1^(u);
    BLS12_Fp12_pow_X_compress_lazy(&t3,&t2);//t3:=t2^(u);
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    Fp12_mul_lazy(&t3,&t1,&t3);//t3:=t1*t3;
    Fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    Fp12_frobenius_map_p3_lazy(&t1,&t1);//t1:=t1^(p^3);
    Fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    Fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    BLS12_Fp12_pow_X_compress_lazy(&t2,&t3);//t2:=t3^(u);
    Fp12_mul_lazy(&t2,&t2,&t0);//t2:=t2*t0;
    Fp12_mul_lazy(&t2,&t2,&tmp);//t2:=t2*f;
    Fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    Fp12_frobenius_map_p1_lazy(&t2,&t3);//t2:=t3^p;
    Fp12_mul_lazy(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
    
}
