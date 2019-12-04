#include <ELiPS/bls12_final_exp.h>

/*----------------------------------------------------------------------------*/
//bls12
void bls12_final_exp_plain(fp12_t *ANS,fp12_t *A){
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

void bls12_fp12_pow_X(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
	fp12_init(&tmp);
	fp12_init(&A_inv);
	fp12_frobenius_map_p6(&A_inv,A);
	
	fp12_set(&tmp,A);
	for(i=bls12_X_length-1; i>=0; i--){
		switch(bls12_X_binary[i]){
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
void bls12_fp12_pow_X_compress(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_t tmp6,tmp9,tmp11,tmp77;
    
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6(&A_inv,A);
    fp12_set(&tmp,A);
	
    for(i=0;i<77;i++){
        if(i<=10)	fp12_sqr_GS(&tmp,&tmp);
	else		fp12_sqr_compressed(&tmp,&tmp);
	if(i==5){
	fp12_set(&tmp6,&tmp);
	}
	if(i==8){
	fp12_set(&tmp9,&tmp);
	}
	if(i==10){
	fp12_set(&tmp11,&tmp);
	}
    }
	fp12_sqr_recover_g1(&tmp77,&tmp);
	fp12_sqr_recover_g0(&tmp77,&tmp77);
    
    fp12_mul(&tmp,&tmp77,&tmp11);
    fp12_inv(&tmp9,&tmp9);
    fp12_inv(&tmp6,&tmp6);
    fp12_mul(&tmp,&tmp,&tmp9);
    fp12_mul(&tmp,&tmp,&tmp6);
    fp12_set(ANS,&tmp);
}
//TODO;fast monty trick
void bls12_fp12_pow_X_compress2(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t buf,tmp,A_inv;
    fp12_t tmp6,tmp9,tmp11,tmp77;
    
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6(&A_inv,A);
    fp12_set(&tmp,A);
	
    for(i=0;i<77;i++){
    fp12_sqr_compressed(&tmp,&tmp);
	if(i==5){
	fp12_set(&tmp6,&tmp);
	}
	if(i==8){
	fp12_set(&tmp9,&tmp);
	}
	if(i==10){
	fp12_set(&tmp11,&tmp);
	}
    }
	fp12_sqr_recover_g1_noninv(&tmp6,&tmp6);
	fp12_sqr_recover_g1_noninv(&tmp9,&tmp9);
	fp12_sqr_recover_g1_noninv(&tmp11,&tmp11);
	fp12_sqr_recover_g1_noninv(&tmp77,&tmp);
	
    fp12_sqr_recover_g1_montrick(&tmp6,&tmp9,&tmp11,&tmp77);
	
	fp12_sqr_recover_g0(&tmp6,&tmp6);
	fp12_sqr_recover_g0(&tmp9,&tmp9);
	fp12_sqr_recover_g0(&tmp11,&tmp11);
	fp12_sqr_recover_g0(&tmp77,&tmp77);
    
    fp12_mul(&tmp,&tmp77,&tmp11);
    fp12_mul(&buf,&tmp6,&tmp9);
    fp12_frobenius_map_p6(&buf,&buf);
    fp12_mul(&tmp,&tmp,&buf);
    fp12_set(ANS,&tmp);
}

void bls12_fp12_pow_X_lazy(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
	fp12_init(&tmp);
	fp12_init(&A_inv);
	fp12_frobenius_map_p6(&A_inv,A);
	
	fp12_set(&tmp,A);
	for(i=bls12_X_length-1; i>=0; i--){
		switch(bls12_X_binary[i]){
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
void bls12_fp12_pow_X_compress_lazy(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_t tmp6,tmp9,tmp11,tmp77;
    
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6(&A_inv,A);
    fp12_set(&tmp,A);
	
    for(i=0;i<77;i++){
        if(i<=10)	fp12_sqr_GS_lazy(&tmp,&tmp);
	else		fp12_sqr_compressed_lazy(&tmp,&tmp);
	if(i==5){
	fp12_set(&tmp6,&tmp);
	}
	if(i==8){
	fp12_set(&tmp9,&tmp);
	}
	if(i==10){
	fp12_set(&tmp11,&tmp);
	}
    }
	fp12_sqr_recover_g1_lazy(&tmp77,&tmp);
	fp12_sqr_recover_g0_lazy(&tmp77,&tmp77);
    
    fp12_mul_lazy(&tmp,&tmp77,&tmp11);
    fp12_frobenius_map_p6(&tmp9,&tmp9);
    fp12_frobenius_map_p6(&tmp6,&tmp6);
    fp12_mul_lazy(&tmp,&tmp,&tmp9);
    fp12_mul_lazy(&tmp,&tmp,&tmp6);
    fp12_set(ANS,&tmp);
}
void bls12_fp12_pow_X_compress_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_t tmp6,tmp9,tmp11,tmp77;
    
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6_montgomery(&A_inv,A);
    fp12_set(&tmp,A);
	
    for(i=0;i<77;i++){
        if(i<=10)	fp12_sqr_GS_lazy_montgomery(&tmp,&tmp);
		else		fp12_sqr_compressed_lazy_montgomery(&tmp,&tmp);
		if(i==5){
		fp12_set(&tmp6,&tmp);
		}
		if(i==8){
		fp12_set(&tmp9,&tmp);
		}
		if(i==10){
		fp12_set(&tmp11,&tmp);
		}
    }
    
	fp12_sqr_recover_g1_lazy_montgomery(&tmp77,&tmp);
	fp12_sqr_recover_g0_lazy_montgomery(&tmp77,&tmp77);
    
    fp12_mul_lazy_montgomery(&tmp,&tmp77,&tmp11);
    fp12_frobenius_map_p6_montgomery(&tmp9,&tmp9);
    fp12_frobenius_map_p6_montgomery(&tmp6,&tmp6);
    fp12_mul_lazy_montgomery(&tmp,&tmp,&tmp9);
    fp12_mul_lazy_montgomery(&tmp,&tmp,&tmp6);
    fp12_set(ANS,&tmp);
}
void bls12_fp12_pow_X2(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6(&A_inv,A);
    
    fp12_set(&tmp,A);
    for(i=bls12_X2_length-1; i>=0; i--){
        switch(bls12_X2_binary[i]){
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

void bls12_fp12_pow_X2_compress(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_t tmp5,tmp8,tmp10,tmp76;
    
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6(&A_inv,A);
    fp12_set(&tmp,A);
	
    for(i=0;i<76;i++){
        if(i<=10)	fp12_sqr_GS(&tmp,&tmp);
	else		fp12_sqr_compressed(&tmp,&tmp);
	if(i==4){
	fp12_set(&tmp5,&tmp);
	}
	if(i==7){
	fp12_set(&tmp8,&tmp);
	}
	if(i==9){
	fp12_set(&tmp10,&tmp);
	}
    }
	fp12_sqr_recover_g1(&tmp76,&tmp);
	fp12_sqr_recover_g0(&tmp76,&tmp76);
    
    fp12_mul(&tmp,&tmp76,&tmp10);
    fp12_inv(&tmp8,&tmp8);
    fp12_inv(&tmp5,&tmp5);
    fp12_mul(&tmp,&tmp,&tmp8);
    fp12_mul(&tmp,&tmp,&tmp5);
    fp12_set(ANS,&tmp);
}
void bls12_fp12_pow_X2_lazy(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6(&A_inv,A);
    
    fp12_set(&tmp,A);
    for(i=bls12_X2_length-1; i>=0; i--){
        switch(bls12_X2_binary[i]){
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
void bls12_fp12_pow_X2_compress_lazy(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_t tmp5,tmp8,tmp10,tmp76;
    
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6(&A_inv,A);
    fp12_set(&tmp,A);
	
    for(i=0;i<76;i++){
        if(i<=10)	fp12_sqr_GS_lazy(&tmp,&tmp);
	else		fp12_sqr_compressed_lazy(&tmp,&tmp);
	if(i==4){
	fp12_set(&tmp5,&tmp);
	}
	if(i==7){
	fp12_set(&tmp8,&tmp);
	}
	if(i==9){
	fp12_set(&tmp10,&tmp);
	}
    }
	fp12_sqr_recover_g1_lazy(&tmp76,&tmp);
	fp12_sqr_recover_g0_lazy(&tmp76,&tmp76);
    
    fp12_mul_lazy(&tmp,&tmp76,&tmp10);
    fp12_frobenius_map_p6(&tmp8,&tmp8);
    fp12_frobenius_map_p6(&tmp5,&tmp5);
    fp12_mul_lazy(&tmp,&tmp,&tmp8);
    fp12_mul_lazy(&tmp,&tmp,&tmp5);
    fp12_set(ANS,&tmp);
}
void bls12_fp12_pow_X2_compress_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    int i;
    fp12_t tmp,A_inv;
    fp12_t tmp5,tmp8,tmp10,tmp76;
    
    fp12_init(&tmp);
    fp12_init(&A_inv);
    fp12_frobenius_map_p6_montgomery(&A_inv,A);
    fp12_set(&tmp,A);
	
    for(i=0;i<76;i++){
        if(i<=10)	fp12_sqr_GS_lazy_montgomery(&tmp,&tmp);
		else	fp12_sqr_compressed_lazy_montgomery(&tmp,&tmp);
		if(i==4){
		fp12_set(&tmp5,&tmp);
		}
		if(i==7){
		fp12_set(&tmp8,&tmp);
		}
		if(i==9){
		fp12_set(&tmp10,&tmp);
		}
    }
    
	fp12_sqr_recover_g1_lazy_montgomery(&tmp76,&tmp);
	fp12_sqr_recover_g0_lazy_montgomery(&tmp76,&tmp76);
    
    fp12_mul_lazy_montgomery(&tmp,&tmp76,&tmp10);
    fp12_frobenius_map_p6_montgomery(&tmp8,&tmp8);
    fp12_frobenius_map_p6_montgomery(&tmp5,&tmp5);
    fp12_mul_lazy_montgomery(&tmp,&tmp,&tmp8);
    fp12_mul_lazy_montgomery(&tmp,&tmp,&tmp5);
    fp12_set(ANS,&tmp);
}
void bls12_final_exp_optimal(fp12_t *ANS,fp12_t *A){
    fp12_t tmp,t0,t1,t2,t3,t4,t5, test;
    fp12_init(&tmp);
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t5);
    fp12_init(&t4);
    fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    fp12_inv(&t1,A);//f^-1
    fp12_mul(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    fp12_mul(&tmp,&t0,&tmp);//f^(p^2)*f
    
    //HARD PART
    fp12_sqr_cyclotomic(&t0, &tmp);
    bls12_fp12_pow_X(&t1, &t0);
    
    bls12_fp12_pow_X2(&t2,&t1);//t2:=t1^(u2);
    fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    fp12_mul(&t1,&t3,&t1);//t1:=t3*t1;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    bls12_fp12_pow_X(&t2,&t1);//t2:=t1^(u);
    bls12_fp12_pow_X(&t3,&t2);//t3:=t2^(u);
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    fp12_mul(&t3,&t1,&t3);//t3:=t1*t3;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_frobenius_map_p3(&t1,&t1);//t1:=t1^(p^3);
    fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    bls12_fp12_pow_X(&t2,&t3);//t2:=t3^(u);
    fp12_mul(&t2,&t2,&t0);//t2:=t2*t0;
    fp12_mul(&t2,&t2,&tmp);//t2:=t2*f;
    fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    fp12_frobenius_map_p1(&t2,&t3);//t2:=t3^p;
    fp12_mul(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
}

void bls12_final_exp_optimal_compress(fp12_t *ANS,fp12_t *A){
    fp12_t tmp,t0,t1,t2,t3,t4,t5, test;
    fp12_init(&tmp);
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t5);
    fp12_init(&t4);
    fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    fp12_inv(&t1,A);//f^-1
    fp12_mul(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    fp12_mul(&tmp,&t0,&tmp);//f^(p^2)*f
    
    //HARD PART
    fp12_sqr_cyclotomic(&t0, &tmp);
    bls12_fp12_pow_X_compress2(&t1, &t0);
    
    bls12_fp12_pow_X2_compress(&t2,&t1);//t2:=t1^(u2);
    fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    fp12_mul(&t1,&t3,&t1);//t1:=t3*t1;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    bls12_fp12_pow_X_compress(&t2,&t1);//t2:=t1^(u);
    bls12_fp12_pow_X_compress(&t3,&t2);//t3:=t2^(u);
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    fp12_mul(&t3,&t1,&t3);//t3:=t1*t3;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_frobenius_map_p3(&t1,&t1);//t1:=t1^(p^3);
    fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    bls12_fp12_pow_X_compress(&t2,&t3);//t2:=t3^(u);
    fp12_mul(&t2,&t2,&t0);//t2:=t2*t0;
    fp12_mul(&t2,&t2,&tmp);//t2:=t2*f;
    fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    fp12_frobenius_map_p1(&t2,&t3);//t2:=t3^p;
    fp12_mul(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
}
void bls12_final_exp_optimal_lazy(fp12_t *ANS,fp12_t *A){
    fp12_t tmp,t0,t1,t2,t3,t4,t5, test;
    fp12_init(&tmp);
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t5);
    fp12_init(&t4);
    fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    fp12_inv_lazy(&t1,A);//f^-1
    fp12_mul_lazy(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    fp12_mul_lazy(&tmp,&t0,&tmp);//f^(p^2)*f
    
    
    //HARDPART
    fp12_sqr_cyclotomic_lazy(&t0, &tmp);
    bls12_fp12_pow_X_lazy(&t1, &t0);
    
    bls12_fp12_pow_X2_lazy(&t2,&t1);//t2:=t1^(u2);
    fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    fp12_mul_lazy(&t1,&t3,&t1);//t1:=t3*t1;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    bls12_fp12_pow_X_lazy(&t2,&t1);//t2:=t1^(u);
    bls12_fp12_pow_X_lazy(&t3,&t2);//t3:=t2^(u);
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    fp12_mul_lazy(&t3,&t1,&t3);//t3:=t1*t3;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_frobenius_map_p3_lazy(&t1,&t1);//t1:=t1^(p^3);
    fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    bls12_fp12_pow_X_lazy(&t2,&t3);//t2:=t3^(u);
    fp12_mul_lazy(&t2,&t2,&t0);//t2:=t2*t0;
    fp12_mul_lazy(&t2,&t2,&tmp);//t2:=t2*f;
    fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    fp12_frobenius_map_p1_lazy(&t2,&t3);//t2:=t3^p;
    fp12_mul_lazy(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
    
}
void bls12_final_exp_optimal_compress_lazy(fp12_t *ANS,fp12_t *A){
    fp12_t tmp,t0,t1,t2,t3,t4,t5, test;
    fp12_init(&tmp);
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t5);
    fp12_init(&t4);
    fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    fp12_inv_lazy(&t1,A);//f^-1
    fp12_mul_lazy(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    fp12_frobenius_map_p2(&t0,&tmp);//f^(p^2)
    fp12_mul_lazy(&tmp,&t0,&tmp);//f^(p^2)*f
    
    
    //HARDPART
    fp12_sqr_GS_lazy(&t0, &tmp);
    bls12_fp12_pow_X_compress_lazy(&t1, &t0);
    
    bls12_fp12_pow_X2_compress_lazy(&t2,&t1);//t2:=t1^(u2);
    fp12_frobenius_map_p6(&t3,&tmp);//t3:=f^(-1);
    fp12_mul_lazy(&t1,&t3,&t1);//t1:=t3*t1;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    bls12_fp12_pow_X_compress_lazy(&t2,&t1);//t2:=t1^(u);
    bls12_fp12_pow_X_compress_lazy(&t3,&t2);//t3:=t2^(u);
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    
    fp12_mul_lazy(&t3,&t1,&t3);//t3:=t1*t3;
    fp12_frobenius_map_p6(&t1,&t1);//t1:=t1^(-1);
    fp12_frobenius_map_p3_lazy(&t1,&t1);//t1:=t1^(p^3);
    fp12_frobenius_map_p2(&t2,&t2);//t2:=t2^(p^2);
    
    
    fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    bls12_fp12_pow_X_compress_lazy(&t2,&t3);//t2:=t3^(u);
    fp12_mul_lazy(&t2,&t2,&t0);//t2:=t2*t0;
    fp12_mul_lazy(&t2,&t2,&tmp);//t2:=t2*f;
    fp12_mul_lazy(&t1,&t1,&t2);//t1:=t1*t2;
    
    fp12_frobenius_map_p1_lazy(&t2,&t3);//t2:=t3^p;
    fp12_mul_lazy(ANS,&t1,&t2);//t1:=t1*t2;

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
    
}
void bls12_final_exp_optimal_compress_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    fp12_t tmp,t0,t1,t2,t3,t4,t5, test,At;
    fp12_init(&tmp);
    fp12_init(&t0);
    fp12_init(&t1);
    fp12_init(&t2);
    fp12_init(&t3);
    fp12_init(&t5);
    fp12_init(&t4);
    fp12_init(&test);
    fp12_init(&At);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    
    fp12_to_montgomery(&At,A);
    
    //EASY PART
    mpz_set(positive_X,X_z);
    mpz_tdiv_q_ui(positive_X2,positive_X,2);
    
    //f←f^(p^6)*f^-1
    fp12_frobenius_map_p6_montgomery(&t0,&At);//f^(p^6)
    fp12_inv_lazy_montgomery(&t1,&At);//f^-1
    fp12_mul_lazy_montgomery(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    fp12_frobenius_map_p2_montgomery(&t0,&tmp);//f^(p^2)
    fp12_mul_lazy_montgomery(&tmp,&t0,&tmp);//f^(p^2)*f
    
    
    //HARDPART
    fp12_sqr_GS_lazy_montgomery(&t0, &tmp);
    bls12_fp12_pow_X_compress_lazy_montgomery(&t1, &t0);
    
    bls12_fp12_pow_X2_compress_lazy_montgomery(&t2,&t1);//t2:=t1^(u2);
    fp12_frobenius_map_p6_montgomery(&t3,&tmp);//t3:=f^(-1);
    fp12_mul_lazy_montgomery(&t1,&t3,&t1);//t1:=t3*t1;
    fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);
    fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;
    
    bls12_fp12_pow_X_compress_lazy_montgomery(&t2,&t1);//t2:=t1^(u);
    bls12_fp12_pow_X_compress_lazy_montgomery(&t3,&t2);//t3:=t2^(u);
    fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);
    
    fp12_mul_lazy_montgomery(&t3,&t1,&t3);//t3:=t1*t3;
    fp12_frobenius_map_p6_montgomery(&t1,&t1);//t1:=t1^(-1);
    fp12_frobenius_map_p3_lazy_montgomery(&t1,&t1);//t1:=t1^(p^3);
    fp12_frobenius_map_p2_montgomery(&t2,&t2);//t2:=t2^(p^2);
    
    
    fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;
    bls12_fp12_pow_X_compress_lazy_montgomery(&t2,&t3);//t2:=t3^(u);
    fp12_mul_lazy_montgomery(&t2,&t2,&t0);//t2:=t2*t0;
    fp12_mul_lazy_montgomery(&t2,&t2,&tmp);//t2:=t2*f;
    fp12_mul_lazy_montgomery(&t1,&t1,&t2);//t1:=t1*t2;
    
    fp12_frobenius_map_p1_lazy_montgomery(&t2,&t3);//t2:=t3^p;
    fp12_mul_lazy_montgomery(ANS,&t1,&t2);//t1:=t1*t2;
    
    fp12_mod_montgomery(ANS,ANS);

    mpz_clear(positive_X);
    mpz_clear(positive_X2);
    
}