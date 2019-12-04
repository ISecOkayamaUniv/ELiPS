#include <ELiPS/fp.h>

void fp_init(fp_t *A){
    mpn_zero(A->x0,FPLIMB);
}

void fp_printf(char *str,fp_t *A){
    gmp_printf("%s%Nu",str,A->x0,FPLIMB);
}
void fp_println(char *str,fp_t *A){
    gmp_printf("%s%Nu\n",str,A->x0,FPLIMB);
}
void fp_printf_montgomery(char *str,fp_t *A){
    static fp_t out;
    fp_mod_montgomery(&out,A);
    gmp_printf("%s%Nu",str,out.x0,FPLIMB);
}
void fp_set(fp_t *ANS,fp_t *A){
    mpn_copyd(ANS->x0,A->x0,FPLIMB);
}

void fp_set_ui(fp_t *ANS,unsigned long int UI){
    mpn_set_ui(ANS->x0,FPLIMB,UI);
}
void fp_set_mpn(fp_t *ANS,mp_limb_t *A){
    mpn_copyd(ANS->x0,A,FPLIMB);
}

void fp_set_neg(fp_t *ANS,fp_t *A){
    mpn_sub_n(ANS->x0,prime,A->x0,FPLIMB);
    fp_mod(ANS,ANS->x0,FPLIMB);
}
void fp_lshift(fp_t *ANS,fp_t *A,unsigned long int UI){
    mpn_lshift(ANS->x0,A->x0,FPLIMB,UI);
    fp_mod(ANS,ANS->x0,FPLIMB);
}
void fp_lshift2(fp_t *ANS,fp_t *A){
    static mp_limb_t buf[FPLIMB];    
    mpn_lshift(buf,A->x0,FPLIMB,1);
	if(mpn_cmp(buf,prime,FPLIMB)>=0)mpn_sub_n(ANS->x0,buf,prime,FPLIMB);
}
void fp_hlv(fp_t *ANS,fp_t *A)
{
	static mp_limb_t buf[FPLIMB];
	if(A->x0[0] & 1) mpn_add_n(buf,A->x0,prime,FPLIMB);
	else mpn_copyd(buf,A->x0,FPLIMB);
	mpn_rshift(ANS->x0,buf,FPLIMB,1);
}
void fp_set_random(fp_t *ANS,gmp_randstate_t state){
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_urandomm(tmp,state,prime_z);
    mpn_set_mpz(ANS->x0,tmp);
    
    //mpn_random(buf,FPLIMB);
    //fp_mod(ANS,buf,FPLIMB);
    
    mpz_clear(tmp);
}

void pre_montgomery(){
	mp_limb_t tmp1[FPLIMB+1],tmp2[FPLIMB2+2];
    mpz_t N2;
    mpz_t NN;
    mpz_init(N2);
    mpz_init(NN);
    mpz_ui_pow_ui(NN,2,FPLIMB_BITS);
    mpz_invert(N2,prime_z,NN);
    mpz_sub(N2,NN,N2);
    mpn_set_mpz(u,N2);
    mpn_zero(N,FPLIMB2);
    for(int i=0;i<FPLIMB;i++){
        N[i]=prime[i];
    }
    mpn_set_mpz(tmp1,NN);
    mpn_mod(RmodP,tmp1,FPLIMB+1);
    mpn_mul_n(tmp2,tmp1,tmp1,FPLIMB+1);
    mpn_mod(R2,tmp2,FPLIMB2+2);
    mpn_mul(tmp2,R2,FPLIMB,tmp1,FPLIMB+1);
    mpn_mod(R3,tmp2,FPLIMB2+2);
    mpz_clear(N2);
    mpz_clear(NN);
}
void fp_mulmod_montgomery(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
	unsigned long int carry;
	mp_limb_t r;
    static mp_limb_t T[FPLIMB2];
    static unsigned int index=0; 
	static unsigned long int c;
    static int i,j;
    mpn_zero(T,FPLIMB2);
    
    mpn_mul_n(T,A->x0,B->x0,FPLIMB);
    index=0;
    for (i = 0; i < FPLIMB; i++,index++) {
		r = (mp_limb_t)(T[index] * u[0]);
		T[index] = mpn_addmul_1(T+index,prime,FPLIMB,r);
	}
	carry = mpn_add_n(ANS->x0, T+FPLIMB, T, FPLIMB);
	if (carry || (mpn_cmp(ANS->x0, prime, FPLIMB) != -1)) {
		carry = mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
	}
}
void mpn_mulmod_montgomery(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
    
	unsigned long int carry;
	mp_limb_t r;
    static mp_limb_t T[FPLIMB2];
    static unsigned int index=0; 
	static unsigned long int c;
    static int i,j;
    mpn_zero(T,FPLIMB2);
    
    mpn_mul(T,A,A_size,B,B_size);
    index=0;
    for (i = 0; i < FPLIMB; i++,index++) {
		r = (mp_limb_t)(T[index] * u[0]);
		T[index] = mpn_addmul_1(T+index,prime,FPLIMB,r);
	}
	carry = mpn_add_n(ANS, T+FPLIMB, T, FPLIMB);
	if (carry || (mpn_cmp(ANS, prime, FPLIMB) != -1)) {
		carry = mpn_sub_n(ANS,ANS,prime,FPLIMB);
	}
}
void fp_mod_montgomery(fp_t *ANS,fp_t *A){
	unsigned long int carry;
	mp_limb_t r;
    static mp_limb_t T[FPLIMB2];
    static unsigned int index=0; 
	static unsigned long int c;
    static int i,j;
    mpn_zero(T,FPLIMB2);
    
    mpn_copyd(T,A->x0,FPLIMB);
    index=0;
    for (i = 0; i < FPLIMB; i++,index++) {
		r = (mp_limb_t)(T[index] * u[0]);
		T[index] = mpn_addmul_1(T+index,prime,FPLIMB,r);
	}
	carry = mpn_add_n(ANS->x0, T+FPLIMB, T, FPLIMB);
	if (carry || (mpn_cmp(ANS->x0, prime, FPLIMB) != -1)) {
		carry = mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
	}
}
void mpn_mod_montgomery(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size){
    
	unsigned long int carry;
	mp_limb_t r;
    static mp_limb_t T[FPLIMB2];
    static unsigned int index=0; 
	static unsigned long int c;
    static int i,j;
    mpn_zero(T,FPLIMB2);
    
    mpn_copyd(T,A,A_size);
    index=0;
    for (i = 0; i < FPLIMB; i++,index++) {
		r = (mp_limb_t)(T[index] * u[0]);
		T[index] = mpn_addmul_1(T+index,prime,FPLIMB,r);
	}
	carry = mpn_add_n(ANS, T+FPLIMB, T, FPLIMB);
	if (carry || (mpn_cmp(ANS, prime, FPLIMB) != -1)) {
		carry = mpn_sub_n(ANS,ANS,prime,FPLIMB);
	}
}
void fp_to_montgomery(fp_t *ANS, fp_t *A){
	static int i;
	static mp_limb_t tmp[FPLIMB2];
    mpn_zero(tmp,FPLIMB2);
    for(i=FPLIMB;i<FPLIMB2;i++)tmp[i]=A->x0[i-FPLIMB];
    mpn_mod(ANS->x0,tmp,FPLIMB2);
}
void mpn_to_montgomery(mp_limb_t *ANS, mp_limb_t *A){
	static int i;
	static mp_limb_t tmp[FPLIMB2];
    mpn_zero(tmp,FPLIMB2);
    for(i=FPLIMB;i<FPLIMB2;i++)tmp[i]=A[i-FPLIMB];
    mpn_mod(ANS,tmp,FPLIMB2);
}

void fp_mod(fp_t *ans,mp_limb_t *a,mp_size_t size_a){
	mp_limb_t dumy[size_a];
	mpn_tdiv_qr(dumy,ans->x0,0,a,size_a,prime,FPLIMB);
}

void fp_mod_ui(fp_t *ans,mp_limb_t *a,mp_size_t size_a,unsigned long int UI){
	mp_limb_t dumy[size_a];

	mpn_set_ui(buf,FPLIMB,UI);
	mpn_tdiv_qr(dumy,ans->x0,0,a,size_a,buf,1);
}
void fp_mul(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
    static mp_limb_t tmp_mul[FPLIMB2];
    if(mpn_zero_p(A->x0,FPLIMB)==1||mpn_zero_p(B->x0,FPLIMB)==1){
	    mpn_set_ui(ANS->x0,FPLIMB,0);
    }else if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
        if(mpn_cmp(B->x0,prime,FPLIMB)>0)fp_mod(ANS,B->x0,FPLIMB);
	    else fp_set(ANS,B);
    }else if(mpn_cmp_ui(B->x0,FPLIMB,1)==0){
        if(mpn_cmp(A->x0,prime,FPLIMB)>0)fp_mod(ANS,A->x0,FPLIMB);
	    else fp_set(ANS,A);
    }else{
        mpn_mul_n(tmp_mul,A->x0,B->x0,FPLIMB);
        fp_mod(ANS,tmp_mul,FPLIMB2);
    }
}
void fp_mul_lazy(mp_limb_t *ANS,mp_limb_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
    if(mpn_zero_p(A,FPLIMB)==1||mpn_zero_p(B,FPLIMB)==1){
	mpn_set_ui(ANS,FPLIMB2,0);
    }else if(mpn_cmp_ui(A,FPLIMB,1)==0){
	mpn_zero(ANS,FPLIMB2);
	mpn_copyd(ANS,B,FPLIMB);
    }else if(mpn_cmp_ui(B,FPLIMB,1)==0){
	mpn_zero(ANS,FPLIMB2);
	mpn_copyd(ANS,A,FPLIMB);
    }else{
        mpn_mul_n(ANS,A,B,FPLIMB);
    }
}
void fp_mul_final(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
	static mp_limb_t tmp_mul[FPLIMB2];
    if(mpn_zero_p(A->x0,FPLIMB)==1||mpn_zero_p(B->x0,FPLIMB)==1){
	mpn_set_ui(ANS->x0,FPLIMB,0);
    }else if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
	fp_mod(ANS,B->x0,FPLIMB);
    }else if(mpn_cmp_ui(B->x0,FPLIMB,1)==0){
	fp_mod(ANS,A->x0,FPLIMB);
    }else{
        mpn_mul_n(tmp_mul,A->x0,B->x0,FPLIMB);
        fp_mod(ANS,tmp_mul,FPLIMB2);
    }
}
void fp_mul_ui(fp_t *ANS,fp_t *A,unsigned long int UI){
    #ifdef DEBUG_COST_A
    cost_mul_ui++;
    #endif
    static mp_limb_t tmp_mul[FPLIMB2];
    mpn_mul_ui(tmp_mul,A->x0,FPLIMB,UI);
    fp_mod(ANS,tmp_mul,FPLIMB2);
}
void fp_mul_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
    static mp_limb_t tmp_mul[FPLIMB2];
    if(mpn_zero_p(A->x0,FPLIMB)==1||mpn_zero_p(B,FPLIMB)==1){
	mpn_set_ui(ANS->x0,FPLIMB,0);
    }else if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
	mpn_copyd(ANS->x0,B,FPLIMB);
    }else if(mpn_cmp_ui(B,FPLIMB,1)==0){
        fp_set(ANS,A);
    }else{
        mpn_mul_n(tmp_mul,A->x0,B,FPLIMB);
        fp_mod(ANS,tmp_mul,FPLIMB2);
    }
}
void fp_sqr(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_sqr++;
    #endif
    static mp_limb_t tmp_sqr[FPLIMB2];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	mpn_set_ui(ANS->x0,FPLIMB,0);
    }else if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
	mpn_set_ui(ANS->x0,FPLIMB,1);
    }else{
        mpn_sqr(tmp_sqr,A->x0,FPLIMB);
        fp_mod(ANS,tmp_sqr,FPLIMB2);
    }
}
void fp_sqr_lazy(mp_limb_t *ANS,mp_limb_t *A){
    #ifdef DEBUG_COST_A
    cost_sqr++;
    #endif
    if(mpn_zero_p(A,FPLIMB)==1){
	mpn_set_ui(ANS,FPLIMB2,0);
    }else if(mpn_cmp_ui(A,FPLIMB,1)==0){
	mpn_set_ui(ANS,FPLIMB2,1);
    }else{
        mpn_sqr(ANS,A,FPLIMB);
    }
}
void fp_add(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	fp_set(ANS,B);
    }else if(mpn_zero_p(B->x0,FPLIMB)==1){
	fp_set(ANS,A);
    }else{
	mpn_add_n(ANS->x0,A->x0,B->x0,FPLIMB);
	if(mpn_cmp(ANS->x0,prime,FPLIMB)>=0)mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
    }
}
void fp_add_lazy(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    if(mpn_zero_p(A,A_size)==1){
	mpn_copyd(ANS,B,B_size);
    }else if(mpn_zero_p(B,B_size)==1){
	mpn_copyd(ANS,A,A_size);
    }else{
	mpn_add_n(ANS,A,B,ANS_size);
    }
}
void fp_add_lazy_mod(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    if(mpn_zero_p(A,A_size)==1){
	mpn_copyd(ANS,B,B_size);
    }else if(mpn_zero_p(B,B_size)==1){
	mpn_copyd(ANS,A,A_size);
    }else{
	mpn_add_n(ANS,A,B,ANS_size);
	//TODO ANS<-buf
	if(mpn_cmp(ANS,prime,FPLIMB)>=0)mpn_sub_n(ANS,ANS,prime,FPLIMB);
    }
}
void fp_add_final(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	fp_mod(ANS,B->x0,FPLIMB);
    }else if(mpn_zero_p(B->x0,FPLIMB)==1){
	fp_mod(ANS,A->x0,FPLIMB);
    }else{
	mpn_add_n(buf,A->x0,B->x0,FPLIMB);
	fp_mod(ANS,buf,FPLIMB);
    }

}
void fp_add_ui(fp_t *ANS,fp_t *A,unsigned long int UI){
    #ifdef DEBUG_COST_A
    cost_add_ui++;
    #endif
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	fp_set_ui(ANS,UI);
    }else if(UI==0){
	fp_set(ANS,A);
    }else{
	mpn_add_ui(buf,A->x0,FPLIMB,UI);
	if(mpn_cmp(buf,prime,FPLIMB)>0)mpn_sub_n(ANS->x0,buf,prime,FPLIMB);
	else mpn_copyd(ANS->x0,buf,FPLIMB);
    }
}
void fp_add_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	mpn_copyd(ANS->x0,B,FPLIMB);
    }else if(mpn_zero_p(B,FPLIMB)==1){
	fp_set(ANS,A);
    }else{
	mpn_add_n(buf,A->x0,B,FPLIMB);
	if(mpn_cmp(buf,prime,FPLIMB)>0)mpn_sub_n(ANS->x0,buf,prime,FPLIMB);
	else mpn_copyd(ANS->x0,buf,FPLIMB);
    }
}
void fp_sub(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_sub++;
    #endif
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(B->x0,FPLIMB)==1){
	    fp_set(ANS,A);
    }else{
        if(mpn_cmp(A->x0,B->x0,FPLIMB)<0){
    	    mpn_sub_n(buf,A->x0,B->x0,FPLIMB);
    	    mpn_add_n(ANS->x0,prime,buf,FPLIMB);
        }else{
	    mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB);
	}
    }
}
void fp_sub_lazy(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size){
    #ifdef DEBUG_COST_A
    cost_sub++;
    #endif
    static mp_limb_t buf[FPLIMB],bufL[FPLIMB2];
    //assert(A_size>=B_size);
    //Only:A_size>=B_size
    
    if(mpn_zero_p(B,B_size)==1){
		//mpn_zero(ANS,ANS_size);
		mpn_copyd(ANS,A,A_size);
    }else if(mpn_zero_p(A,A_size)==1){
	    mpn_zero(ANS,ANS_size);
	    mpn_mod(buf,B,B_size);
    	mpn_sub_n(ANS,prime,buf,FPLIMB);
    }else{
		if(A_size==FPLIMB2 && B_size==FPLIMB){
	    	if(mpn_chk_limb(A,FPLIMB,FPLIMB2)==0){
	    	    mpn_sub(ANS,A,A_size,B,B_size);
        	}else{
				if(mpn_cmp(A,B,FPLIMB)<0){
			    	mpn_sub_n(ANS,A,B,FPLIMB);
    				while(mpn_cmp(ANS,prime,FPLIMB)>=0)mpn_add_n(ANS,ANS,prime,FPLIMB);
    				
    				//TODO:If there bug, shoud use this.
			    	//mpn_sub_n(ANS,B,A,FPLIMB);
    				//while(mpn_cmp(ANS,prime,FPLIMB)>=0)mpn_sub_n(ANS,ANS,prime,FPLIMB);
    				//mpn_sub_n(ANS,prime,ANS,FPLIMB);
        		}else{
        			//TODO:bug
			    	mpn_zero(ANS,ANS_size);
			    	mpn_sub_n(ANS,A,B,FPLIMB);
	    	    }
	    	}
		}else if(A_size==FPLIMB2 && B_size==FPLIMB2){
		    if(mpn_cmp(A,B,ANS_size)<0){
    			mpn_sub_n(ANS,A,B,ANS_size);
				mpn_add_n(ANS + FPLIMB, ANS + FPLIMB,prime, FPLIMB);
    	    }else{
				mpn_sub_n(ANS,A,B,ANS_size);
		    }
		}else if(A_size==FPLIMB && B_size==FPLIMB){
			if(ANS_size==FPLIMB2){
				printf("error:limb size");
				exit(1);
			}
		    if(mpn_cmp(A,B,FPLIMB)<0){
    			mpn_sub_n(ANS,A,B,FPLIMB);
    			while(mpn_cmp(ANS,prime,FPLIMB)>=0)mpn_add_n(ANS,ANS,prime,FPLIMB);
    			
    			//TODO:If there bug, shoud use this.
    			//mpn_sub_n(ANS,B,A,FPLIMB);
    			//while(mpn_cmp(ANS,prime,FPLIMB)>=0)mpn_sub_n(ANS,ANS,prime,FPLIMB);
    			//mpn_sub_n(ANS,prime,ANS,FPLIMB);
    	    }else{
				mpn_sub_n(ANS,A,B,ANS_size);
		    }
		}else{
			printf("error:limb size");
			exit(1);
		}
    }
}
void fp_sub_lazy_mod(fp_t *ANS,mp_limb_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_sub++;
    #endif
    //Only:A_size=B_size=FPLIMB2
    static mp_limb_t buf[FPLIMB2];
    if(mpn_zero_p(B,FPLIMB2)==1){
	fp_mod(ANS,A,FPLIMB2);
    }else{
	if(mpn_cmp(A,B,FPLIMB2)<0){
    	    mpn_sub_n(buf,B,A,FPLIMB2);
	    	fp_mod(ANS,buf,FPLIMB2);
    	    mpn_sub_n(ANS->x0,prime,ANS->x0,FPLIMB);
        }else{
	    mpn_sub_n(buf,A,B,FPLIMB2);
	    fp_mod(ANS,buf,FPLIMB2);
	}
    }
}
void fp_sub_final(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_sub++;
    #endif
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(B->x0,FPLIMB)==1){
	fp_mod(ANS,A->x0,FPLIMB);
    }else{
        if(mpn_cmp(A->x0,B->x0,FPLIMB)<0){
    	    mpn_sub_n(buf,B->x0,A->x0,FPLIMB);
	    fp_mod(ANS,buf,FPLIMB);
    	    mpn_sub_n(ANS->x0,prime,ANS->x0,FPLIMB);
        }else{
	    mpn_sub_n(buf,A->x0,B->x0,FPLIMB);
	    fp_mod(ANS,buf,FPLIMB);
	}
    }
}
void fp_sub_ui(fp_t *ANS,fp_t *A,unsigned long int UI){
    #ifdef DEBUG_COST_A
    cost_sub_ui++;
    #endif
	static mp_limb_t buf[FPLIMB];
    if(UI==0){
	fp_set(ANS,A);
    }else{
	mpn_sub_ui(ANS->x0,A->x0,FPLIMB,UI);
    }
}
void fp_sub_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_sub++;
    #endif
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(B,FPLIMB)==1){
	fp_set(ANS,A);
    }else{
        if(mpn_cmp(A->x0,B,FPLIMB)<0){
    	    mpn_sub_n(buf,B,A->x0,FPLIMB);
    	    mpn_sub_n(buf,prime,buf,FPLIMB);
	    fp_mod(ANS,buf,FPLIMB);
        }else{
	    mpn_sub_n(buf,A->x0,B,FPLIMB);
	    fp_mod(ANS,buf,FPLIMB);
	}
    }
}
void fp_inv(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_inv++;
    #endif
    static mp_limb_t prime_tmp[FPLIMB],gp[FPLIMB],sp[FPLIMB],tmp[FPLIMB],buf[FPLIMB];
    static mp_size_t buf_size;

    //mpn_init(gp,FPLIMB);
    mpn_init(sp,FPLIMB);
    //mpn_init(tmp,FPLIMB);
    //mpn_init(prime_tmp,FPLIMB);

    mpn_copyd(prime_tmp,prime,FPLIMB);
	
    mpn_add_n(buf,A->x0,prime,FPLIMB);
    mpn_gcdext(gp,sp,&buf_size,buf,FPLIMB,prime_tmp,FPLIMB);
	
    if( buf_size < 0 ){
        mpn_sub_n(tmp, prime,sp,FPLIMB);
    }else{	
	    mpn_copyd(tmp,sp,FPLIMB);
    }
    
    fp_mod(ANS,tmp,FPLIMB);
}
void fp_inv_montgomery(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
	fp_inv(ANS,A);	
	mpn_mulmod_montgomery(ANS->x0,FPLIMB,ANS->x0,FPLIMB,R3,FPLIMB);
	//gmp_printf("R2=%Nu\n",R2,FPLIMB);
}

int  fp_legendre(fp_t *A){

    int i;
    mpz_t tmp1,tmp2;
    fp_t tmp1_fp;
    mpz_init(tmp1);
    mpz_init(tmp2);
    fp_init(&tmp1_fp);
    
    mpz_sub_ui(tmp1,prime_z,1);
    mpz_tdiv_q_ui(tmp2,tmp1,2);
    fp_pow(&tmp1_fp,A,tmp2);
	
    if(mpn_cmp_ui(tmp1_fp.x0,FPLIMB,1)==0)		i=1;
    else if(mpn_cmp_ui(tmp1_fp.x0,FPLIMB,0)==0)	i=0;
    else					i=-1;
    
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    
    return i;
}

int  fp_isCNR(fp_t *A){
    fp_t tmp;
    fp_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_sub_ui(exp,prime_z,1);
    mpz_tdiv_q_ui(exp,exp,3);
    fp_pow(&tmp,A,exp);
    
    if(fp_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}
void fp_sqrt(fp_t *ANS,fp_t *A){
    fp_t x,y,t,k,n,tmp;
    fp_init(&x);
    fp_init(&y);
    fp_init(&t);
    fp_init(&k);
    fp_init(&n);
    fp_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state1;
    gmp_randinit_default (state1);
    gmp_randseed_ui(state1,(unsigned long)time(NULL));
    fp_set_random(&n,state1);
    
    while(fp_legendre(&n)!=-1){
        fp_set_random(&n,state1);
    }
    mpz_sub_ui(q,prime_z,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    fp_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp_pow(&x,A,exp);
    fp_mul(&tmp,&x,&x);
    fp_mul(&k,&tmp,A);
    fp_mul(&x,&x,A);
    while(mpn_cmp_ui(k.x0,FPLIMB,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fp_pow(&tmp,&k,exp);
        while(mpn_cmp_ui(tmp.x0,FPLIMB,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fp_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fp_pow(&t,&y,result);
        fp_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fp_mul(&x,&x,&t);
        fp_mul(&k,&k,&y);
    }
    fp_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}
void fp_pow(fp_t *ANS,fp_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp_t tmp;
    fp_init(&tmp);
    
    fp_set(&tmp,A);
    
    for(i=1;i<length; i++){
        fp_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            fp_mul(&tmp,A,&tmp);
        }
    }
    fp_set(ANS,&tmp);
    
}
void fp_pow_mpn(fp_t *ans,fp_t *a,mp_limb_t *r,mp_size_t n){
    fp_t Temp,tmp;
    mp_limb_t bit[FPLIMB],bit_copy[FPLIMB];
    size_t bit_size;
    mp_size_t size;
    int i,cnt;
    size=FPLIMB;
	
    fp_init(&Temp);
    fp_init(&tmp);
    mpn_init(bit,size);
    mpn_init(bit_copy,size);
	
    if(mpn_zero_p(r,size)==1){
    	fp_set_ui(ans,1);
    	return;
}
    bit_size=mpn_sizeinbase(r,n,2);
    bit_size--;
	
    mpn_set_ui(bit,1,size);
    mpn_lshift_ext(bit,bit,size,bit_size);
	
    //scm
    mpn_copyd(Temp.x0,a->x0,size);
    while(bit_size>0){
	mpn_copyd(tmp.x0,Temp.x0,size);
	fp_mul(&Temp,&tmp,&tmp);
	//bit
	mpn_rshift(bit,bit,size,1);
	mpn_and_n(bit_copy,r,bit,size);
	//bit=1 -> 
	if(mpn_zero_p(bit_copy,size)==0){
		mpn_copyd(tmp.x0,Temp.x0,size);
		fp_mul(&Temp,&tmp,a);
	}
	mpn_init(bit_copy,size);
	bit_size--;
	
    }
    mpn_copyd(ans->x0,Temp.x0,size);
}

int  fp_cmp(fp_t *A,fp_t *B){
    if(mpn_cmp(A->x0,B->x0,FPLIMB)==0){
        return 0;   
    }
    return 1;
}

int  fp_cmp_ui(fp_t *A,unsigned long int UI){
    if(mpn_cmp_ui(A->x0,FPLIMB,UI)==0){
        return 0;
    }
    return 1;
}

int  fp_cmp_mpn(fp_t *A,mp_limb_t *B){
    if(mpn_cmp(A->x0,B,FPLIMB)==0){
        return 0;
    }
    return 1;
}

int  fp_cmp_zero(fp_t *A){
    if(mpn_cmp_ui(A->x0,FPLIMB,0)==0){
        return 0;
    }
    return 1;
}

int  fp_cmp_one(fp_t *A){
    if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
        return 0;
    }
    return 1;
}

int fp_montgomery_trick(fp_t *A_inv,fp_t *A,int n){
    int i;
    fp_t ANS[n],ALL_inv;
	fp_set(&ANS[0],&A[0]);
	fp_t check;
	
	for(i=1;i<n;i++){
	fp_mul(&ANS[i],&ANS[i-1],&A[i]);
	}
	fp_inv(&ALL_inv,&ANS[n-1]);	
	for(i=n-1;i>0;i--){
    fp_mul(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    fp_mul(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    fp_set(&A_inv[0],&ALL_inv);
    /*
    for(i=0;i<n;i++){
    fp_mul(&check,&A[i],&A_inv[i]);
    printf("check:%d",i);	
	fp_println("=",&check);
    }
    */
    return 0;
}
int fp_montgomery_trick_montgomery(fp_t *A_inv,fp_t *A,int n){
    int i;
    fp_t ANS[n],ALL_inv;
	fp_set(&ANS[0],&A[0]);
	fp_t check;
	
	for(i=1;i<n;i++){
	fp_mulmod_montgomery(&ANS[i],&ANS[i-1],&A[i]);
	}
	fp_inv_montgomery(&ALL_inv,&ANS[n-1]);
	for(i=n-1;i>0;i--){
    fp_mulmod_montgomery(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    fp_mulmod_montgomery(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    fp_set(&A_inv[0],&ALL_inv);
    /*
    for(i=0;i<n;i++){
    fp_mul(&check,&A[i],&A_inv[i]);
    printf("check:%d",i);	
	fp_println("=",&check);
    }
    */
    return 0;
}