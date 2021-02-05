#include <ELiPS/fp.h>

void fp_init(fp_t *A){
    mpn_zero(A->x0,FPLIMB);
}

void fpd_init(fpd_t *A){
    mpn_zero(A->x0,FPLIMB2);
}

void fp_printf(char *str,fp_t *A){
    gmp_printf("%s%Nu",str,A->x0,FPLIMB);
}

void fpd_printf(char *str,fpd_t *A){
    gmp_printf("%s%Nu",str,A->x0,FPLIMB2);
}

void fp_println(char *str,fp_t *A){
    gmp_printf("%s%Nu\n",str,A->x0,FPLIMB);
}

//new
void fpd_println(char *str,fpd_t *A){
    gmp_printf("%s%Nu\n",str,A->x0,FPLIMB2);
}

void fp_printf_montgomery(char *str,fp_t *A){
    static fp_t out;
    fp_mod_montgomery(&out,A);
    gmp_printf("%s%Nu",str,out.x0,FPLIMB);
}

//new
void fp_println_montgomery(char *str,fp_t *A){
    static fp_t out;
    fp_mod_montgomery(&out,A);
    gmp_printf("%s%Nu\n",str,out.x0,FPLIMB);
}

void fp_set(fp_t *ANS,fp_t *A){
    mpn_copyd(ANS->x0,A->x0,FPLIMB);
}

void fpd_set(fpd_t *ANS,fpd_t *A){
    mpn_copyd(ANS->x0,A->x0,FPLIMB2);
}

void fp_set_ui(fp_t *ANS,unsigned long int UI){
    mpn_set_ui(ANS->x0,FPLIMB,UI);
}

void fp_set_mpn(fp_t *ANS,mp_limb_t *A){
    mpn_copyd(ANS->x0,A,FPLIMB);
}

void fp_set_neg(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_ASSERT
    assert(mpn_cmp(A->x0,prime,FPLIMB)>0)
    #endif
    if(fp_cmp_zero(A)==0)
        fp_set(ANS,A);
    else
        mpn_sub_n(ANS->x0,prime,A->x0,FPLIMB);
}

void fp_set_neg_montgomery(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_ASSERT
    assert(mpn_cmp(A->x0,prime,FPLIMB)>0)
    #endif
    if(fp_cmp_zero(A)==0)
        fp_set(ANS,A);
    else
        mpn_sub_n(ANS->x0,prime,A->x0,FPLIMB);
}


void fp_lshift(fp_t *ANS,fp_t *A,unsigned long int UI){
    mpn_lshift(ANS->x0,A->x0,FPLIMB,UI);
    fp_mod(ANS,ANS->x0,FPLIMB);
}

void fp_l1shift(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    mpn_lshift(ANS->x0,A->x0,FPLIMB,1);
    if(mpn_cmp(ANS->x0,prime,FPLIMB)>=0)mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
}

void fp_r1shift(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    if(A->x0[0] & 1) mpn_add_n(ANS->x0,A->x0,prime,FPLIMB);
    else mpn_copyd(ANS->x0,A->x0,FPLIMB);
    mpn_rshift(ANS->x0,ANS->x0,FPLIMB,1);
}

void fp_set_random(fp_t *ANS,gmp_randstate_t state){
    mpz_t tmp;
    mpz_init(tmp);
    mpz_urandomm(tmp,state,prime_z);
    mpn_set_mpz(ANS->x0,tmp);
    mpz_clear(tmp);
}

void fp_set_random_montgomery(fp_t *ANS,gmp_randstate_t state){
    mpz_t tmp;
    mpz_init(tmp);
    mpz_urandomm(tmp,state,prime_z);
    mpn_set_mpz(ANS->x0,tmp);
    mpz_clear(tmp);
    fp_to_montgomery(ANS,ANS);
}

void pre_montgomery(){
    mp_limb_t tmp1[FPLIMB+1],tmp2[FPLIMB2+2];
    mpz_t tmp_z;
    mpz_t R;
    mpz_t R3_z;
    mp_limb_t R2[FPLIMB2+2];

    mpz_init(tmp_z);
    mpz_init(R);
    mpz_init(R3_z);

    for(int i=0;i<FPLIMB;i++) N[i]=prime[i];
    mpz_ui_pow_ui(R,2,FPLIMB_BITS);
    mpz_invert(tmp_z,prime_z,R);
    mpz_sub(tmp_z,R,tmp_z);
    mpn_set_mpz(tmp1,tmp_z);
    Ni_neg = tmp1[0];

    mpn_set_mpz(tmp1,R);
    mpn_mod(tmp1,tmp1,FPLIMB+1);
    mpn_copyd(RmodP,tmp1,FPLIMB);

    mpz_pow_ui(R3_z,R,3);
    mpz_mod(R3_z,R3_z,prime_z);
    mpn_set_mpz(R3,R3_z);

    mpz_clear(tmp_z);
    mpz_clear(R);
    mpz_clear(R3_z);
}

void fp_mulmod_montgomery(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    cost_mod++;
    #endif
    static mp_limb_t T[FPLIMB2];
    mpn_zero(T,FPLIMB2);

    mpn_mul_n(T,A->x0,B->x0,FPLIMB);
    for (int i = 0; i < FPLIMB; i++)
        T[i] = mpn_addmul_1(&T[i],prime,FPLIMB,T[i] * Ni_neg);

    mpn_add_n(ANS->x0, T+FPLIMB, T, FPLIMB);
    if (mpn_cmp(ANS->x0, prime, FPLIMB) != -1) mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
}

void fp_sqrmod_montgomery(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_sqr++;
    cost_mod++;
    #endif
    static mp_limb_t T[FPLIMB2];
    mpn_zero(T,FPLIMB2);

    mpn_sqr(T,A->x0,FPLIMB);
    for (int i = 0; i < FPLIMB; i++)
        T[i] = mpn_addmul_1(&T[i],prime,FPLIMB,T[i] * Ni_neg);

    mpn_add_n(ANS->x0, T+FPLIMB, T, FPLIMB);
    if (mpn_cmp(ANS->x0, prime, FPLIMB) != -1) mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
}

void fp_mod_montgomery(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_mod++;
    #endif
    static mp_limb_t T[FPLIMB2];
    mpn_zero(T,FPLIMB2);

    mpn_copyd(T,A->x0,FPLIMB);
    for (int i = 0; i < FPLIMB; i++)
        T[i] = mpn_addmul_1(&T[i],prime,FPLIMB,T[i] * Ni_neg);

    mpn_add_n(ANS->x0, T+FPLIMB, T, FPLIMB);
    if (mpn_cmp(ANS->x0, prime, FPLIMB) != -1) mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
}

void fp_to_montgomery(fp_t *ANS, fp_t *A){
    #ifdef DEBUG_COST_A
    //cost_mod++;
    cost_mod_nomal++;
    #endif
    static int i;
    static mp_limb_t tmp[FPLIMB2];
    mpn_zero(tmp,FPLIMB2);
    for(i=FPLIMB;i<FPLIMB2;i++)tmp[i]=A->x0[i-FPLIMB];
    mpn_mod(ANS->x0,tmp,FPLIMB2);
}

void fp_mod(fp_t *ans,mp_limb_t *a,mp_size_t size_a){
    #ifdef DEBUG_COST_A
    cost_mod_nomal++;
    #endif
    mp_limb_t dumy[size_a];
    mpn_tdiv_qr(dumy,ans->x0,0,a,size_a,prime,FPLIMB);
}

void fp_mod_ui(fp_t *ans,mp_limb_t *a,mp_size_t size_a,unsigned long int UI){
    mp_limb_t dumy[size_a];

    mpn_set_ui(buf,FPLIMB,UI);
    mpn_tdiv_qr(dumy,ans->x0,0,a,size_a,buf,1);
}

void fp_mul(fp_t *ANS,fp_t *A,fp_t *B){
    static mp_limb_t tmp_mul[FPLIMB2];
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif

    mpn_mul_n(tmp_mul,A->x0,B->x0,FPLIMB);
    fp_mod(ANS,tmp_mul,FPLIMB2);
}

void fp_mul_nonmod(fpd_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
    mpn_mul_n(ANS->x0,A->x0,B->x0,FPLIMB);
}

void fp_mul_ui(fp_t *ANS,fp_t *A,unsigned long int UI){
    static mp_limb_t tmp_mul[FPLIMB2];
    mpn_mul_ui(tmp_mul,A->x0,FPLIMB,UI);
    fp_mod(ANS,tmp_mul,FPLIMB2);
}
void fp_mul_ui_nonmod_single(fp_t *ANS,fp_t *A,unsigned long int UI){
    mpn_mul_ui(ANS->x0,A->x0,FPLIMB,UI);
}

void fp_mul_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_mul++;
    #endif
    static mp_limb_t tmp_mul[FPLIMB2];
    mpn_mul_n(tmp_mul,A->x0,B,FPLIMB);
    fp_mod(ANS,tmp_mul,FPLIMB2);
}

void fp_sqr(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_sqr++;
    #endif
    static mp_limb_t tmp_sqr[FPLIMB2];
    mpn_sqr(tmp_sqr,A->x0,FPLIMB);
    fp_mod(ANS,tmp_sqr,FPLIMB2);
}

void fp_sqr_nonmod(fpd_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_sqr++;
    #endif
    mpn_sqr(ANS->x0,A->x0,FPLIMB);
}

void fp_add(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    mpn_add_n(ANS->x0,A->x0,B->x0,FPLIMB);
    if(mpn_cmp(ANS->x0,prime,FPLIMB)>=0) mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
}

void fp_add_nonmod_single(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_add_nonmod++;
    #endif
    mpn_add_n(ANS->x0,A->x0,B->x0,FPLIMB);
}

void fp_add_nonmod_double(fpd_t *ANS,fpd_t *A,fpd_t *B){
    #ifdef DEBUG_COST_A
    cost_add_nonmod_double++;
    #endif
    mpn_add_n(ANS->x0,A->x0,B->x0,FPLIMB2);
}

void fp_add_ui(fp_t *ANS,fp_t *A,unsigned long int UI){
    #ifdef DEBUG_COST_A
    cost_add_ui++;
    #endif
    mpn_add_ui(ANS->x0,A->x0,FPLIMB,UI);
    if(mpn_cmp(ANS->x0,prime,FPLIMB)>0) mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
}

void fp_add_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_add++;
    #endif
    mpn_add_n(ANS->x0,A->x0,B,FPLIMB);
    if(mpn_cmp(ANS->x0,prime,FPLIMB)>0) mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
}

void fp_sub(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_sub++;
    #endif
    static mp_limb_t buf[FPLIMB];

    if(mpn_cmp(A->x0,B->x0,FPLIMB)<0){
        mpn_sub_n(buf,A->x0,B->x0,FPLIMB);
        mpn_add_n(ANS->x0,prime,buf,FPLIMB);
    }else{
        mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB);
    }
}

void fp_sub_nonmod_single(fp_t *ANS,fp_t *A,fp_t *B){
    #ifdef DEBUG_COST_A
    cost_sub_nonmod++;
    #endif

    if(mpn_cmp(A->x0,B->x0,FPLIMB)<0){
        mpn_sub_n(ANS->x0,B->x0,A->x0,FPLIMB);
        while(mpn_cmp(ANS->x0,prime,FPLIMB)>=0){
            mpn_sub_n(ANS->x0,ANS->x0,prime,FPLIMB);
        }
        mpn_sub_n(ANS->x0,prime,ANS->x0,FPLIMB);
    }else{
       mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB);
    }
}

void fp_sub_nonmod_double(fpd_t *ANS,fpd_t *A,fpd_t *B){
    #ifdef DEBUG_COST_A
    cost_sub_nonmod_double++;
    #endif
    static mp_limb_t buf[FPLIMB2];

    if(mpn_cmp(A->x0,B->x0,FPLIMB2)<0){
    		mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB2);
			mpn_add_n(ANS->x0 + FPLIMB, ANS->x0 + FPLIMB,prime, FPLIMB);
    }else{
       mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB2);
    }
}
void fp_sub_ui(fp_t *ANS,fp_t *A,unsigned long int UI){
    #ifdef DEBUG_COST_A
    cost_sub_ui++;
    #endif
    if(UI==0)   fp_set(ANS,A);
    else   mpn_sub_ui(ANS->x0,A->x0,FPLIMB,UI);
}
void fp_sub_mpn(fp_t *ANS,fp_t *A,mp_limb_t *B){
    #ifdef DEBUG_COST_A
    cost_sub++;
    #endif
    static mp_limb_t buf[FPLIMB];

    if(mpn_cmp(A->x0,B,FPLIMB)<0){
        mpn_sub_n(buf,A->x0,B,FPLIMB);
        mpn_add_n(ANS->x0,prime,buf,FPLIMB);
    }else{
        mpn_sub_n(ANS->x0,A->x0,B,FPLIMB);
    }
}
//remove fp_mod
void fp_inv(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_inv++;
    #endif
    static mp_limb_t prime_tmp[FPLIMB],gp[FPLIMB],sp[FPLIMB],buf[FPLIMB];
    static mp_size_t buf_size;

    mpn_init(sp,FPLIMB);
    mpn_copyd(prime_tmp,prime,FPLIMB);

    mpn_add_n(buf,A->x0,prime,FPLIMB);
    mpn_gcdext(gp,sp,&buf_size,buf,FPLIMB,prime_tmp,FPLIMB);

    if( buf_size < 0 ){
        mpn_sub_n(ANS->x0, prime,sp,FPLIMB);
    }else{
        mpn_copyd(ANS->x0,sp,FPLIMB);
    }
}

void fp_inv_montgomery(fp_t *ANS,fp_t *A){
    #ifdef DEBUG_COST_A
    cost_mul--;
    cost_mod--;
    #endif
    fp_inv(ANS,A);
    mpn_mulmod_montgomery(ANS->x0,FPLIMB,ANS->x0,FPLIMB,R3,FPLIMB);
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

    if(mpn_cmp_ui(tmp1_fp.x0,FPLIMB,1)==0)        i=1;
    else if(mpn_cmp_ui(tmp1_fp.x0,FPLIMB,0)==0)    i=0;
    else                    i=-1;

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

//TODO: To montgomery
void fp_pow(fp_t *ANS,fp_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp_t tmp;
    fp_init(&tmp);//not need?

    fp_set(&tmp,A);

    for(i=1;i<length; i++){
        fp_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            fp_mul(&tmp,A,&tmp);
        }
    }
    fp_set(ANS,&tmp);

}

//TODO: consider modification "return mpn_cmp()"
int  fp_cmp(fp_t *A,fp_t *B){
    if(!mpn_cmp(A->x0,B->x0,FPLIMB)) return 0;
    else return 1;
}

int  fp_cmp_ui(fp_t *A,unsigned long int UI){
    if(!mpn_cmp_ui(A->x0,FPLIMB,UI)) return 0;
    else return 1;
}

int  fp_cmp_mpn(fp_t *A,mp_limb_t *B){
    if(!mpn_cmp(A->x0,B,FPLIMB)) return 0;
    else return 1;
}

int  fp_cmp_zero(fp_t *A){
    if(!mpn_cmp_ui(A->x0,FPLIMB,0)) return 0;
    else return 1;
}

int  fp_cmp_one(fp_t *A){
    if(!mpn_cmp_ui(A->x0,FPLIMB,1)) return 0;
    else return 1;
}

int fp_montgomery_trick_montgomery(fp_t *A_inv,fp_t *A,int n){
    int i;
    fp_t ANS[n],ALL_inv;

    fp_set(ANS,A);

    for(i=1;i<n;i++){
        fp_mulmod_montgomery(&ANS[i],&ANS[i-1],&A[i]);
    }
    fp_inv_montgomery(&ALL_inv,&ANS[n-1]);
    for(i=n-1;i>0;i--){
        fp_mulmod_montgomery(&A_inv[i],&ALL_inv,&ANS[i-1]);
        fp_mulmod_montgomery(&ALL_inv,&ALL_inv,&A[i]);
    }

    fp_set(A_inv,&ALL_inv);
    return 0;
}

void fp_lshift_ui_nonmod_single(fp_t *ANS,fp_t *A,int s){
    #ifdef DEBUG_COST_A
    cost_add_ui++;
    #endif
    mpn_lshift(ANS->x0,A->x0,FPLIMB,s);
}

void fp_lshift_ui_nonmod_double(fpd_t *ANS,fpd_t *A,int s){
    #ifdef DEBUG_COST_A
    cost_sqr++;
    #endif
    mpn_lshift(ANS->x0,A->x0,FPLIMB2,s);
}

int fp_legendre_sqrt(fp_t *ANS,fp_t *A){
    //need to 4|(p+1)
    fp_t C,D,A_tmp;
    int i;

    //legendre
    fp_pow(&C,A,sqrt_power_z);
    fp_mul(&D,&C,&C);
    fp_mul(&D,&D,A);

    if(mpn_cmp_ui(D.x0,FPLIMB,1)==0)        i=1;
    else if(mpn_cmp_ui(D.x0,FPLIMB,0)==0)    return 0;
    else                    return -1;

    //sqrt
    fp_mul(ANS,&C,A);
    return 1;
}