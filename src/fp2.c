#include <ELiPS/fp2.h>
//fp2
void fp2_init(fp2_t *A){
    fp_init(&A->x0);
    fp_init(&A->x1);
}

void fp2_printf(char *str,fp2_t *A){
    gmp_printf("%s(",str);
    fp_printf("",&A->x0);
    gmp_printf(",");
    fp_printf("",&A->x1);
    gmp_printf(")");
}

void fp2_println(char *str,fp2_t *A){
    gmp_printf("%s(",str);
    fp_printf("",&A->x0);
    gmp_printf(",");
    fp_printf("",&A->x1);
    gmp_printf(")\n");
}

void fpd2_printf(char *str,fpd2_t *A){
    gmp_printf("%s(",str);
    fpd_printf("",&A->x0);
    gmp_printf(",");
    fpd_printf("",&A->x1);
    gmp_printf(")");
}

void fpd2_println(char *str,fpd2_t *A){
    gmp_printf("%s(",str);
    fpd_printf("",&A->x0);
    gmp_printf(",");
    fpd_printf("",&A->x1);
    gmp_printf(")\n");
}

void fp2_printf_montgomery(char *str,fp2_t *A){
    gmp_printf("%s(",str);
    fp_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp_printf_montgomery("",&A->x1);
    gmp_printf(")");
}

void fp2_println_montgomery(char *str,fp2_t *A){
    gmp_printf("%s(",str);
    fp_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp_printf_montgomery("",&A->x1);
    gmp_printf(")\n");
}

void fp2_set(fp2_t *ANS,fp2_t *A){
    fp_set(&ANS->x0,&A->x0);
    fp_set(&ANS->x1,&A->x1);
}
void fpd2_set(fpd2_t *ANS,fpd2_t *A){
    fpd_set(&ANS->x0,&A->x0);
    fpd_set(&ANS->x1,&A->x1);
}

void fp2_set_ui(fp2_t *ANS,unsigned long int UI){
    fp_set_ui(&ANS->x0,UI);
    fp_set_ui(&ANS->x1,0);
}
void fp2_set_ui_ui(fp2_t *ANS,unsigned long int UI){
    fp_set_ui(&ANS->x0,UI);
    fp_set_ui(&ANS->x1,UI);
}

void fp2_set_mpn(fp2_t *ANS,mp_limb_t *A){
    fp_set_mpn(&ANS->x0,A);
    fp_set_ui(&ANS->x1,0);
}

void fp2_set_neg(fp2_t *ANS,fp2_t *A){
    fp_set_neg(&ANS->x0,&A->x0);
    fp_set_neg(&ANS->x1,&A->x1);
}
void fp2_to_montgomery(fp2_t *ANS,fp2_t *A){
    fp_to_montgomery(&ANS->x0,&A->x0);
    fp_to_montgomery(&ANS->x1,&A->x1);
}
void fp2_mod_montgomery(fp2_t *ANS,fp2_t *A){
    fp_mod_montgomery(&ANS->x0,&A->x0);
    fp_mod_montgomery(&ANS->x1,&A->x1);
}
void fp2_mod_montgomery_double(fp2_t *ANS,fpd2_t *A){
    mpn_mod_montgomery(ANS->x0.x0,FPLIMB,A->x0.x0,FPLIMB2);
    mpn_mod_montgomery(ANS->x1.x0,FPLIMB,A->x1.x0,FPLIMB2);
}
void fp2_lshift(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_lshift(&ANS->x0,&A->x0,UI);
    fp_lshift(&ANS->x1,&A->x1,UI);
}
void fp2_l1shift(fp2_t *ANS,fp2_t *A){
    fp_l1shift(&ANS->x0,&A->x0);
    fp_l1shift(&ANS->x1,&A->x1);
}
void fp2_r1shift(fp2_t *ANS,fp2_t *A)
{
    fp_r1shift(&ANS->x0, &A->x0);
    fp_r1shift(&ANS->x1, &A->x1);
}

void fp2_set_random(fp2_t *ANS,gmp_randstate_t state){
    fp_set_random(&ANS->x0,state);
    fp_set_random(&ANS->x1,state);
}
void fp2_mul(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static fp_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;

    //set
    fp_mul(&tmp1_fp,&A->x0,&B->x0);//a*c
    fp_mul(&tmp2_fp,&A->x1,&B->x1);//b*d
    fp_add(&tmp3_fp,&A->x0,&A->x1);//a+b
    fp_add(&tmp4_fp,&B->x0,&B->x1);//c+d
    //x0
    fp_sub(&ANS->x0,&tmp1_fp,&tmp2_fp);//a*c+b*d*v
    //x1
    fp_mul(&ANS->x1,&tmp3_fp,&tmp4_fp);//(a+b)(c+d)
    fp_sub(&ANS->x1,&ANS->x1,&tmp1_fp);
    fp_sub(&ANS->x1,&ANS->x1,&tmp2_fp);
}

void fp2_mul_lazy(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static fpd_t buf1,buf2;
    static fpd_t tmp1,tmp2,tmptmp;
    static fp_t tmp3,tmp4;

    //set
    fp_mul_nonmod(&tmp1,&A->x0,&B->x0);//a*c
    fp_mul_nonmod(&tmp2,&A->x1,&B->x1);//b*d

    fp_add_nonmod_single(&tmp3,&A->x0,&A->x1);
    fp_add_nonmod_single(&tmp4,&B->x0,&B->x1);

    //x0
    fp_sub_nonmod_double(&tmptmp,&tmp1,&tmp2);
    fp_mod(&ANS->x0,tmptmp.x0,FPLIMB2);

    //x1
    fp_mul_nonmod(&buf1,&tmp3,&tmp4);
    fp_sub_nonmod_double(&buf2,&buf1,&tmp1);
    fp_sub_nonmod_double(&tmptmp,&buf2,&tmp2);
    fp_mod(&ANS->x1,tmptmp.x0,FPLIMB2);
}
void fp2_mul_lazy_montgomery(fp2_t *ANS,fp2_t *A,fp2_t *B){
    static fpd_t buf1,buf2;
    static fpd_t tmp1,tmp2,tmptmp;
    static fp_t tmp3,tmp4;

    //set
    fp_mul_nonmod(&tmp1,&A->x0,&B->x0);//a*c
    fp_mul_nonmod(&tmp2,&A->x1,&B->x1);//b*d

    fp_add_nonmod_single(&tmp3,&A->x0,&A->x1);
    fp_add_nonmod_single(&tmp4,&B->x0,&B->x1);

    //x0
    fp_sub_nonmod_double(&tmptmp,&tmp1,&tmp2);
    mpn_mod_montgomery(ANS->x0.x0,FPLIMB,tmptmp.x0,FPLIMB2);

    //x1
    fp_mul_nonmod(&buf1,&tmp3,&tmp4);
    fp_sub_nonmod_double(&buf2,&buf1,&tmp1);
    fp_sub_nonmod_double(&tmptmp,&buf2,&tmp2);
    mpn_mod_montgomery(ANS->x1.x0,FPLIMB,tmptmp.x0,FPLIMB2);
}

void fp2_mul_nonmod_montgomery(fpd2_t *ANS,fp2_t *A,fp2_t *B){
    static fpd_t buf1,buf2;
    static fpd_t tmp1,tmp2,tmptmp;
    static fp_t tmp3,tmp4;

    //set
    fp_mul_nonmod(&tmp1,&A->x0,&B->x0);//a*c
    fp_mul_nonmod(&tmp2,&A->x1,&B->x1);//b*d

    fp_add_nonmod_single(&tmp3,&A->x0,&A->x1);
    fp_add_nonmod_single(&tmp4,&B->x0,&B->x1);

    //x0
    fp_sub_nonmod_double(&ANS->x0,&tmp1,&tmp2);

    //x1
    fp_mul_nonmod(&buf1,&tmp3,&tmp4);
    fp_sub_nonmod_double(&buf2,&buf1,&tmp1);
    fp_sub_nonmod_double(&ANS->x1,&buf2,&tmp2);
}
void fp2_mul_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_mul_ui(&ANS->x0,&A->x0,UI);
    fp_mul_ui(&ANS->x1,&A->x1,UI);
}
void fp2_mul_ui_nonmod_single(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_mul_ui_nonmod_single(&ANS->x0,&A->x0,UI);
    fp_mul_ui_nonmod_single(&ANS->x1,&A->x1,UI);
}
void fp2_mul_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    fp_mul_mpn(&ANS->x0,&A->x0,B);
    fp_mul_mpn(&ANS->x1,&A->x1,B);
}
void fp2_mul_mpn_montgomery(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    mpn_mulmod_montgomery(ANS->x0.x0,FPLIMB,A->x0.x0,FPLIMB,B,FPLIMB);
    mpn_mulmod_montgomery(ANS->x1.x0,FPLIMB,A->x1.x0,FPLIMB,B,FPLIMB);
}
void fp2_mul_basis(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp;
    fp_set(&tmp1_fp,&A->x0);

    fp_sub(&ANS->x0,&tmp1_fp,&A->x1);
    fp_add(&ANS->x1,&tmp1_fp,&A->x1);
}
//nonmod
void fp2_mul_basis_nonmod_single(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1;
    mpn_copyd(tmp1.x0,A->x0.x0,FPLIMB);

    fp_sub_nonmod_single(&ANS->x0,&tmp1,&A->x1);
    fp_add_nonmod_single(&ANS->x1,&tmp1,&A->x1);
}
void fp2_mul_basis_nonmod_double(fpd2_t *ANS,fpd2_t *A){
    static fpd_t tmp1;
    mpn_copyd(tmp1.x0,A->x0.x0,FPLIMB2);

    fp_sub_nonmod_double(&ANS->x0,&tmp1,&A->x1);
    fp_add_nonmod_double(&ANS->x1,&tmp1,&A->x1);
}

void fp2_inv_basis(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp,tmp2_fp;
    fp_set(&tmp1_fp,&A->x0);
    fp_set(&tmp2_fp,&A->x1);

    fp_add(&ANS->x0,&tmp1_fp,&tmp2_fp);
    fp_mul_mpn(&ANS->x0,&ANS->x0,Alpha_1_inv.x0.x0);
    fp_sub(&ANS->x1,&tmp2_fp,&tmp1_fp);
    fp_mul_mpn(&ANS->x1,&ANS->x1,Alpha_1_inv.x0.x0);
}

//Karat
void fp2_sqr(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp,tmp2_fp;
    fp_add(&tmp1_fp,&A->x0,&A->x1);
    fp_sub(&tmp2_fp,&A->x0,&A->x1);
    //x1
    fp_mul(&ANS->x1,&A->x0,&A->x1);
    fp_add(&ANS->x1,&ANS->x1,&ANS->x1);
    //x0
    fp_mul(&ANS->x0,&tmp1_fp,&tmp2_fp);
}

void fp2_sqr_lazy(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1,tmp2,tmp3,buf;
    static fp_t chk1,chk2;

    fp_add_nonmod_single(&tmp1,&A->x0,&A->x1);
    fp_sub_nonmod_single(&tmp2,&A->x0,&A->x1);

    //x1
    fp_mul(&tmp3,&A->x0,&A->x1);
    fp_add(&ANS->x1,&tmp3,&tmp3);

    //x0
    fp_mul(&ANS->x0,&tmp1,&tmp2);
}
void fp2_sqr_lazy_montgomery(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1,tmp2,tmp3,buf;
    static fp_t chk1,chk2;

    fp_add_nonmod_single(&tmp1,&A->x0,&A->x1);
    fp_sub_nonmod_single(&tmp2,&A->x0,&A->x1);

    //x1
    fp_mulmod_montgomery(&tmp3,&A->x0,&A->x1);
    fp_add(&ANS->x1,&tmp3,&tmp3);

    //x0
    fp_mulmod_montgomery(&ANS->x0,&tmp1,&tmp2);
}

void fp2_add(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_add(&ANS->x0,&A->x0,&B->x0);
    fp_add(&ANS->x1,&A->x1,&B->x1);
}

void fp2_add_nonmod_single(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_add_nonmod_single(&ANS->x0,&A->x0,&B->x0);
    fp_add_nonmod_single(&ANS->x1,&A->x1,&B->x1);
}
void fp2_add_nonmod_double(fpd2_t *ANS,fpd2_t *A,fpd2_t *B){
    fp_add_nonmod_double(&ANS->x0,&A->x0,&B->x0);
    fp_add_nonmod_double(&ANS->x1,&A->x1,&B->x1);
}
void fp2_add_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_add_ui(&ANS->x0,&A->x0,UI);
    fp_add_ui(&ANS->x1,&A->x1,0);
}
void fp2_add_ui_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_add_ui(&ANS->x0,&A->x0,UI);
    fp_add_ui(&ANS->x1,&A->x1,UI);
}
void fp2_add_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    fp_add_mpn(&ANS->x0,&A->x0,B);
    fp_add_mpn(&ANS->x1,&A->x1,B);
}
void fp2_sub(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_sub(&ANS->x0,&A->x0,&B->x0);
    fp_sub(&ANS->x1,&A->x1,&B->x1);
}
void fp2_sub_nonmod_single(fp2_t *ANS,fp2_t *A,fp2_t *B){
    fp_sub_nonmod_single(&ANS->x0,&A->x0,&B->x0);
    fp_sub_nonmod_single(&ANS->x1,&A->x1,&B->x1);
}
void fp2_sub_nonmod_double(fpd2_t *ANS,fpd2_t *A,fpd2_t *B){
    fp_sub_nonmod_double(&ANS->x0,&A->x0,&B->x0);
    fp_sub_nonmod_double(&ANS->x1,&A->x1,&B->x1);
}
void fp2_sub_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_sub_ui(&ANS->x0,&A->x0,UI);
    fp_sub_ui(&ANS->x1,&A->x1,0);
}
void fp2_sub_ui_ui(fp2_t *ANS,fp2_t *A,unsigned long int UI){
    fp_sub_ui(&ANS->x0,&A->x0,UI);
    fp_sub_ui(&ANS->x1,&A->x1,UI);
}
void fp2_sub_mpn(fp2_t *ANS,fp2_t *A,mp_limb_t *B){
    fp_sub_mpn(&ANS->x0,&A->x0,B);
    fp_sub_mpn(&ANS->x1,&A->x1,B);
}

void fp2_inv(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp,tmp2_fp,tmp3_fp,tmp4_fp;
    fp_set(&tmp1_fp,&A->x0);
    fp_set_neg(&tmp2_fp,&A->x1);

    fp_sqr(&tmp3_fp,&tmp1_fp);
    fp_mul(&tmp4_fp,&tmp2_fp,&A->x1);
    fp_sub(&tmp3_fp,&tmp3_fp,&tmp4_fp);
    fp_inv(&tmp3_fp,&tmp3_fp);
    fp_mul(&ANS->x0,&tmp1_fp,&tmp3_fp);
    fp_mul(&ANS->x1,&tmp2_fp,&tmp3_fp);
}
void fp2_inv_lazy(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp,tmp2_fp;
    static fp_t tmp3;
    static fp_t tmp1,tmp2;
    fp_set(&tmp1_fp,&A->x0);
    fp_set_neg(&tmp2_fp,&A->x1);

    fp_sqr(&tmp1,&tmp1_fp);
    fp_mul(&tmp2,&tmp2_fp,&A->x1);
    fp_sub(&tmp3,&tmp1,&tmp2);
    fp_inv(&tmp3,&tmp3);
    fp_mul(&ANS->x0,&tmp1_fp,&tmp3);
    fp_mul(&ANS->x1,&tmp2_fp,&tmp3);
}
void fp2_inv_lazy_montgomery(fp2_t *ANS,fp2_t *A){
    static fp_t tmp1_fp,tmp2_fp;
    static fp_t tmp3;
    static fp_t tmp1,tmp2;
    fp_set(&tmp1_fp,&A->x0);
    fp_set_neg(&tmp2_fp,&A->x1);

    fp_sqrmod_montgomery(&tmp1,&tmp1_fp);
    fp_mulmod_montgomery(&tmp2,&tmp2_fp,&A->x1);
    fp_sub_nonmod_single(&tmp3,&tmp1,&tmp2);
    fp_inv_montgomery(&tmp3,&tmp3);
    fp_mulmod_montgomery(&ANS->x0,&tmp1_fp,&tmp3);
    fp_mulmod_montgomery(&ANS->x1,&tmp2_fp,&tmp3);
}
int  fp2_legendre(fp2_t *A){
    fp2_t tmp;
    fp2_init(&tmp);

    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime_z,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp2_pow(&tmp,A,exp);

    if(fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

int  fp2_isCNR(fp2_t *A){
    fp2_t tmp;
    fp2_init(&tmp);
    mpz_t exp;
    mpz_init(exp);

    mpz_pow_ui(exp,prime_z,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    fp2_pow(&tmp,A,exp);

    if(fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }

}

void fp2_sqrt(fp2_t *ANS,fp2_t *A){
    fp2_t x,y,t,k,n,tmp;
    fp2_init(&x);
    fp2_init(&y);
    fp2_init(&t);
    fp2_init(&k);
    fp2_init(&n);
    fp2_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);

    fp2_set_random(&n,state);
    while(fp2_legendre(&n)!=-1){
        fp2_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,2);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    fp2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp2_pow(&x,A,exp);
    fp2_mul(&tmp,&x,&x);
    fp2_mul(&k,&tmp,A);
    fp2_mul(&x,&x,A);
    while(fp2_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fp2_pow(&tmp,&k,exp);
        while(fp2_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fp2_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fp2_pow(&t,&y,result);
        fp2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fp2_mul(&x,&x,&t);
        fp2_mul(&k,&k,&y);
    }
    fp2_set(ANS,&x);

    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}

void fp2_pow(fp2_t *ANS,fp2_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp2_t tmp;
    fp2_init(&tmp);
    fp2_set(&tmp,A);

    for(i=1;i<length; i++){
        fp2_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            fp2_mul(&tmp,A,&tmp);
        }
    }

    fp2_set(ANS,&tmp);
}

int  fp2_cmp(fp2_t *A,fp2_t *B){
    if(!fp_cmp(&A->x0,&B->x0) && !fp_cmp(&A->x1,&B->x1)) return 0;
    else return 1;
}

int  fp2_cmp_ui(fp2_t *A,unsigned long int UI){
    if(!fp_cmp_ui(&A->x0,UI) && !fp_cmp_ui(&A->x1,UI)) return 0;
    return 1;
}

int  fp2_cmp_mpn(fp2_t *A,mp_limb_t *B){
    if(!fp_cmp_mpn(&A->x0,B) && !fp_cmp_mpn(&A->x1,B)) return 0;
    else return 1;
}

int  fp2_cmp_zero(fp2_t *A){
    if(!fp_cmp_zero(&A->x0) && !fp_cmp_zero(&A->x1)) return 0;
    return 1;
}

int  fp2_cmp_one(fp2_t *A){
    if(!fp_cmp_one(&A->x0) && !fp_cmp_zero(&A->x1)) return 0;
    else return 1;
}
int fp2_montgomery_trick(fp2_t *A_inv,fp2_t *A,int n){
    int i;
    fp2_t ANS[n],ALL_inv;
    fp2_set(ANS,A);

    for(i=1;i<n;i++){
        fp2_mul_lazy(&ANS[i],&ANS[i-1],&A[i]);
    }
    fp2_inv_lazy(&ALL_inv,&ANS[n-1]);
    for(i=n-1;i>0;i--){
        fp2_mul_lazy(&A_inv[i],&ALL_inv,&ANS[i-1]);
        fp2_mul_lazy(&ALL_inv,&ALL_inv,&A[i]);
    }

    fp2_set(A_inv,&ALL_inv);
    return 0;
}
int fp2_montgomery_trick_montgomery(fp2_t *A_inv,fp2_t *A,int n){
    int i;
    fp2_t ANS[n],ALL_inv;
    fp2_set(ANS,A);

    for(i=1;i<n;i++){
        fp2_mul_lazy_montgomery(&ANS[i],&ANS[i-1],&A[i]);
    }
    fp2_inv_lazy_montgomery(&ALL_inv,&ANS[n-1]);
    for(i=n-1;i>0;i--){
        fp2_mul_lazy_montgomery(&A_inv[i],&ALL_inv,&ANS[i-1]);
        fp2_mul_lazy_montgomery(&ALL_inv,&ALL_inv,&A[i]);
    }

    fp2_set(A_inv,&ALL_inv);
    return 0;
}

void fp2_sqr_nonmod_montgomery(fpd2_t *ANS,fp2_t *A){//add
    static fp_t tmp1_fp,tmp2_fp;
    fp_add_nonmod_single(&tmp1_fp,&A->x0,&A->x1);
    fp_sub_nonmod_single(&tmp2_fp,&A->x0,&A->x1);
    //x1
    fp_mul_nonmod(&ANS->x1,&A->x0,&A->x1);
    fp_add_nonmod_double(&ANS->x1,&ANS->x1,&ANS->x1);
    //x0
    fp_mul_nonmod(&ANS->x0,&tmp1_fp,&tmp2_fp);
}

void fp2_lshift_ui_nonmod_single(fp2_t *ANS,fp2_t *A,int s){
    fp_lshift_ui_nonmod_single(&ANS->x0,&A->x0,s);
    fp_lshift_ui_nonmod_single(&ANS->x1,&A->x1,s);
}
void fp2_lshift_ui_nonmod_double(fpd2_t *ANS,fpd2_t *A,int s){
    fp_lshift_ui_nonmod_double(&ANS->x0,&A->x0,s);
    fp_lshift_ui_nonmod_double(&ANS->x1,&A->x1,s);
}

void fp2_mul_3_twist_b(fp2_t *ANS,fp2_t *A){
    static fp2_t tmp;
    fp2_lshift_ui_nonmod_single(&tmp,A,3);

    fp_sub_nonmod_single(&ANS->x0,&tmp.x0,&tmp.x1);
    fp_add_nonmod_single(&ANS->x1,&tmp.x1,&tmp.x0);
}
