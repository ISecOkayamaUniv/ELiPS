#include "ELiPS/bls12.h"

void test_g1_scm_jsf(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    efp_t temp;
    efp_jacobian_t out_tmpJ_P,next_tmpJ_P,tmpJ_P,tmpJ_P_neg,tmpJ_P_4x,tmpJ_P_4x_neg;
    efp_jacobian_init(&out_tmpJ_P);
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_init(&tmp_P_neg);
    efp_init(&tmp_P_4x);
    efp_init(&tmp_P_4x_neg);

    efp_jacobian_init(&next_tmpJ_P);
    efp_jacobian_init(&tmpJ_P);
    efp_jacobian_init(&tmpJ_P_neg);
    efp_jacobian_init(&tmpJ_P_4x);
    efp_jacobian_init(&tmpJ_P_4x_neg);

    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_jacobian_t table[9];
    for(i=0; i<9; i++){
        efp_jacobian_init(&table[i]);
    }

    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
	efp_to_montgomery(&tmp_P,&tmp_P);//add
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg //montfomery
    efp_skew_frobenius_map_p2_montgomery(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg

    //set jacobian
    efp_affine_to_jacobian_montgomery(&tmpJ_P,&tmp_P);
    efp_affine_to_jacobian_montgomery(&tmpJ_P_neg,&tmp_P_neg);
    efp_affine_to_jacobian_montgomery(&tmpJ_P_4x,&tmp_P_4x);
    efp_affine_to_jacobian_montgomery(&tmpJ_P_4x_neg,&tmp_P_4x_neg);

    //set table
    table[0].infinity=1;                        //00
    efp_affine_to_jacobian_montgomery(&table[1],&tmp_P);                //01
    efp_affine_to_jacobian_montgomery(&table[2],&tmp_P_4x);                //10
    efp_eca_jacobian_lazy_montgomery(&table[3],&tmpJ_P_4x,&tmpJ_P);       //11
    efp_affine_to_jacobian_montgomery(&table[4],&tmp_P_neg);            //0-1
    efp_affine_to_jacobian_montgomery(&table[5],&tmp_P_4x_neg);            //-10
    efp_eca_jacobian_lazy_montgomery(&table[6],&tmpJ_P_4x_neg,&tmpJ_P_neg);    //-1-1
    efp_eca_jacobian_lazy_montgomery(&table[7],&tmpJ_P_4x,&tmpJ_P_neg);        //1-1
    efp_eca_jacobian_lazy_montgomery(&table[8],&tmpJ_P_4x_neg,&tmpJ_P);        //-11

    fp_t a,b,c,d;
    fp_t ai,bi,ci,di;
    fp_t ab,cd;
    fp_t all,buf_fp;

    fp_set(&a,&table[3].z);
    fp_set(&b,&table[6].z);
    fp_set(&c,&table[7].z);
    fp_set(&d,&table[8].z);

    fp_mulmod_montgomery(&ab,&a,&b);
    fp_mulmod_montgomery(&cd,&c,&d);

	fp_mulmod_montgomery(&all,&ab,&cd);
	fp_inv_montgomery(&all,&all);

	fp_mulmod_montgomery(&buf_fp,&b,&cd);
	fp_mulmod_montgomery(&ai,&buf_fp,&all);
	fp_mulmod_montgomery(&buf_fp,&a,&cd);
	fp_mulmod_montgomery(&bi,&buf_fp,&all);

	fp_mulmod_montgomery(&buf_fp,&c,&ab);
	fp_mulmod_montgomery(&di,&buf_fp,&all);
	fp_mulmod_montgomery(&buf_fp,&d,&ab);
	fp_mulmod_montgomery(&ci,&buf_fp,&all);

	efp_jacobian_to_mixture_noninv_montgomery(&table[3],&table[3],&ai);
	efp_jacobian_to_mixture_noninv_montgomery(&table[6],&table[6],&bi);
	efp_jacobian_to_mixture_noninv_montgomery(&table[7],&table[7],&ci);
	efp_jacobian_to_mixture_noninv_montgomery(&table[8],&table[8],&di);

    //s0,s1
    mpz_neg(buf,X_z);
    mpz_pow_ui(buf,buf,2);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);

    //get loop_length
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //jsf
    int jsf_length;
    int jsf_binary[2][loop_length+1];
    char check[5];
    for(i=0; i<loop_length; i++){
        jsf_binary[0][i]=0;
        jsf_binary[1][i]=0;
    }
    int *jsf_pointer[2];
    jsf_pointer[0]=jsf_binary[0];
    jsf_pointer[1]=jsf_binary[1];
    Joint_sparse_form(jsf_pointer,s,&jsf_length);
    int binary[jsf_length+1];

    for(i=jsf_length; i>=0; i--){
        if(jsf_binary[1][i]==0 && jsf_binary[0][i]==0)         binary[i]=0;
        else if(jsf_binary[1][i]==0 && jsf_binary[0][i]==1)     binary[i]=1;
        else if(jsf_binary[1][i]==1 && jsf_binary[0][i]==0)     binary[i]=2;
        else if(jsf_binary[1][i]==1 && jsf_binary[0][i]==1)    binary[i]=3;
        else if(jsf_binary[1][i]==0 && jsf_binary[0][i]==-1)    binary[i]=4;
        else if(jsf_binary[1][i]==-1 && jsf_binary[0][i]==0)    binary[i]=5;
        else if(jsf_binary[1][i]==-1 && jsf_binary[0][i]==-1)    binary[i]=6;
        else if(jsf_binary[1][i]==1 && jsf_binary[0][i]==-1)    binary[i]=7;
        else if(jsf_binary[1][i]==-1 && jsf_binary[0][i]==1)    binary[i]=8;
    }
    efp_jacobian_set(&next_tmpJ_P,&table[binary[jsf_length]]);
    //scm
    for(i=jsf_length-1; i>=0; i--){
        //printf("binary[i]=%d:\n",binary[i]);
        efp_ecd_jacobian_lazy_montgomery(&next_tmpJ_P,&next_tmpJ_P);
        //efp_eca_jacobian_lazy(&out_tmpJ_P,&next_tmpJ_P,&table[binary[i]]);
        //efp_jacobian_printf("jacobian=",&out_tmpJ_P);printf("\n");
        efp_eca_mixture_lazy_montgomery(&next_tmpJ_P,&next_tmpJ_P,&table[binary[i]]);
        //efp_jacobian_printf("mixture=",&next_tmpJ_P);printf("\n");
        //getchar();
    }

    efp_jacobian_to_affine_montgomery(&next_tmp_P,&next_tmpJ_P);
	efp_mod_montgomery(&next_tmp_P,&next_tmp_P);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;

    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
