#include <ELiPS/bls12_g1_scm.h>
void bls12_g1_scm_basic(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    efp_t tmp_P;
	efp_init(&tmp_P);

	efp12_to_efp(&tmp_P,P);
	efp_scm(&tmp_P,&tmp_P,scalar);
	efp_to_efp12(ANS,&tmp_P);

}

void bls12_g1_scm(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P;
    efp_jacobian_t next_tmpJ_P,tmpJ_P[8],tmpJ_P_neg[8],tmpJ_P_4x[8],tmpJ_P_4x_neg[8],tmpJ_2P;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_jacobian_init(&tmpJ_2P);
    for(i=0;i<8;i++){
    efp_jacobian_init(&tmpJ_P[i]);
    efp_jacobian_init(&tmpJ_P_neg[i]);
    efp_jacobian_init(&tmpJ_P_4x[i]);
    efp_jacobian_init(&tmpJ_P_4x_neg[i]);
    }
    efp_jacobian_init(&next_tmpJ_P);

    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_jacobian_t table0[17],table1[17];
    for(i=0; i<17; i++){
        efp_jacobian_init(&table0[i]);
        efp_jacobian_init(&table1[i]);
    }

    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
	efp_to_montgomery(&tmp_P,&tmp_P);
    efp_affine_to_jacobian_montgomery(&tmpJ_P[0],&tmp_P);
    efp_ecd_jacobian_lazy_montgomery(&tmpJ_2P,&tmpJ_P[0]);
    for(i=1;i<8;i++){
	    efp_eca_jacobian_lazy_montgomery(&tmpJ_P[i],&tmpJ_P[i-1],&tmpJ_2P);
    }

	fp_t point_table[8],inv_table[8];
	for(i=0;i<8;i++)    fp_set(&point_table[i],&tmpJ_P[i].z);
    fp_montgomery_trick_montgomery(inv_table,point_table,8);
	for(i=0;i<8;i++)     efp_jacobian_to_mixture_noninv_montgomery(&tmpJ_P[i],&tmpJ_P[i],&inv_table[i]);

	for(i=0;i<8;i++){
	    efp_jacobian_set_neg(&tmpJ_P_neg[i],&tmpJ_P[i]);            //tmp_P_neg
    	efp_jacobian_skew_frobenius_map_p2_montgomery(&tmpJ_P_4x[i],&tmpJ_P[i]);        //tmp_P_4x
    	efp_jacobian_set_neg(&tmpJ_P_4x_neg[i],&tmpJ_P_4x[i]);        //tmp_P_4x_neg
    }
    //set table
    table0[0].infinity=1;                        //0
    table1[0].infinity=1;                        //0

    for(i=0;i<8;i++){
    efp_jacobian_set(&table0[i+1],&tmpJ_P[i]);                //[1]P
    efp_jacobian_set(&table0[i+9],&tmpJ_P_neg[i]);                //[-1]P
    efp_jacobian_set(&table1[i+1],&tmpJ_P_4x[i]);                //[1]P'
    efp_jacobian_set(&table1[i+9],&tmpJ_P_4x_neg[i]);                //[-1]P'
    }
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
    //naf
    int naf_length,naf_length0,naf_length1;
    int naf_binary[2][loop_length+1];
    for(i=0; i<loop_length+1; i++){
        naf_binary[0][i]=0;
        naf_binary[1][i]=0;
    }
    int *naf_pointer[2];
    naf_pointer[0]=naf_binary[0];
    naf_pointer[1]=naf_binary[1];

    naf_length0 = w_naf(naf_binary[0],s[0],5);
    naf_length1 = w_naf(naf_binary[1],s[1],5);
    if(naf_length0<naf_length1)	naf_length=naf_length1;
    else naf_length=naf_length0;

    //naf_length=loop_length-1;
    int binary0[naf_length+1],binary1[naf_length+1];

     for(i=naf_length; i>=0; i--){
        if(naf_binary[0][i]==0)         binary0[i]=0;
     	else if(naf_binary[0][i]>0)     	binary0[i]=(naf_binary[0][i]+1)>>1;
        else	binary0[i]=((17-(naf_binary[0][i]+16))>>1)+8;

        if(naf_binary[1][i]==0)         binary1[i]=0;
     	else if(naf_binary[1][i]>0)     	binary1[i]=(naf_binary[1][i]+1)>>1;
        else	binary1[i]=((17-(naf_binary[1][i]+16))>>1)+8;
    }
   if(naf_length0==naf_length1)	efp_eca_jacobian_lazy_montgomery(&next_tmpJ_P,&table0[binary0[naf_length]],&table1[binary1[naf_length]]);
   else if(naf_length0<naf_length1)	efp_jacobian_set(&next_tmpJ_P,&table1[binary1[naf_length]]);
    else efp_jacobian_set(&next_tmpJ_P,&table0[binary0[naf_length]]);

    //scm
    for(i=naf_length-1; i>=0; i--){
        efp_ecd_jacobian_lazy_montgomery(&next_tmpJ_P,&next_tmpJ_P);
        if(binary0[i]!=0)        efp_eca_mixture_lazy_montgomery(&next_tmpJ_P,&next_tmpJ_P,&table0[binary0[i]]);
        if(binary1[i]!=0)        efp_eca_mixture_lazy_montgomery(&next_tmpJ_P,&next_tmpJ_P,&table1[binary1[i]]);
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
