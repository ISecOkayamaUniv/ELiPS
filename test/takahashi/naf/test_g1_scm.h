void test_g1_scm_basic(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    efp_t tmp_P;
	efp_init(&tmp_P);

	efp12_to_efp(&tmp_P,P);
	efp_scm(&tmp_P,&tmp_P,scalar);
	efp_to_efp12(ANS,&tmp_P);

}

void test_g1_scm_w_naf(efp12_t *ANS,efp12_t *P,mpz_t scalar,int w){

    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    int window_size=w;
    int naf_point_cnt=2;
    int naf_table_size=1;
    for(int i=1;i<window_size-1;i++){
        naf_point_cnt*=2;
        naf_table_size*=2;
    }
    // printf("naf_point_cnt=%d\n",naf_point_cnt);
    // printf("naf_table_size=%d\n",naf_table_size);
    efp_t next_tmp_P,tmp_P;
    efp_jacobian_t next_tmpJ_P,tmpJ_P[naf_table_size],tmpJ_P_neg[naf_table_size],tmpJ_P_4x[naf_table_size],tmpJ_P_4x_neg[naf_table_size],tmpJ_2P;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_jacobian_init(&tmpJ_2P);
    for(i=0;i<naf_table_size;i++){
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
    efp_jacobian_t table0[naf_point_cnt+1],table1[naf_point_cnt+1];
    for(i=0; i<naf_point_cnt+1; i++){
        efp_jacobian_init(&table0[i]);
        efp_jacobian_init(&table1[i]);
    }

    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
	efp_to_montgomery(&tmp_P,&tmp_P);
    efp_affine_to_jacobian_montgomery(&tmpJ_P[0],&tmp_P);
    efp_ecd_jacobian_lazy_montgomery(&tmpJ_2P,&tmpJ_P[0]);
    for(i=1;i<naf_table_size;i++){
	    efp_eca_jacobian_lazy_montgomery(&tmpJ_P[i],&tmpJ_P[i-1],&tmpJ_2P);
    }

	fp_t point_table[naf_table_size],inv_table[naf_table_size];
	for(i=0;i<naf_table_size;i++)    fp_set(&point_table[i],&tmpJ_P[i].z);
    fp_montgomery_trick_montgomery(inv_table,point_table,naf_table_size);
	for(i=0;i<naf_table_size;i++)     efp_jacobian_to_mixture_noninv_montgomery(&tmpJ_P[i],&tmpJ_P[i],&inv_table[i]);

	for(i=0;i<naf_table_size;i++){
	    efp_jacobian_set_neg(&tmpJ_P_neg[i],&tmpJ_P[i]);            //tmp_P_neg
    	efp_jacobian_skew_frobenius_map_p2_montgomery(&tmpJ_P_4x[i],&tmpJ_P[i]);        //tmp_P_4x
    	efp_jacobian_set_neg(&tmpJ_P_4x_neg[i],&tmpJ_P_4x[i]);        //tmp_P_4x_neg
    }
    //set table
    table0[0].infinity=1;                        //0
    table1[0].infinity=1;                        //0

    for(i=0;i<naf_table_size;i++){
    efp_jacobian_set(&table0[i+1],&tmpJ_P[i]);                //[1]P
    efp_jacobian_set(&table0[i+naf_table_size+1],&tmpJ_P_neg[i]);                //[-1]P
    efp_jacobian_set(&table1[i+1],&tmpJ_P_4x[i]);                //[1]P'
    efp_jacobian_set(&table1[i+naf_table_size+1],&tmpJ_P_4x_neg[i]);                //[-1]P'
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

    naf_length0 = w_naf(naf_binary[0],s[0],window_size);
    naf_length1 = w_naf(naf_binary[1],s[1],window_size);
    if(naf_length0<naf_length1)	naf_length=naf_length1;
    else naf_length=naf_length0;

    //naf_length=loop_length-1;
    int binary0[naf_length+1],binary1[naf_length+1];

     for(i=naf_length; i>=0; i--){
        if(naf_binary[0][i]==0)         binary0[i]=0;
     	else if(naf_binary[0][i]>0)     	binary0[i]=(naf_binary[0][i]+1)>>1;
        else	binary0[i]=((naf_point_cnt+1-(naf_binary[0][i]+naf_point_cnt))>>1)+naf_table_size;

        if(naf_binary[1][i]==0)         binary1[i]=0;
     	else if(naf_binary[1][i]>0)     	binary1[i]=(naf_binary[1][i]+1)>>1;
        else	binary1[i]=((naf_point_cnt+1-(naf_binary[1][i]+naf_point_cnt))>>1)+naf_table_size;
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

void test_g1_scm(efp12_t *ANS,efp12_t *P,mpz_t scalar){

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
