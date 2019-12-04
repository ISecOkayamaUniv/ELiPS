#include <ELiPS/bls12_g1_scm.h>
void bls12_g1_scm_plain(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    efp_t tmp_P;
	efp_init(&tmp_P);
	
	efp12_to_efp(&tmp_P,P);
	efp_scm(&tmp_P,&tmp_P,scalar);
	efp_to_efp12(ANS,&tmp_P);
	
}
void bls12_g1_scm_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    efp_t tmp_P;
	efp_init(&tmp_P);
	
	efp12_to_efp(&tmp_P,P);
	efp_scm_lazy(&tmp_P,&tmp_P,scalar);
	efp_to_efp12(ANS,&tmp_P);
	
}
void bls12_g1_scm_2split(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_4x;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_init(&tmp_P_4x);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_t table[4];
    for(i=0; i<4; i++){
        efp_init(&table[i]);
    }
    
    //set
    efp12_to_efp(&tmp_P,P);                //tmp_P
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);//tmp_P_4x
    
    //set table
    table[0].infinity=1;                        //00
    efp_set(&table[1],&tmp_P);            //01
    efp_set(&table[2],&tmp_P_4x);            //10
    efp_eca(&table[3],&tmp_P,&tmp_P_4x);    //11
    
    //s0,s1
    mpz_set(buf,X_z);//neg->set
    mpz_pow_ui(buf,buf,2);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    
    //binary
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        //printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }

    }
    //set binary
    char binary_s[2][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<2; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }else{
            char binary_buf[length_s[i]];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
        binary[i]=strtol(str,&e,2);
    }
    efp_set(&next_tmp_P,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        efp_ecd(&next_tmp_P,&next_tmp_P);
        efp_eca(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
    }
    
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}

void bls12_g1_scm_2split_2naf(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_init(&tmp_P_neg);
    efp_init(&tmp_P_4x);
    efp_init(&tmp_P_4x_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_t table[9];
    for(i=0; i<9; i++){
        efp_init(&table[i]);
    }
    
    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg
    
    //set table
    table[0].infinity=1;                        //00
    efp_set(&table[1],&tmp_P);                //01
    efp_set(&table[2],&tmp_P_4x);                //10
    efp_eca(&table[3],&tmp_P_4x,&tmp_P);        //11
    efp_set(&table[4],&tmp_P_neg);            //0-1
    efp_set(&table[5],&tmp_P_4x_neg);            //-10
    efp_eca(&table[6],&tmp_P_4x_neg,&tmp_P_neg);    //-1-1
    efp_eca(&table[7],&tmp_P_4x,&tmp_P_neg);        //1-1
    efp_eca(&table[8],&tmp_P_4x_neg,&tmp_P);        //-11

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
    int naf_length;
    int naf_binary[2][loop_length+1];
    char check[5];
    for(i=0; i<loop_length+1; i++){
        naf_binary[0][i]=0;
        naf_binary[1][i]=0;
    }
    int *naf_pointer[2];
    naf_pointer[0]=naf_binary[0];
    naf_pointer[1]=naf_binary[1];
    w_naf(naf_binary[0],s[0],2);
    w_naf(naf_binary[1],s[1],2);
    
    naf_length=loop_length;
    int binary[naf_length+1];
    for(i=naf_length; i>=0; i--){
        if(naf_binary[1][i]==0 && naf_binary[0][i]==0)         binary[i]=0;
        else if(naf_binary[1][i]==0 && naf_binary[0][i]==1)     binary[i]=1;
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==0)     binary[i]=2;
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==1)    binary[i]=3;
        else if(naf_binary[1][i]==0 && naf_binary[0][i]==-1)    binary[i]=4;
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==0)    binary[i]=5;
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==-1)    binary[i]=6;
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==-1)    binary[i]=7;
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==1)    binary[i]=8;
    }
    efp_set(&next_tmp_P,&table[binary[naf_length]]);
        //scm
    for(i=naf_length-1; i>=0; i--){
        efp_ecd(&next_tmp_P,&next_tmp_P);
        efp_eca(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
    }
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_3naf_shamia(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    efp_t tmp_3P,tmp_3P_neg,tmp_3P_4x,tmp_3P_4x_neg;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_init(&tmp_P_neg);
    efp_init(&tmp_P_4x);
    efp_init(&tmp_P_4x_neg);
    
    efp_init(&tmp_3P);
    efp_init(&tmp_3P_neg);
    efp_init(&tmp_3P_4x);
    efp_init(&tmp_3P_4x_neg);
    
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_t table[25];
    for(i=0; i<25; i++){
        efp_init(&table[i]);
    }
    
    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg
    
    efp_ecd(&tmp_3P,&tmp_P);
    efp_eca(&tmp_3P,&tmp_3P,&tmp_P);                    //tmp_3P
    efp_set_neg(&tmp_3P_neg,&tmp_3P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_3P_4x,&tmp_3P);        //tmp_P_4x
    efp_set_neg(&tmp_3P_4x_neg,&tmp_3P_4x);        //tmp_P_4x_neg
    
    //set table
    table[0].infinity=1;                        //00
    efp_set(&table[1],&tmp_P);                //01
    efp_set(&table[2],&tmp_3P);                //03
    efp_set(&table[3],&tmp_P_neg);                //0-1
    efp_set(&table[4],&tmp_3P_neg);                //0-3
    
    efp_set(&table[5],&tmp_P_4x);                //10
    efp_eca(&table[6],&tmp_P_4x,&tmp_P);        //11
    efp_eca(&table[7],&tmp_P_4x,&tmp_3P);        //13
    efp_eca(&table[8],&tmp_P_4x,&tmp_P_neg);    //1-1
    efp_eca(&table[9],&tmp_P_4x,&tmp_3P_neg);        //1-3
    
    efp_set(&table[10],&tmp_P_4x_neg);                //-10
    efp_set_neg(&table[11],&table[8]);         //-11
    efp_set_neg(&table[12],&table[9]);         //-13
    efp_set_neg(&table[13],&table[6]);         //-1-1
    efp_set_neg(&table[14],&table[7]);        //-1-3
    
    efp_set(&table[15],&tmp_3P_4x);                //30
    efp_eca(&table[16],&tmp_3P_4x,&tmp_P);        //31
    efp_eca(&table[17],&tmp_3P_4x,&tmp_3P);        //33
    efp_eca(&table[18],&tmp_3P_4x,&tmp_P_neg);    //3-1
    efp_eca(&table[19],&tmp_3P_4x,&tmp_3P_neg);        //3-3
    
    efp_set(&table[20],&tmp_3P_4x_neg);                //-30
    efp_set_neg(&table[21],&table[18]);         //-31
    efp_set_neg(&table[22],&table[19]);         //-33
    efp_set_neg(&table[23],&table[16]);    //-3-1
    efp_set_neg(&table[24],&table[17]);        //-3-3

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
    
    /*
    struct timeval tv_A,tv_B;
    gettimeofday(&tv_A,NULL);
    naf_length0 = w_naf(naf_binary[0],s[0],3);
    naf_length1 = w_naf(naf_binary[1],s[1],3);
    gettimeofday(&tv_B,NULL);
    printf("naf     : %.4f[ms]\n",timedifference_msec(tv_A,tv_B));
    */
    naf_length0 = w_naf(naf_binary[0],s[0],3);
    naf_length1 = w_naf(naf_binary[1],s[1],3);
    if(naf_length0<naf_length1)	naf_length=naf_length1;
    else naf_length=naf_length0;
    
    //naf_length=loop_length-1;
    int binary[naf_length+1];
    for(i=naf_length; i>=0; i--){
        if(naf_binary[1][i]==0 && naf_binary[0][i]==0)         binary[i]=0;
        else if(naf_binary[1][i]==0 && naf_binary[0][i]==1)     binary[i]=1;
        else if(naf_binary[1][i]==0 && naf_binary[0][i]==3)     binary[i]=2;
        else if(naf_binary[1][i]==0 && naf_binary[0][i]==-1)     binary[i]=3;
        else if(naf_binary[1][i]==0 && naf_binary[0][i]==-3)     binary[i]=4;
        
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==0)     binary[i]=5;
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==1)     binary[i]=6;
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==3)     binary[i]=7;
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==-1)     binary[i]=8;
        else if(naf_binary[1][i]==1 && naf_binary[0][i]==-3)     binary[i]=9;
        
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==0)     binary[i]=10;
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==1)     binary[i]=11;
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==3)     binary[i]=12;
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==-1)     binary[i]=13;
        else if(naf_binary[1][i]==-1 && naf_binary[0][i]==-3)     binary[i]=14;
        
        else if(naf_binary[1][i]==3 && naf_binary[0][i]==0)     binary[i]=15;
        else if(naf_binary[1][i]==3 && naf_binary[0][i]==1)     binary[i]=16;
        else if(naf_binary[1][i]==3 && naf_binary[0][i]==3)     binary[i]=17;
        else if(naf_binary[1][i]==3 && naf_binary[0][i]==-1)     binary[i]=18;
        else if(naf_binary[1][i]==3 && naf_binary[0][i]==-3)     binary[i]=19;
        
        else if(naf_binary[1][i]==-3 && naf_binary[0][i]==0)     binary[i]=20;
        else if(naf_binary[1][i]==-3 && naf_binary[0][i]==1)     binary[i]=21;
        else if(naf_binary[1][i]==-3 && naf_binary[0][i]==3)     binary[i]=22;
        else if(naf_binary[1][i]==-3 && naf_binary[0][i]==-1)     binary[i]=23;
        else if(naf_binary[1][i]==-3 && naf_binary[0][i]==-3)     binary[i]=24;
    }
    efp_set(&next_tmp_P,&table[binary[naf_length]]);
        //scm
    //int ADD=0;
    for(i=naf_length-1; i>=0; i--){
        efp_ecd(&next_tmp_P,&next_tmp_P);
        efp_eca(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
        //if(binary[i]!=0)ADD++;
    }
    //printf("naf_shamre ADD:%d\n",ADD);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_3naf_interleaving(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    efp_t tmp_3P,tmp_3P_neg,tmp_3P_4x,tmp_3P_4x_neg;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_init(&tmp_P_neg);
    efp_init(&tmp_P_4x);
    efp_init(&tmp_P_4x_neg);
    
    efp_init(&tmp_3P);
    efp_init(&tmp_3P_neg);
    efp_init(&tmp_3P_4x);
    efp_init(&tmp_3P_4x_neg);
    
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_t table0[5],table1[5];
    for(i=0; i<5; i++){
        efp_init(&table0[i]);
        efp_init(&table1[i]);
    }
    
    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg
    
    efp_ecd(&tmp_3P,&tmp_P);
    efp_eca(&tmp_3P,&tmp_3P,&tmp_P);                    //tmp_3P
    efp_set_neg(&tmp_3P_neg,&tmp_3P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_3P_4x,&tmp_3P);        //tmp_P_4x
    efp_set_neg(&tmp_3P_4x_neg,&tmp_3P_4x);        //tmp_P_4x_neg
    
    //set table
    table0[0].infinity=1;                        //0
    efp_set(&table0[1],&tmp_P);                //[1]P
    efp_set(&table0[2],&tmp_3P);                //[3]P
    efp_set(&table0[3],&tmp_P_neg);                //[-1]P
    efp_set(&table0[4],&tmp_3P_neg);                //[-3]P
    
    table1[0].infinity=1;                        //0
    efp_set(&table1[1],&tmp_P_4x);                //[1]P'
    efp_set(&table1[2],&tmp_3P_4x);                //[3]P'
    efp_set(&table1[3],&tmp_P_4x_neg);                //[-1]P'
    efp_set(&table1[4],&tmp_3P_4x_neg);                //[-3]P'

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
    /*
    struct timeval tv_A,tv_B;
    gettimeofday(&tv_A,NULL);
    naf_length0 = w_naf(naf_binary[0],s[0],3);
    naf_length1 = w_naf(naf_binary[1],s[1],3);
    gettimeofday(&tv_B,NULL);
    printf("naf     : %.4f[ms]\n",timedifference_msec(tv_A,tv_B));
    */
    naf_length0 = w_naf(naf_binary[0],s[0],3);
    naf_length1 = w_naf(naf_binary[1],s[1],3);
    if(naf_length0<naf_length1)	naf_length=naf_length1;
    else naf_length=naf_length0;
    
    //naf_length=loop_length-1;
    int binary0[naf_length+1],binary1[naf_length+1];
    for(i=naf_length; i>=0; i--){
        if(naf_binary[0][i]==0)         binary0[i]=0;
        else if(naf_binary[0][i]==1)     binary0[i]=1;
        else if(naf_binary[0][i]==3)     binary0[i]=2;
        else if(naf_binary[0][i]==-1)     binary0[i]=3;
        else if(naf_binary[0][i]==-3)     binary0[i]=4;
        if(naf_binary[1][i]==0)         binary1[i]=0;
        else if(naf_binary[1][i]==1)     binary1[i]=1;
        else if(naf_binary[1][i]==3)     binary1[i]=2;
        else if(naf_binary[1][i]==-1)     binary1[i]=3;
        else if(naf_binary[1][i]==-3)     binary1[i]=4;
    }
   if(naf_length0==naf_length1)	efp_eca(&next_tmp_P,&table0[binary0[naf_length]],&table1[binary1[naf_length]]);
   else if(naf_length0<naf_length1)	efp_set(&next_tmp_P,&table1[binary1[naf_length]]);
    else efp_set(&next_tmp_P,&table0[binary0[naf_length]]);
    
    //int ADD=0;
        //scm
    for(i=naf_length-1; i>=0; i--){
        efp_ecd(&next_tmp_P,&next_tmp_P);
        efp_eca(&next_tmp_P,&next_tmp_P,&table0[binary0[i]]);
        efp_eca(&next_tmp_P,&next_tmp_P,&table1[binary1[i]]);
        //if(binary0[i]!=0)ADD++;
        //if(binary1[i]!=0)ADD++;
    }
    //printf("naf_interleaving ADD:%d\n",ADD);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_5naf_interleaving(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_2P;
    efp_t tmp_P[8],tmp_P_neg[8],tmp_P_4x[8],tmp_P_4x_neg[8];
    efp_init(&next_tmp_P);
    efp_init(&tmp_2P);
    for(i=0;i<8;i++){
    efp_init(&tmp_P[i]);
    efp_init(&tmp_P_neg[i]);
    efp_init(&tmp_P_4x[i]);
    efp_init(&tmp_P_4x_neg[i]);
    }    
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_t table0[17],table1[17];
    for(i=0; i<17; i++){
        efp_init(&table0[i]);
        efp_init(&table1[i]);
    }
    //set
    efp12_to_efp(&tmp_P[0],P);                    //tmp_P
    efp_ecd(&tmp_2P,&tmp_P[0]);
    efp_set_neg(&tmp_P_neg[0],&tmp_P[0]);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x[0],&tmp_P[0]);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg[0],&tmp_P_4x[0]);        //tmp_P_4x_neg
    for(i=1;i<8;i++){
    efp_eca(&tmp_P[i],&tmp_P[i-1],&tmp_2P);                    //tmp_3P
    efp_set_neg(&tmp_P_neg[i],&tmp_P[i]);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x[i],&tmp_P[i]);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg[i],&tmp_P_4x[i]);        //tmp_P_4x_neg
    }
    //set table
    table0[0].infinity=1;                        //0
    table1[0].infinity=1;                        //0
    
    for(i=0;i<8;i++){
    efp_set(&table0[i+1],&tmp_P[i]);                //[1]P
    efp_set(&table0[i+9],&tmp_P_neg[i]);                //[-1]P
    efp_set(&table1[i+1],&tmp_P_4x[i]);                //[1]P'
    efp_set(&table1[i+9],&tmp_P_4x_neg[i]);                //[-1]P'
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
    /*
    struct timeval tv_A,tv_B;
    gettimeofday(&tv_A,NULL);
    naf_length0 = w_naf(naf_binary[0],s[0],3);
    naf_length1 = w_naf(naf_binary[1],s[1],3);
    gettimeofday(&tv_B,NULL);
    printf("naf     : %.4f[ms]\n",timedifference_msec(tv_A,tv_B));
    */
    naf_length0 = w_naf(naf_binary[0],s[0],5);
    naf_length1 = w_naf(naf_binary[1],s[1],5);
    if(naf_length0<naf_length1)	naf_length=naf_length1;
    else naf_length=naf_length0;
    
    //naf_length=loop_length-1;
    int binary0[naf_length+1],binary1[naf_length+1];
    
    for(i=naf_length; i>=0; i--){
        if(naf_binary[0][i]==0)         binary0[i]=0;
        else if(naf_binary[0][i]==1)     binary0[i]=1;
        else if(naf_binary[0][i]==3)     binary0[i]=2;
        else if(naf_binary[0][i]==5)     binary0[i]=3;
        else if(naf_binary[0][i]==7)     binary0[i]=4;
        else if(naf_binary[0][i]==9)     binary0[i]=5;
        else if(naf_binary[0][i]==11)     binary0[i]=6;
        else if(naf_binary[0][i]==13)     binary0[i]=7;
        else if(naf_binary[0][i]==15)     binary0[i]=8;
        
        else if(naf_binary[0][i]==-1)     binary0[i]=9;
        else if(naf_binary[0][i]==-3)     binary0[i]=10;
        else if(naf_binary[0][i]==-5)     binary0[i]=11;
        else if(naf_binary[0][i]==-7)     binary0[i]=12;
        else if(naf_binary[0][i]==-9)     binary0[i]=13;
        else if(naf_binary[0][i]==-11)     binary0[i]=14;
        else if(naf_binary[0][i]==-13)     binary0[i]=15;
        else if(naf_binary[0][i]==-15)     binary0[i]=16;
        
        if(naf_binary[1][i]==0)         binary1[i]=0;
        else if(naf_binary[1][i]==1)     binary1[i]=1;
        else if(naf_binary[1][i]==3)     binary1[i]=2;
        else if(naf_binary[1][i]==5)     binary1[i]=3;
        else if(naf_binary[1][i]==7)     binary1[i]=4;
        else if(naf_binary[1][i]==9)     binary1[i]=5;
        else if(naf_binary[1][i]==11)     binary1[i]=6;
        else if(naf_binary[1][i]==13)     binary1[i]=7;
        else if(naf_binary[1][i]==15)     binary1[i]=8;
        
        else if(naf_binary[1][i]==-1)     binary1[i]=9;
        else if(naf_binary[1][i]==-3)     binary1[i]=10;
        else if(naf_binary[1][i]==-5)     binary1[i]=11;
        else if(naf_binary[1][i]==-7)     binary1[i]=12;
        else if(naf_binary[1][i]==-9)     binary1[i]=13;
        else if(naf_binary[1][i]==-11)     binary1[i]=14;
        else if(naf_binary[1][i]==-13)     binary1[i]=15;
        else if(naf_binary[1][i]==-15)     binary1[i]=16;
    }
   if(naf_length0==naf_length1)	efp_eca(&next_tmp_P,&table0[binary0[naf_length]],&table1[binary1[naf_length]]);
   else if(naf_length0<naf_length1)	efp_set(&next_tmp_P,&table1[binary1[naf_length]]);
    else efp_set(&next_tmp_P,&table0[binary0[naf_length]]);
    
    //int ADD=0;
        //scm
    for(i=naf_length-1; i>=0; i--){
        efp_ecd(&next_tmp_P,&next_tmp_P);
        efp_eca(&next_tmp_P,&next_tmp_P,&table0[binary0[i]]);
        efp_eca(&next_tmp_P,&next_tmp_P,&table1[binary1[i]]);
        //if(binary0[i]!=0)ADD++;
        //if(binary1[i]!=0)ADD++;
    }
    //printf("naf_interleaving ADD:%d\n",ADD);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_jsf(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_init(&tmp_P_neg);
    efp_init(&tmp_P_4x);
    efp_init(&tmp_P_4x_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_t table[9];
    for(i=0; i<9; i++){
        efp_init(&table[i]);
    }
    
    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg
    //set table
    table[0].infinity=1;                        //00
    efp_set(&table[1],&tmp_P);                //01
    efp_set(&table[2],&tmp_P_4x);                //10
    efp_eca(&table[3],&tmp_P_4x,&tmp_P);        //11
    efp_set(&table[4],&tmp_P_neg);            //0-1
    efp_set(&table[5],&tmp_P_4x_neg);            //-10
    efp_eca(&table[6],&tmp_P_4x_neg,&tmp_P_neg);    //-1-1
    efp_eca(&table[7],&tmp_P_4x,&tmp_P_neg);        //1-1
    efp_eca(&table[8],&tmp_P_4x_neg,&tmp_P);        //-11
    
    //print table
//    for(i=0;i<9;i++){
//	printf("table[%d]:",i);efp_printf(&table[i],"");printf("\n");
    //}

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
    
    /*
    struct timeval tv_A,tv_B;
    gettimeofday(&tv_A,NULL);
    Joint_sparse_form(jsf_pointer,s,&jsf_length);
    gettimeofday(&tv_B,NULL);
    printf("jsf     : %.4f[ms]\n",timedifference_msec(tv_A,tv_B));
    */
    
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
    efp_set(&next_tmp_P,&table[binary[jsf_length]]);
    //scm
    //int ADD=0;
    for(i=jsf_length-1; i>=0; i--){
        efp_ecd(&next_tmp_P,&next_tmp_P);
        efp_eca(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
        //if(binary[i]!=0)ADD++;
    }
    //printf("jsf ADD:%d\n",ADD);
    
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_jsf_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_init(&tmp_P_neg);
    efp_init(&tmp_P_4x);
    efp_init(&tmp_P_4x_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    efp_t table[9];
    for(i=0; i<9; i++){
        efp_init(&table[i]);
    }
    
    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg
    //set table
    table[0].infinity=1;                        //00
    efp_set(&table[1],&tmp_P);                //01
    efp_set(&table[2],&tmp_P_4x);                //10
    efp_eca_lazy(&table[3],&tmp_P_4x,&tmp_P);        //11
    efp_set(&table[4],&tmp_P_neg);            //0-1
    efp_set(&table[5],&tmp_P_4x_neg);            //-10
    efp_eca_lazy(&table[6],&tmp_P_4x_neg,&tmp_P_neg);    //-1-1
    efp_eca_lazy(&table[7],&tmp_P_4x,&tmp_P_neg);        //1-1
    efp_eca_lazy(&table[8],&tmp_P_4x_neg,&tmp_P);        //-11
    
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
    efp_set(&next_tmp_P,&table[binary[jsf_length]]);
    //scm
    for(i=jsf_length-1; i>=0; i--){
        efp_ecd_lazy(&next_tmp_P,&next_tmp_P);
        efp_eca_lazy(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
    }
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_jsf_jacobian_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    efp_t temp;
    efp_jacobian_t next_tmpJ_P,tmpJ_P,tmpJ_P_neg,tmpJ_P_4x,tmpJ_P_4x_neg;
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
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg

    //set jacobian
    efp_affine_to_jacobian(&tmpJ_P,&tmp_P);
    efp_affine_to_jacobian(&tmpJ_P_neg,&tmp_P_neg);
    efp_affine_to_jacobian(&tmpJ_P_4x,&tmp_P_4x);
    efp_affine_to_jacobian(&tmpJ_P_4x_neg,&tmp_P_4x_neg);

    //set table
    table[0].infinity=1;                        //00
    efp_affine_to_jacobian(&table[1],&tmp_P);                //01
    efp_affine_to_jacobian(&table[2],&tmp_P_4x);                //10
    efp_eca_jacobian_lazy(&table[3],&tmpJ_P_4x,&tmpJ_P);       //11
    efp_affine_to_jacobian(&table[4],&tmp_P_neg);            //0-1
    efp_affine_to_jacobian(&table[5],&tmp_P_4x_neg);            //-10
    efp_eca_jacobian_lazy(&table[6],&tmpJ_P_4x_neg,&tmpJ_P_neg);    //-1-1
    efp_eca_jacobian_lazy(&table[7],&tmpJ_P_4x,&tmpJ_P_neg);        //1-1
    efp_eca_jacobian_lazy(&table[8],&tmpJ_P_4x_neg,&tmpJ_P);        //-11
    

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
        efp_ecd_jacobian_lazy(&next_tmpJ_P,&next_tmpJ_P);
        efp_eca_jacobian_lazy(&next_tmpJ_P,&next_tmpJ_P,&table[binary[i]]);
    }
    
    efp_jacobian_to_affine(&next_tmp_P,&next_tmpJ_P);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}

void bls12_g1_scm_2split_5naf_interleaving_mixture(efp12_t *ANS,efp12_t *P,mpz_t scalar){
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
    efp_affine_to_jacobian(&tmpJ_P[0],&tmp_P);
    efp_ecd_jacobian(&tmpJ_2P,&tmpJ_P[0]);
    for(i=1;i<8;i++){
	    efp_eca_jacobian(&tmpJ_P[i],&tmpJ_P[i-1],&tmpJ_2P);
    }
	
	fp_t point_table[8],inv_table[8];
	for(i=0;i<8;i++)    fp_set(&point_table[i],&tmpJ_P[i].z);
    fp_montgomery_trick(inv_table,point_table,8);
	for(i=0;i<8;i++)     efp_mix(&tmpJ_P[i],&tmpJ_P[i],&inv_table[i]);
	
	for(i=0;i<8;i++){
	    efp_jacobian_set_neg(&tmpJ_P_neg[i],&tmpJ_P[i]);            //tmp_P_neg
    	efp_jacobian_skew_frobenius_map_p2(&tmpJ_P_4x[i],&tmpJ_P[i]);        //tmp_P_4x
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
    /*
    struct timeval tv_A,tv_B;
    gettimeofday(&tv_A,NULL);
    naf_length0 = w_naf(naf_binary[0],s[0],3);
    naf_length1 = w_naf(naf_binary[1],s[1],3);
    gettimeofday(&tv_B,NULL);
    printf("naf     : %.4f[ms]\n",timedifference_msec(tv_A,tv_B));
    */
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
   if(naf_length0==naf_length1)	efp_eca_jacobian(&next_tmpJ_P,&table0[binary0[naf_length]],&table1[binary1[naf_length]]);
   else if(naf_length0<naf_length1)	efp_jacobian_set(&next_tmpJ_P,&table1[binary1[naf_length]]);
    else efp_jacobian_set(&next_tmpJ_P,&table0[binary0[naf_length]]);
    
    //scm
    for(i=naf_length-1; i>=0; i--){
        efp_ecd_jacobian(&next_tmpJ_P,&next_tmpJ_P);
        if(binary0[i]!=0)        efp_eca_mixture(&next_tmpJ_P,&next_tmpJ_P,&table0[binary0[i]]);
        if(binary1[i]!=0)        efp_eca_mixture(&next_tmpJ_P,&next_tmpJ_P,&table1[binary1[i]]);
    }
    
    efp_jacobian_to_affine(&next_tmp_P,&next_tmpJ_P);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_jsf_mixture_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar){
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
    efp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    efp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    efp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg

    //set jacobian
    efp_affine_to_jacobian(&tmpJ_P,&tmp_P);
    efp_affine_to_jacobian(&tmpJ_P_neg,&tmp_P_neg);
    efp_affine_to_jacobian(&tmpJ_P_4x,&tmp_P_4x);
    efp_affine_to_jacobian(&tmpJ_P_4x_neg,&tmp_P_4x_neg);

    //set table
    table[0].infinity=1;                        //00
    efp_affine_to_jacobian(&table[1],&tmp_P);                //01
    efp_affine_to_jacobian(&table[2],&tmp_P_4x);                //10
    efp_eca_jacobian_lazy(&table[3],&tmpJ_P_4x,&tmpJ_P);       //11
    efp_affine_to_jacobian(&table[4],&tmp_P_neg);            //0-1
    efp_affine_to_jacobian(&table[5],&tmp_P_4x_neg);            //-10
    efp_eca_jacobian_lazy(&table[6],&tmpJ_P_4x_neg,&tmpJ_P_neg);    //-1-1
    efp_eca_jacobian_lazy(&table[7],&tmpJ_P_4x,&tmpJ_P_neg);        //1-1
    efp_eca_jacobian_lazy(&table[8],&tmpJ_P_4x_neg,&tmpJ_P);        //-11
    
    fp_t a,b,c,d;
    fp_t ai,bi,ci,di;
    fp_t ab,cd;
    fp_t all,buf_fp;
    
    fp_set(&a,&table[3].z);
    fp_set(&b,&table[6].z);
    fp_set(&c,&table[7].z);
    fp_set(&d,&table[8].z);
    
    fp_mul(&ab,&a,&b);
    fp_mul(&cd,&c,&d);
    
	fp_mul(&all,&ab,&cd);
	fp_inv(&all,&all);
	
	fp_mul(&buf_fp,&b,&cd);
	fp_mul(&ai,&buf_fp,&all);
	fp_mul(&buf_fp,&a,&cd);
	fp_mul(&bi,&buf_fp,&all);
	
	fp_mul(&buf_fp,&c,&ab);
	fp_mul(&di,&buf_fp,&all);
	fp_mul(&buf_fp,&d,&ab);
	fp_mul(&ci,&buf_fp,&all);
	
	efp_mix(&table[3],&table[3],&ai);
	efp_mix(&table[6],&table[6],&bi);
	efp_mix(&table[7],&table[7],&ci);
	efp_mix(&table[8],&table[8],&di);
	
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
        efp_ecd_jacobian_lazy(&next_tmpJ_P,&next_tmpJ_P);
        //efp_eca_jacobian_lazy(&out_tmpJ_P,&next_tmpJ_P,&table[binary[i]]);
        //efp_jacobian_printf("jacobian=",&out_tmpJ_P);printf("\n");
        efp_eca_mixture_lazy(&next_tmpJ_P,&next_tmpJ_P,&table[binary[i]]);
        //efp_jacobian_printf("mixture=",&next_tmpJ_P);printf("\n");
        //getchar();
    }
    
    efp_jacobian_to_affine(&next_tmp_P,&next_tmpJ_P);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_5naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar){
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
    efp_affine_to_jacobian(&tmpJ_P[0],&tmp_P);
    efp_ecd_jacobian_lazy(&tmpJ_2P,&tmpJ_P[0]);
    for(i=1;i<8;i++){
	    efp_eca_jacobian_lazy(&tmpJ_P[i],&tmpJ_P[i-1],&tmpJ_2P);
    }
	
	fp_t point_table[8],inv_table[8];
	for(i=0;i<8;i++)    fp_set(&point_table[i],&tmpJ_P[i].z);
    fp_montgomery_trick(inv_table,point_table,8);
	for(i=0;i<8;i++)     efp_mix(&tmpJ_P[i],&tmpJ_P[i],&inv_table[i]);
	
	for(i=0;i<8;i++){
	    efp_jacobian_set_neg(&tmpJ_P_neg[i],&tmpJ_P[i]);            //tmp_P_neg
    	efp_jacobian_skew_frobenius_map_p2(&tmpJ_P_4x[i],&tmpJ_P[i]);        //tmp_P_4x
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
    /*
    struct timeval tv_A,tv_B;
    gettimeofday(&tv_A,NULL);
    naf_length0 = w_naf(naf_binary[0],s[0],3);
    naf_length1 = w_naf(naf_binary[1],s[1],3);
    gettimeofday(&tv_B,NULL);
    printf("naf     : %.4f[ms]\n",timedifference_msec(tv_A,tv_B));
    */
    naf_length0 = w_naf(naf_binary[0],s[0],5);
    naf_length1 = w_naf(naf_binary[1],s[1],5);
    if(naf_length0<naf_length1)	naf_length=naf_length1;
    else naf_length=naf_length0;
    
    //naf_length=loop_length-1;
    int binary0[naf_length+1],binary1[naf_length+1];
    /*
    for(i=naf_length; i>=0; i--){
        if(naf_binary[0][i]==0)         binary0[i]=0;
        else if(naf_binary[0][i]==1)     binary0[i]=1;
        else if(naf_binary[0][i]==3)     binary0[i]=2;
        else if(naf_binary[0][i]==5)     binary0[i]=3;
        else if(naf_binary[0][i]==7)     binary0[i]=4;
        else if(naf_binary[0][i]==9)     binary0[i]=5;
        else if(naf_binary[0][i]==11)     binary0[i]=6;
        else if(naf_binary[0][i]==13)     binary0[i]=7;
        else if(naf_binary[0][i]==15)     binary0[i]=8;
        
        else if(naf_binary[0][i]==-1)     binary0[i]=9;
        else if(naf_binary[0][i]==-3)     binary0[i]=10;
        else if(naf_binary[0][i]==-5)     binary0[i]=11;
        else if(naf_binary[0][i]==-7)     binary0[i]=12;
        else if(naf_binary[0][i]==-9)     binary0[i]=13;
        else if(naf_binary[0][i]==-11)     binary0[i]=14;
        else if(naf_binary[0][i]==-13)     binary0[i]=15;
        else if(naf_binary[0][i]==-15)     binary0[i]=16;
        
        if(naf_binary[1][i]==0)         binary1[i]=0;
        else if(naf_binary[1][i]==1)     binary1[i]=1;
        else if(naf_binary[1][i]==3)     binary1[i]=2;
        else if(naf_binary[1][i]==5)     binary1[i]=3;
        else if(naf_binary[1][i]==7)     binary1[i]=4;
        else if(naf_binary[1][i]==9)     binary1[i]=5;
        else if(naf_binary[1][i]==11)     binary1[i]=6;
        else if(naf_binary[1][i]==13)     binary1[i]=7;
        else if(naf_binary[1][i]==15)     binary1[i]=8;
        
        else if(naf_binary[1][i]==-1)     binary1[i]=9;
        else if(naf_binary[1][i]==-3)     binary1[i]=10;
        else if(naf_binary[1][i]==-5)     binary1[i]=11;
        else if(naf_binary[1][i]==-7)     binary1[i]=12;
        else if(naf_binary[1][i]==-9)     binary1[i]=13;
        else if(naf_binary[1][i]==-11)     binary1[i]=14;
        else if(naf_binary[1][i]==-13)     binary1[i]=15;
        else if(naf_binary[1][i]==-15)     binary1[i]=16;
    }
    */
     for(i=naf_length; i>=0; i--){
        if(naf_binary[0][i]==0)         binary0[i]=0;
     	else if(naf_binary[0][i]>0)     	binary0[i]=(naf_binary[0][i]+1)>>1;
        else	binary0[i]=((17-(naf_binary[0][i]+16))>>1)+8;
        
        if(naf_binary[1][i]==0)         binary1[i]=0;
     	else if(naf_binary[1][i]>0)     	binary1[i]=(naf_binary[1][i]+1)>>1;
        else	binary1[i]=((17-(naf_binary[1][i]+16))>>1)+8;
    }
   if(naf_length0==naf_length1)	efp_eca_jacobian_lazy(&next_tmpJ_P,&table0[binary0[naf_length]],&table1[binary1[naf_length]]);
   else if(naf_length0<naf_length1)	efp_jacobian_set(&next_tmpJ_P,&table1[binary1[naf_length]]);
    else efp_jacobian_set(&next_tmpJ_P,&table0[binary0[naf_length]]);
    
    //scm
    for(i=naf_length-1; i>=0; i--){
        efp_ecd_jacobian_lazy(&next_tmpJ_P,&next_tmpJ_P);
        if(binary0[i]!=0)        efp_eca_mixture_lazy(&next_tmpJ_P,&next_tmpJ_P,&table0[binary0[i]]);
        if(binary1[i]!=0)        efp_eca_mixture_lazy(&next_tmpJ_P,&next_tmpJ_P,&table1[binary1[i]]);
    }
    
    efp_jacobian_to_affine(&next_tmp_P,&next_tmpJ_P);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_5naf_interleaving_mixture_lazy_montgomery(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    
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
	for(i=0;i<8;i++)     efp_mix_montgomery(&tmpJ_P[i],&tmpJ_P[i],&inv_table[i]);
	
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
    
    efp_jacobian_montgomery(&next_tmp_P,&next_tmpJ_P);
    efp_mod_montgomery(&next_tmp_P,&next_tmp_P);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g1_scm_2split_7naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    efp_t next_tmp_P,tmp_P;
    efp_jacobian_t next_tmpJ_P,tmpJ_P[32],tmpJ_P_neg[32],tmpJ_P_4x[32],tmpJ_P_4x_neg[32],tmpJ_2P;
    efp_init(&next_tmp_P);
    efp_init(&tmp_P);
    efp_jacobian_init(&tmpJ_2P);
    for(i=0;i<32;i++){
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
    efp_jacobian_t table0[65],table1[65];
    for(i=0; i<65; i++){
        efp_jacobian_init(&table0[i]);
        efp_jacobian_init(&table1[i]);
    }
    
    //set
    efp12_to_efp(&tmp_P,P);                    //tmp_P
    efp_affine_to_jacobian(&tmpJ_P[0],&tmp_P);
    efp_ecd_jacobian_lazy(&tmpJ_2P,&tmpJ_P[0]);
    for(i=1;i<32;i++){
	    efp_eca_jacobian_lazy(&tmpJ_P[i],&tmpJ_P[i-1],&tmpJ_2P);
    }
	
	fp_t point_table[32],inv_table[32];
	for(i=0;i<32;i++)    fp_set(&point_table[i],&tmpJ_P[i].z);
    
    fp_montgomery_trick(inv_table,point_table,32);
		
	for(i=0;i<32;i++)     efp_mix(&tmpJ_P[i],&tmpJ_P[i],&inv_table[i]);
	
	for(i=0;i<32;i++){
	    efp_jacobian_set_neg(&tmpJ_P_neg[i],&tmpJ_P[i]);            //tmp_P_neg
    	efp_jacobian_skew_frobenius_map_p2(&tmpJ_P_4x[i],&tmpJ_P[i]);        //tmp_P_4x
    	efp_jacobian_set_neg(&tmpJ_P_4x_neg[i],&tmpJ_P_4x[i]);        //tmp_P_4x_neg
    }
    //set table
    table0[0].infinity=1;                        //0
    table1[0].infinity=1;                        //0
    
    for(i=0;i<32;i++){
    efp_jacobian_set(&table0[i+1],&tmpJ_P[i]);                //[1]P
    efp_jacobian_set(&table0[i+33],&tmpJ_P_neg[i]);                //[-1]P
    efp_jacobian_set(&table1[i+1],&tmpJ_P_4x[i]);                //[1]P'
    efp_jacobian_set(&table1[i+33],&tmpJ_P_4x_neg[i]);                //[-1]P'
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
    
    naf_length0 = w_naf(naf_binary[0],s[0],7);
    naf_length1 = w_naf(naf_binary[1],s[1],7);
    if(naf_length0<naf_length1)	naf_length=naf_length1;
    else naf_length=naf_length0;
    
    //naf_length=loop_length-1;
    int binary0[naf_length+1],binary1[naf_length+1];
     for(i=naf_length; i>=0; i--){
        if(naf_binary[0][i]==0)         binary0[i]=0;
     	else if(naf_binary[0][i]>0)     	binary0[i]=(naf_binary[0][i]+1)>>1;
        else	binary0[i]=((65-(naf_binary[0][i]+64))>>1)+32;
        
        if(naf_binary[1][i]==0)         binary1[i]=0;
     	else if(naf_binary[1][i]>0)     	binary1[i]=(naf_binary[1][i]+1)>>1;
        else	binary1[i]=((65-(naf_binary[1][i]+64))>>1)+32;
    }
   if(naf_length0==naf_length1)	efp_eca_jacobian_lazy(&next_tmpJ_P,&table0[binary0[naf_length]],&table1[binary1[naf_length]]);
   else if(naf_length0<naf_length1)	efp_jacobian_set(&next_tmpJ_P,&table1[binary1[naf_length]]);
    else efp_jacobian_set(&next_tmpJ_P,&table0[binary0[naf_length]]);
    
    //scm
    for(i=naf_length-1; i>=0; i--){
        efp_ecd_jacobian_lazy(&next_tmpJ_P,&next_tmpJ_P);
        efp_eca_mixture_lazy(&next_tmpJ_P,&next_tmpJ_P,&table0[binary0[i]]);
        efp_eca_mixture_lazy(&next_tmpJ_P,&next_tmpJ_P,&table1[binary1[i]]);
    }
    
    efp_jacobian_to_affine(&next_tmp_P,&next_tmpJ_P);
    efp_to_efp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
}