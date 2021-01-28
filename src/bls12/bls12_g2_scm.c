#include <ELiPS/bls12_g2_scm.h>
void bls12_g2_scm_basic(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    efp2_t tmp_Q;
	efp2_init(&tmp_Q);
	
	efp12_to_efp2(&tmp_Q,Q);
	efp2_scm(&tmp_Q,&tmp_Q,scalar);
	efp2_to_efp12(ANS,&tmp_Q);
}
void bls12_g2_scm(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,j,length_s[4],loop_length;
    efp2_t next_twisted_Q,twisted_Q,twisted_2Q;
    efp2_jacobian_t next_twistedJ_Q,twistedJ_2Q;
    efp2_jacobian_t twistedJ_Q[8],twistedJ_Q_x[8],twistedJ_Q_2x[8],twistedJ_Q_3x[8];
    efp2_jacobian_t twistedJ_Q_neg[8],twistedJ_Q_x_neg[8],twistedJ_Q_2x_neg[8],twistedJ_Q_3x_neg[8];
    
    efp2_init(&twisted_Q);
    efp2_init(&twisted_2Q);
    efp2_init(&next_twisted_Q);
    efp2_jacobian_init(&next_twistedJ_Q);
    efp2_jacobian_init(&twistedJ_2Q);
    for(i=0;i<8;i++){
    efp2_jacobian_init(&twistedJ_Q[i]);
    efp2_jacobian_init(&twistedJ_Q_x[i]);
    efp2_jacobian_init(&twistedJ_Q_2x[i]);
    efp2_jacobian_init(&twistedJ_Q_3x[i]);
    efp2_jacobian_init(&twistedJ_Q_neg[i]);
    efp2_jacobian_init(&twistedJ_Q_x_neg[i]);
    efp2_jacobian_init(&twistedJ_Q_2x_neg[i]);
    efp2_jacobian_init(&twistedJ_Q_3x_neg[i]);
    }
    
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    efp2_jacobian_t table[4][17];
    for(i=0; i<17; i++){
        efp2_jacobian_init(&table[0][i]);
        efp2_jacobian_init(&table[1][i]);
        efp2_jacobian_init(&table[2][i]);
        efp2_jacobian_init(&table[3][i]);
    }
    
    //set
    efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    //TODO:lazy
    efp2_ecd(&twisted_2Q,&twisted_Q);
	efp2_to_montgomery(&twisted_Q,&twisted_Q);
	efp2_to_montgomery(&twisted_2Q,&twisted_2Q);
    efp2_affine_to_jacobian_montgomery(&twistedJ_Q[0],&twisted_Q);
    efp2_affine_to_jacobian_montgomery(&twistedJ_2Q,&twisted_2Q);
    for(i=1;i<8;i++){
	    efp2_eca_jacobian_lazy_montgomery(&twistedJ_Q[i],&twistedJ_Q[i-1],&twistedJ_2Q);
    }
    
	fp2_t point_table[8],inv_table[8];
	for(i=0;i<8;i++)    fp2_set(&point_table[i],&twistedJ_Q[i].z);
    fp2_montgomery_trick_montgomery(inv_table,point_table,8);
	for(i=0;i<8;i++)     efp2_mix_montgomery(&twistedJ_Q[i],&twistedJ_Q[i],&inv_table[i]);
	
	
	for(i=0;i<8;i++){
	    efp2_jacobian_set_neg(&twistedJ_Q_neg[i],&twistedJ_Q[i]);            //twisted_P_neg
    #ifdef X_PLUS
        efp2_jacobian_skew_frobenius_map_p1_montgomery(&twistedJ_Q_x[i],&twistedJ_Q[i]);        //twisted_Q_x
    #endif
    #ifdef X_MINUS
        efp2_jacobian_skew_frobenius_map_p1_montgomery(&twistedJ_Q_x[i],&twistedJ_Q[i]);        //twisted_Q_x
        efp2_jacobian_set_neg(&twistedJ_Q_x[i],&twistedJ_Q_x[i]);
    #endif
        efp2_jacobian_skew_frobenius_map_p2_montgomery(&twistedJ_Q_2x[i],&twistedJ_Q[i]);    //twisted_Q_2x
    #ifdef X_PLUS
        efp2_jacobian_skew_frobenius_map_p3_montgomery(&twistedJ_Q_3x[i],&twistedJ_Q[i]);    //twisted_Q_3x
    #endif
    #ifdef X_MINUS
        efp2_jacobian_skew_frobenius_map_p3_montgomery(&twistedJ_Q_3x[i],&twistedJ_Q[i]);    //twisted_Q_3x
        efp2_jacobian_set_neg(&twistedJ_Q_3x[i],&twistedJ_Q_3x[i]);
    #endif
        efp2_jacobian_set_neg(&twistedJ_Q_x_neg[i],&twistedJ_Q_x[i]);        //twisted_P_4x_neg
    	efp2_jacobian_set_neg(&twistedJ_Q_2x_neg[i],&twistedJ_Q_2x[i]);        //twisted_P_4x_neg
    	efp2_jacobian_set_neg(&twistedJ_Q_3x_neg[i],&twistedJ_Q_3x[i]);        //twisted_P_4x_neg
    }
    
    //set table
    table[0][0].infinity=1;                        //0
    table[1][0].infinity=1;                        //0
    table[2][0].infinity=1;                        //0
    table[3][0].infinity=1;                        //0
    
    for(i=0;i<8;i++){
    	efp2_jacobian_set(&table[0][i+1],&twistedJ_Q[i]);
    	efp2_jacobian_set(&table[0][i+9],&twistedJ_Q_neg[i]);
    	efp2_jacobian_set(&table[1][i+1],&twistedJ_Q_x[i]);
    	efp2_jacobian_set(&table[1][i+9],&twistedJ_Q_x_neg[i]);
    	efp2_jacobian_set(&table[2][i+1],&twistedJ_Q_2x[i]);
    	efp2_jacobian_set(&table[2][i+9],&twistedJ_Q_2x_neg[i]);
    	efp2_jacobian_set(&table[3][i+1],&twistedJ_Q_3x[i]); 
    	efp2_jacobian_set(&table[3][i+9],&twistedJ_Q_3x_neg[i]);
    }
    
    //set
    //s0,s1,s2,s3
    mpz_set(x_1,X_z);
    mpz_mul(x_2,x_1,x_1);
    mpz_tdiv_qr(B,A,scalar,x_2);
    mpz_tdiv_qr(s[1],s[0],A,x_1);
    mpz_tdiv_qr(s[3],s[2],B,x_1);
    
    //binary
    loop_length=0;
    for(i=0; i<4; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    
    //naf
    int naf_length[5];
    int naf_binary[4][loop_length+1];
    for(i=0; i<loop_length+1; i++){
        naf_binary[0][i]=0;
        naf_binary[1][i]=0;
        naf_binary[2][i]=0;
        naf_binary[3][i]=0;
    }
    int *naf_pointer[4];
    naf_pointer[0]=naf_binary[0];
    naf_pointer[1]=naf_binary[1];
    naf_pointer[2]=naf_binary[2];
    naf_pointer[3]=naf_binary[3];
    
    naf_length[1] = w_naf(naf_binary[0],s[0],5);
    naf_length[2] = w_naf(naf_binary[1],s[1],5);
    naf_length[3] = w_naf(naf_binary[2],s[2],5);
    naf_length[4] = w_naf(naf_binary[3],s[3],5);
    
    naf_length[0]=naf_length[1];
    for(i=2;i<5;i++){
    	if(naf_length[0]<naf_length[i])naf_length[0]=naf_length[i];
       }
    //naf_length=loop_length-1;
    int binary[4][naf_length[0]+1];
    
     for(i=naf_length[0]; i>=0; i--){
        if(naf_binary[0][i]==0)         binary[0][i]=0;
     	else if(naf_binary[0][i]>0)     	binary[0][i]=(naf_binary[0][i]+1)>>1;
        else	binary[0][i]=((17-(naf_binary[0][i]+16))>>1)+8;
        
        if(naf_binary[1][i]==0)         binary[1][i]=0;
     	else if(naf_binary[1][i]>0)     	binary[1][i]=(naf_binary[1][i]+1)>>1;
        else	binary[1][i]=((17-(naf_binary[1][i]+16))>>1)+8;
        
        if(naf_binary[2][i]==0)         binary[2][i]=0;
     	else if(naf_binary[2][i]>0)     	binary[2][i]=(naf_binary[2][i]+1)>>1;
        else	binary[2][i]=((17-(naf_binary[2][i]+16))>>1)+8;
        
        if(naf_binary[3][i]==0)         binary[3][i]=0;
     	else if(naf_binary[3][i]>0)     	binary[3][i]=(naf_binary[3][i]+1)>>1;
        else	binary[3][i]=((17-(naf_binary[3][i]+16))>>1)+8;
    }
    
	next_twistedJ_Q.infinity=1;
	for(i=1;i<5;i++){
		if(naf_length[0]==naf_length[i]){
			efp2_eca_jacobian_lazy_montgomery(&next_twistedJ_Q,&next_twistedJ_Q,&table[i-1][binary[i-1][naf_length[0]]]);
		}
	}
	
    //scm
    for(i=naf_length[0]-1; i>=0; i--){
        efp2_ecd_jacobian_lazy_montgomery(&next_twistedJ_Q,&next_twistedJ_Q);
        if(binary[0][i]!=0)        efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q,&next_twistedJ_Q,&table[0][binary[0][i]]);
        if(binary[1][i]!=0)        efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q,&next_twistedJ_Q,&table[1][binary[1][i]]);
        if(binary[2][i]!=0)        efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q,&next_twistedJ_Q,&table[2][binary[2][i]]);
        if(binary[3][i]!=0)        efp2_eca_mixture_lazy_montgomery(&next_twistedJ_Q,&next_twistedJ_Q,&table[3][binary[3][i]]);
    }
    
        
    efp2_jacobian_to_affine_montgomery(&next_twisted_Q,&next_twistedJ_Q);
    efp2_mod_montgomery(&next_twisted_Q,&next_twisted_Q);
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}