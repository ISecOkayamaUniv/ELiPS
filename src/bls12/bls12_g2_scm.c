#include <ELiPS/bls12_g2_scm.h>
void bls12_g2_scm_basic(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    efp2_t tmp_Q;
	efp2_init(&tmp_Q);
	
	efp12_to_efp2(&tmp_Q,Q);
	efp2_scm(&tmp_Q,&tmp_Q,scalar);
	efp2_to_efp12(ANS,&tmp_Q);
}
void bls12_g2_scm_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    efp2_t tmp_Q;
	efp2_init(&tmp_Q);
	
	efp12_to_efp2(&tmp_Q,Q);
	efp2_scm_lazy(&tmp_Q,&tmp_Q,scalar);
	efp2_to_efp12(ANS,&tmp_Q);
}

void bls12_g2_scm_2split(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    int i,length_s[2],loop_length;
	efp2_t next_twisted_Q,twisted_Q,skew_Q;
	efp2_init(&next_twisted_Q);
	efp2_init(&twisted_Q);
	efp2_init(&skew_Q);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	efp2_t table[4];
	for(i=0; i<4; i++){
		efp2_init(&table[i]);
	}
	
	//set
	efp12_to_efp2(&twisted_Q,Q);				//twisted_Q
	efp2_skew_frobenius_map_p1(&skew_Q,&twisted_Q);//skew_Q
	
	//set table
	table[0].infinity=1;						//00
	efp2_set(&table[1],&twisted_Q);			//01
	efp2_set(&table[2],&skew_Q);				//10
	efp2_eca(&table[3],&twisted_Q,&skew_Q);		//11
	
	//s0,s1
	mpz_sub_ui(buf,trace_z,1);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
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
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	efp2_set(&next_twisted_Q,&table[binary[0]]);
	
	//scm
	for(i=1; i<loop_length; i++){
		efp2_ecd(&next_twisted_Q,&next_twisted_Q);
		efp2_eca(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	efp2_to_efp12(ANS,&next_twisted_Q);
	ANS->infinity=next_twisted_Q.infinity;
	
	mpz_clear(buf);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
}

void bls12_g2_scm_2split_jsf(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    int i,length_s[2],loop_length;
	efp2_t next_tmp_Q,tmp_Q,tmp_Q_neg,skew_Q,skew_Q_neg;
	efp2_init(&next_tmp_Q);
	efp2_init(&tmp_Q);
	efp2_init(&tmp_Q_neg);
	efp2_init(&skew_Q);
	efp2_init(&skew_Q_neg);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	efp2_t table[9];
	for(i=0; i<9; i++){
		efp2_init(&table[i]);
	}
	
	//set
	efp12_to_efp2(&tmp_Q,Q);					//tmp_Q
	efp2_set_neg(&tmp_Q_neg,&tmp_Q);			//tmp_Q_neg
	efp2_skew_frobenius_map_p1(&skew_Q,&tmp_Q);		//skew_Q
	efp2_set_neg(&skew_Q_neg,&skew_Q);			//skew_Q_neg
	
	//set table
	table[0].infinity=1;						//00
	efp2_set(&table[1],&tmp_Q);				//01
	efp2_set(&table[2],&skew_Q);				//10
	efp2_eca(&table[3],&skew_Q,&tmp_Q);		//11
	efp2_set(&table[4],&tmp_Q_neg);			//0-1
	efp2_set(&table[5],&skew_Q_neg);			//-10
	efp2_eca(&table[6],&skew_Q_neg,&tmp_Q_neg);	//-1-1
	efp2_eca(&table[7],&skew_Q,&tmp_Q_neg);		//1-1
	efp2_eca(&table[8],&skew_Q_neg,&tmp_Q);		//-11
	
	//s0,s1
	mpz_sub_ui(buf,trace_z,1);
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
		if(jsf_binary[1][i]==0 && jsf_binary[0][i]==0) 		binary[i]=0;
		else if(jsf_binary[1][i]==0 && jsf_binary[0][i]==1) 	binary[i]=1;
		else if(jsf_binary[1][i]==1 && jsf_binary[0][i]==0) 	binary[i]=2;
		else if(jsf_binary[1][i]==1 && jsf_binary[0][i]==1)	binary[i]=3;
		else if(jsf_binary[1][i]==0 && jsf_binary[0][i]==-1)	binary[i]=4;
		else if(jsf_binary[1][i]==-1 && jsf_binary[0][i]==0)	binary[i]=5;
		else if(jsf_binary[1][i]==-1 && jsf_binary[0][i]==-1)	binary[i]=6;
		else if(jsf_binary[1][i]==1 && jsf_binary[0][i]==-1)	binary[i]=7;
		else if(jsf_binary[1][i]==-1 && jsf_binary[0][i]==1)	binary[i]=8;
	}
	efp2_set(&next_tmp_Q,&table[binary[jsf_length]]);
	//scm
	for(i=jsf_length-1; i>=0; i--){
		efp2_ecd(&next_tmp_Q,&next_tmp_Q);
		efp2_eca(&next_tmp_Q,&next_tmp_Q,&table[binary[i]]);
	}
	efp2_to_efp12(ANS,&next_tmp_Q);

	mpz_clear(buf);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
}

void bls12_g2_scm_4split(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,length_s[4],loop_length;
    efp2_t next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    efp2_init(&next_twisted_Q);
    efp2_init(&twisted_Q);
    efp2_init(&twisted_Q_x);
    efp2_init(&twisted_Q_2x);
    efp2_init(&twisted_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    efp2_t table[16];
    for(i=0; i<16; i++){
        efp2_init(&table[i]);
    }
    
    //set twisted_Q
    efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    efp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    efp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    efp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    
    //set table
    table[0].infinity=1;                            //0000
    efp2_set(&table[1],&twisted_Q);                //0001
    efp2_set(&table[2],&twisted_Q_x);                //0010
    efp2_eca(&table[3],&twisted_Q_x,&twisted_Q);        //0011
    efp2_set(&table[4],&twisted_Q_2x);                //0100
    efp2_eca(&table[5],&twisted_Q_2x,&twisted_Q);        //0101
    efp2_eca(&table[6],&twisted_Q_2x,&twisted_Q_x);        //0110
    efp2_eca(&table[7],&table[6],&twisted_Q);            //0111
    efp2_set(&table[8],&twisted_Q_3x);                //1000
    efp2_eca(&table[9],&twisted_Q_3x,&twisted_Q);        //1001
    efp2_eca(&table[10],&twisted_Q_3x,&twisted_Q_x);    //1010
    efp2_eca(&table[11],&twisted_Q_3x,&table[3]);        //1011
    efp2_eca(&table[12],&twisted_Q_3x,&twisted_Q_2x);    //1100
    efp2_eca(&table[13],&table[12],&twisted_Q);        //1101
    efp2_eca(&table[14],&table[12],&twisted_Q_x);        //1110
    efp2_eca(&table[15],&table[14],&twisted_Q);        //1111
    
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
        //printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //printf("\n");
    //set binary
    char binary_s[4][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<4; i++){
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
        sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
        binary[i]=strtol(str,&e,2);
    }
    
    efp2_set(&next_twisted_Q,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        efp2_ecd(&next_twisted_Q,&next_twisted_Q);
        efp2_eca(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
    }
    
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;
    
    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}

void bls12_g2_scm_4split_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,length_s[4],loop_length;
    efp2_t next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    efp2_init(&next_twisted_Q);
    efp2_init(&twisted_Q);
    efp2_init(&twisted_Q_x);
    efp2_init(&twisted_Q_2x);
    efp2_init(&twisted_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    efp2_t table[16];
    for(i=0; i<16; i++){
        efp2_init(&table[i]);
    }
    
    //set twisted_Q
    efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    efp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    efp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    efp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    
    //set table
    table[0].infinity=1;                            //0000
    efp2_set(&table[1],&twisted_Q);                //0001
    efp2_set(&table[2],&twisted_Q_x);                //0010
    efp2_eca_lazy(&table[3],&twisted_Q_x,&twisted_Q);        //0011
    efp2_set(&table[4],&twisted_Q_2x);                //0100
    efp2_eca_lazy(&table[5],&twisted_Q_2x,&twisted_Q);        //0101
    efp2_eca_lazy(&table[6],&twisted_Q_2x,&twisted_Q_x);        //0110
    efp2_eca_lazy(&table[7],&table[6],&twisted_Q);            //0111
    efp2_set(&table[8],&twisted_Q_3x);                //1000
    efp2_eca_lazy(&table[9],&twisted_Q_3x,&twisted_Q);        //1001
    efp2_eca_lazy(&table[10],&twisted_Q_3x,&twisted_Q_x);    //1010
    efp2_eca_lazy(&table[11],&twisted_Q_3x,&table[3]);        //1011
    efp2_eca_lazy(&table[12],&twisted_Q_3x,&twisted_Q_2x);    //1100
    efp2_eca_lazy(&table[13],&table[12],&twisted_Q);        //1101
    efp2_eca_lazy(&table[14],&table[12],&twisted_Q_x);        //1110
    efp2_eca_lazy(&table[15],&table[14],&twisted_Q);        //1111
    
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
        //printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //printf("\n");
    //set binary
    char binary_s[4][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<4; i++){
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
        sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
        binary[i]=strtol(str,&e,2);
    }
    
    efp2_set(&next_twisted_Q,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        efp2_ecd_lazy(&next_twisted_Q,&next_twisted_Q);
        efp2_eca_lazy(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
    }
    
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;
    
    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g2_scm_4split_jacobian_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,length_s[4],loop_length;
    efp2_t next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    efp2_jacobian_t next_twistedZ_Q,twistedZ_Q,twistedZ_Q_x,twistedZ_Q_2x,twistedZ_Q_3x;
    efp2_init(&next_twisted_Q);
    efp2_init(&twisted_Q);
    efp2_init(&twisted_Q_x);
    efp2_init(&twisted_Q_2x);
    efp2_init(&twisted_Q_3x);
    efp2_jacobian_init(&next_twistedZ_Q);
    efp2_jacobian_init(&twistedZ_Q);
    efp2_jacobian_init(&twistedZ_Q_x);
    efp2_jacobian_init(&twistedZ_Q_2x);
    efp2_jacobian_init(&twistedZ_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    efp2_jacobian_t table[16];
    for(i=0; i<16; i++){
        efp2_jacobian_init(&table[i]);
    }

    
    //set twisted_Q
    efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    efp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    efp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    efp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    //set jacobian
    efp2_affine_to_jacobian(&twistedZ_Q,&twisted_Q);
    efp2_affine_to_jacobian(&twistedZ_Q_x,&twisted_Q_x);
    efp2_affine_to_jacobian(&twistedZ_Q_2x,&twisted_Q_2x);
    efp2_affine_to_jacobian(&twistedZ_Q_3x,&twisted_Q_3x);
    
    //set table
    table[0].infinity=1;                            //0000
    efp2_affine_to_jacobian(&table[1],&twisted_Q);                //0001
    efp2_affine_to_jacobian(&table[2],&twisted_Q_x);                //0010
    efp2_eca_jacobian_lazy(&table[3],&twistedZ_Q_x,&twistedZ_Q);        //0011
    efp2_affine_to_jacobian(&table[4],&twisted_Q_2x);                //0100
    efp2_eca_jacobian_lazy(&table[5],&twistedZ_Q_2x,&twistedZ_Q);        //0101
    efp2_eca_jacobian_lazy(&table[6],&twistedZ_Q_2x,&twistedZ_Q_x);        //0110
    efp2_eca_jacobian_lazy(&table[7],&table[6],&twistedZ_Q);            //0111
    efp2_affine_to_jacobian(&table[8],&twisted_Q_3x);                //1000
    efp2_eca_jacobian_lazy(&table[9],&twistedZ_Q_3x,&twistedZ_Q);        //1001
    efp2_eca_jacobian_lazy(&table[10],&twistedZ_Q_3x,&twistedZ_Q_x);    //1010
    efp2_eca_jacobian_lazy(&table[11],&twistedZ_Q_3x,&table[3]);        //1011
    efp2_eca_jacobian_lazy(&table[12],&twistedZ_Q_3x,&twistedZ_Q_2x);    //1100
    efp2_eca_jacobian_lazy(&table[13],&table[12],&twistedZ_Q);        //1101
    efp2_eca_jacobian_lazy(&table[14],&table[12],&twistedZ_Q_x);        //1110
    efp2_eca_jacobian_lazy(&table[15],&table[14],&twistedZ_Q);        //1111
    

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
        //printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //printf("\n");
    //set binary
    char binary_s[4][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<4; i++){
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
        sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
        binary[i]=strtol(str,&e,2);
    }
    
    efp2_jacobian_set(&next_twistedZ_Q,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        efp2_ecd_jacobian_lazy(&next_twistedZ_Q,&next_twistedZ_Q);
        efp2_eca_jacobian_lazy(&next_twistedZ_Q,&next_twistedZ_Q,&table[binary[i]]);
    }
    efp2_jacobian_to_affine(&next_twisted_Q,&next_twistedZ_Q);
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g2_scm_4split_mixture_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,j,length_s[4],loop_length;
    efp2_t next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    efp2_jacobian_t next_twistedZ_Q,twistedZ_Q,twistedZ_Q_x,twistedZ_Q_2x,twistedZ_Q_3x;
    efp2_init(&next_twisted_Q);
    efp2_init(&twisted_Q);
    efp2_init(&twisted_Q_x);
    efp2_init(&twisted_Q_2x);
    efp2_init(&twisted_Q_3x);
    efp2_jacobian_init(&next_twistedZ_Q);
    efp2_jacobian_init(&twistedZ_Q);
    efp2_jacobian_init(&twistedZ_Q_x);
    efp2_jacobian_init(&twistedZ_Q_2x);
    efp2_jacobian_init(&twistedZ_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    efp2_jacobian_t table[16];
    for(i=0; i<16; i++){
        efp2_jacobian_init(&table[i]);
    }

    
    //set twisted_Q
    efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    efp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    efp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    efp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    //set jacobian
    efp2_affine_to_jacobian(&twistedZ_Q,&twisted_Q);
    efp2_affine_to_jacobian(&twistedZ_Q_x,&twisted_Q_x);
    efp2_affine_to_jacobian(&twistedZ_Q_2x,&twisted_Q_2x);
    efp2_affine_to_jacobian(&twistedZ_Q_3x,&twisted_Q_3x);
    
    //set table
    table[0].infinity=1;                            //0000
    efp2_affine_to_jacobian(&table[1],&twisted_Q);                //0001
    efp2_affine_to_jacobian(&table[2],&twisted_Q_x);                //0010
    efp2_eca_jacobian_lazy(&table[3],&twistedZ_Q_x,&twistedZ_Q);        //0011
    efp2_affine_to_jacobian(&table[4],&twisted_Q_2x);                //0100
    efp2_eca_jacobian_lazy(&table[5],&twistedZ_Q_2x,&twistedZ_Q);        //0101
    efp2_eca_jacobian_lazy(&table[6],&twistedZ_Q_2x,&twistedZ_Q_x);        //0110
    efp2_eca_jacobian_lazy(&table[7],&table[6],&twistedZ_Q);            //0111
    efp2_affine_to_jacobian(&table[8],&twisted_Q_3x);                //1000
    efp2_eca_jacobian_lazy(&table[9],&twistedZ_Q_3x,&twistedZ_Q);        //1001
    efp2_eca_jacobian_lazy(&table[10],&twistedZ_Q_3x,&twistedZ_Q_x);    //1010
    efp2_eca_jacobian_lazy(&table[11],&twistedZ_Q_3x,&table[3]);        //1011
    efp2_eca_jacobian_lazy(&table[12],&twistedZ_Q_3x,&twistedZ_Q_2x);    //1100
    efp2_eca_jacobian_lazy(&table[13],&table[12],&twistedZ_Q);        //1101
    efp2_eca_jacobian_lazy(&table[14],&table[12],&twistedZ_Q_x);        //1110
    efp2_eca_jacobian_lazy(&table[15],&table[14],&twistedZ_Q);        //1111
    
    fp2_t mul_table[11],chk;
    fp2_t ans_table[11];
    fp2_t inv_table[11];
    fp2_t inv;
    
    fp2_t all;
    
    j=0;
    for(i=0;i<16;i++){
    	if( i != 0 && i != 1 && i != 2 && i != 4 && i != 8){
    	fp2_set(&mul_table[j],&table[i].z);
    	j++;
    	}
    }
    fp2_set(&ans_table[0],&mul_table[0]);
    for(i=1;i<11;i++){
    fp2_mul(&ans_table[i],&ans_table[i-1],&mul_table[i]);
    }
    fp2_inv(&all,&ans_table[10]);
    for(i=10;i>0;i--){
    fp2_mul(&inv_table[i],&all,&ans_table[i-1]);	
    fp2_mul(&all,&all,&mul_table[i]);
    }
    fp2_set(&inv_table[0],&all);
    
	for(i=0;i<11;i++){
    fp2_mul(&chk,&inv_table[i],&mul_table[i]);
    }
	j=0;
    for(i=0;i<16;i++){
    	if( i != 0 && i != 1 && i != 2 && i != 4 && i != 8){
    	efp2_mix(&table[i],&table[i],&inv_table[j]);
    	j++;
    	}
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
        //printf("G2 length_s%d:%d\n",i,length_s[i]);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //printf("\n");
    //set binary
    char binary_s[4][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<4; i++){
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
        sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
        binary[i]=strtol(str,&e,2);
    }
    
    efp2_jacobian_set(&next_twistedZ_Q,&table[binary[0]]);
    
    efp2_jacobian_t out;
    efp2_jacobian_init(&out);
    //scm
    for(i=1; i<loop_length; i++){
        efp2_ecd_jacobian_lazy(&next_twistedZ_Q,&next_twistedZ_Q);        efp2_eca_mixture_lazy(&next_twistedZ_Q,&next_twistedZ_Q,&table[binary[i]]);
    }
    efp2_jacobian_to_affine(&next_twisted_Q,&next_twistedZ_Q);
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
/*
void bls12_g2_scm_4split_2naf_shamir_mixture_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,j,length_s[4],loop_length;
    efp2_t next_twisted_Q,twisted_Q;
    efp2_jacobian_t next_twistedJ_Q,twistedJ_2Q;
    efp2_jacobian_t twistedJ_Q,twistedJ_Q_x,twistedJ_Q_2x,twistedJ_Q_3x;
    efp2_jacobian_t twistedJ_Q_neg,twistedJ_Q_x_neg,twistedJ_Q_2x_neg,twistedJ_Q_3x_neg;
    
    efp2_init(&twisted_Q);
    efp2_init(&next_twisted_Q);
    efp2_jacobian_init(&next_twistedJ_Q);
    efp2_jacobian_init(&twistedJ_2Q);
    efp2_jacobian_init(&twistedJ_Q[i]);
    efp2_jacobian_init(&twistedJ_Q_x[i]);
    efp2_jacobian_init(&twistedJ_Q_2x[i]);
    efp2_jacobian_init(&twistedJ_Q_3x[i]);
    efp2_jacobian_init(&twistedJ_Q_neg[i]);
    efp2_jacobian_init(&twistedJ_Q_x_neg[i]);
    efp2_jacobian_init(&twistedJ_Q_2x_neg[i]);
    efp2_jacobian_init(&twistedJ_Q_3x_neg[i]);
    
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    efp2_jacobian_t table[81];
    for(i=0; i<81; i++){
        efp2_jacobian_init(&table[i]);
    }
    
    //set
    efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    efp2_affine_to_jacobian(&twistedJ_Q[0],&twisted_Q);
    efp2_ecd_jacobian_lazy(&twistedJ_2Q,&twistedJ_Q[0]);
    for(i=1;i<2;i++){
	    efp2_eca_jacobian_lazy(&twistedJ_Q[i],&twistedJ_Q[i-1],&twistedJ_2Q);
    }
	
	fp2_t point_table,inv_table;
	fp2_set(&point_table,&twistedJ_Q.z);
    fp2_inv_lazy(inv_table,point_table);
	efp2_mix(&twistedJ_Q,&twistedJ_Q,&inv_table);
	
	efp2_jacobian_set_neg(&twistedJ_Q_neg,&twistedJ_Q);            //twisted_P_neg
    efp2_jacobian_skew_frobenius_map_p1(&twistedJ_Q_x,&twistedJ_Q);        //twisted_Q_x
    efp2_jacobian_skew_frobenius_map_p2(&twistedJ_Q_2x,&twistedJ_Q);    //twisted_Q_2x
    efp2_jacobian_skew_frobenius_map_p3(&twistedJ_Q_3x,&twistedJ_Q);    //twisted_Q_3x
    efp2_jacobian_set_neg(&twistedJ_Q_x_neg,&twistedJ_Q_x);        //twisted_P_4x_neg
   	efp2_jacobian_set_neg(&twistedJ_Q_2x_neg,&twistedJ_Q_2x);        //twisted_P_4x_neg
    efp2_jacobian_set_neg(&twistedJ_Q_3x_neg,&twistedJ_Q_3x);        //twisted_P_4x_neg
    
    //set table
    table[0].infinity=1;                            //0000
    efp2_jacobian_set(&table[1],&twistedJ_Q);                //0001
    efp2_jacobian_set(&table[2],&twistedJ_Q_neg);                //000-1
    
    efp2_jacobian_set(&table[3],&twisted_Q_x);                //0010
    efp2_eca_mixture_lazy(&table[4],&twistedJ_Q_x,&twistedJ_Q);        //0011
    efp2_eca_mixture_lazy(&table[5],&twistedJ_Q_x,&twistedJ_Q);        //001-1
    
    efp2_jacobian_set(&table[6],&twisted_Q_x_neg);                //00-10
    efp2_set_neg(&table[7],&table[5]);        //00-11
    efp2_set_neg(&table[8],&table[4]);        //00-1-1
    
    efp2_jacobian_set(&table[9],&twisted_Q_2x);                //0100
    efp2_eca_mixture_lazy(&table[10],&twistedZ_Q_2x,&twistedZ_Q);        //0101
    efp2_eca_mixture_lazy(&table[11],&twistedZ_Q_2x,&twistedZ_Q_neg);                 //010-1
    
    efp2_eca_mixture_lazy(&table[12],&twistedZ_Q_2x,&twistedZ_Q_x);        //0110
    efp2_eca_mixture_lazy(&table[13],&table[12],&twistedZ_Q);            //0111
    efp2_eca_mixture_lazy(&table[14],&table[12],&twistedZ_Q_neg);            //011-1
    
    efp2_eca_mixture_lazy(&table[15],&twistedZ_Q_2x,&twistedZ_Q_x_neg);        //01-10
    efp2_eca_mixture_lazy(&table[16],&table[15],&twistedZ_Q);        //01-11
    efp2_eca_mixture_lazy(&table[17],&table[15],&twistedZ_Q_neg);        //01-1-1
    
    efp2_jacobian_set(&table[18],&twisted_Q_2x_neg);                //0-100
    efp2_set_neg(&table[19],&tsble[11]);        //0-101
    efp2_set_neg(&table[20],&tsble[12]);        //0-10-1
    
    efp2_set_neg(&table[21],&ttable[15]);        //0-110
    efp2_set_neg(&table[22],&ttable[17]);            //0-111
    efp2_set_neg(&table[23],&ttable[16]);            //0-11-1
    
    efp2_set_neg(&table[24],&ttable[12]);       //0-1-10
    efp2_set_neg(&table[25],&ttable[14]);        //0-1-11
    efp2_set_neg(&table[26],&ttable[13]);        //0-1-1-1
    
    
    efp2_jacobian_set(&table[1],&twisted_Q);                            //1000
    efp2_jacobian_set(&table[1],&twisted_Q);                //1001
    efp2_jacobian_set(&table[2],&twisted_Q_neg);                //100-1
    
    efp2_jacobian_set(&table[3],&twisted_Q_x);                //1010
    efp2_eca_mixture_lazy(&table[4],&twistedZ_Q_x,&twistedZ_Q);        //1011
    efp2_eca_mixture_lazy(&table[5],&twistedZ_Q_x,&twistedZ_Q);        //101-1
    
    efp2_jacobian_set(&table[6],&twisted_Q_x_neg);                //10-10
    efp2_set_neg(&table[7],&table[5]);        //10-11
    efp2_set_neg(&table[8],&table[4]);        //10-1-1
    
    efp2_jacobian_set(&table[9],&twisted_Q_2x);                //1100
    efp2_eca_mixture_lazy(&table[10],&twistedZ_Q_2x,&twistedZ_Q);        //1101
    efp2_eca_mixture_lazy(&table[11],&twistedZ_Q_2x,&twistedZ_Q_neg);                 //110-1
    
    efp2_eca_mixture_lazy(&table[12],&twistedZ_Q_2x,&twistedZ_Q_x);        //1110
    efp2_eca_mixture_lazy(&table[13],&table[12],&twistedZ_Q);            //1111
    efp2_eca_mixture_lazy(&table[14],&table[12],&twistedZ_Q_neg);            //111-1
    
    efp2_eca_mixture_lazy(&table[15],&twistedZ_Q_2x,&twistedZ_Q_x_neg);        //11-10
    efp2_eca_mixture_lazy(&table[16],&table[15],&twistedZ_Q);        //11-11
    efp2_eca_mixture_lazy(&table[17],&table[15],&twistedZ_Q_neg);        //11-1-1
    
    efp2_jacobian_set(&table[18],&twisted_Q_2x_neg);                //1-100
    efp2_set_neg(&table[19],&tsble[11]);        //1-101
    efp2_set_neg(&table[20],&tsble[12]);        //1-10-1
    
    efp2_set_neg(&table[21],&ttable[15]);        //1-110
    efp2_set_neg(&table[22],&ttable[17]);            //1-111
    efp2_set_neg(&table[23],&ttable[16]);            //1-11-1
    
    efp2_set_neg(&table[24],&ttable[12]);       //1-1-10
    efp2_set_neg(&table[25],&ttable[14]);        //1-1-11
    efp2_set_neg(&table[26],&ttable[13]);        //1-1-1-1
    
    
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
    
    naf_length[1] = w_naf(naf_binary[0],s[0],3);
    naf_length[2] = w_naf(naf_binary[1],s[1],3);
    naf_length[3] = w_naf(naf_binary[2],s[2],3);
    naf_length[4] = w_naf(naf_binary[3],s[3],3);
    
    naf_length[0]=naf_length[1];
    for(i=2;i<5;i++){
    	if(naf_length[0]<naf_length[i])naf_length[0]=naf_length[i];
       }
    //naf_length=loop_length-1;
    int binary[4][naf_length[0]+1];
    
     for(i=naf_length[0]; i>=0; i--){
     	binary[i]=0;
        if(naf_binary[0][i]==0)         binary[i]+=0;
        else if(naf_binary[0][i]==1)         binary[i]+=1;
        else if(naf_binary[0][i]==-1)         binary[i]+=2;
        
        if(naf_binary[1][i]==0)         binary[i]+=0;
        else if(naf_binary[1][i]==1)         binary[i]+=3;
        else if(naf_binary[1][i]==-1)         binary[i]+=6;
        
        if(naf_binary[2][i]==0)         binary[i]+=0;
        else if(naf_binary[2][i]==1)         binary[i]+=9
        else if(naf_binary[2][i]==-1)         binary[i]+=18;
        
        if(naf_binary[3][i]==0)         binary[i]+=0;
        else if(naf_binary[3][i]==1)         binary[i]+=27;
        else if(naf_binary[3][i]==-1)         binary[i]+=54;
    }
	efp2_jacobian_set(&next_twistedJ_Q,&table[binary[naf_length[0]]]);
    //scm
    for(i=naf_length[0]-1; i>=0; i--){
        efp2_ecd_jacobian_lazy(&next_twistedJ_Q,&next_twistedJ_Q);
        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[binary[i]]);
    }
        
    efp2_jacobian_to_affine(&next_twisted_Q,&next_twistedJ_Q);
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
*/
void bls12_g2_scm_4split_3naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,j,length_s[4],loop_length;
    efp2_t next_twisted_Q,twisted_Q;
    efp2_jacobian_t next_twistedJ_Q,twistedJ_2Q;
    efp2_jacobian_t twistedJ_Q[2],twistedJ_Q_x[2],twistedJ_Q_2x[2],twistedJ_Q_3x[2];
    efp2_jacobian_t twistedJ_Q_neg[2],twistedJ_Q_x_neg[2],twistedJ_Q_2x_neg[2],twistedJ_Q_3x_neg[2];
    
    efp2_init(&twisted_Q);
    efp2_init(&next_twisted_Q);
    efp2_jacobian_init(&next_twistedJ_Q);
    efp2_jacobian_init(&twistedJ_2Q);
    for(i=0;i<2;i++){
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
    efp2_jacobian_t table[4][5];
    for(i=0; i<5; i++){
        efp2_jacobian_init(&table[0][i]);
        efp2_jacobian_init(&table[1][i]);
        efp2_jacobian_init(&table[2][i]);
        efp2_jacobian_init(&table[3][i]);
    }
    
    //set
    efp12_to_efp2(&twisted_Q,Q);                    //twisted_Q
    efp2_affine_to_jacobian(&twistedJ_Q[0],&twisted_Q);
    efp2_ecd_jacobian_lazy(&twistedJ_2Q,&twistedJ_Q[0]);
    for(i=1;i<2;i++){
	    efp2_eca_jacobian_lazy(&twistedJ_Q[i],&twistedJ_Q[i-1],&twistedJ_2Q);
    }
	
	fp2_t point_table[2],inv_table[2];
	for(i=0;i<2;i++)    fp2_set(&point_table[i],&twistedJ_Q[i].z);
    fp2_montgomery_trick(inv_table,point_table,2);
	for(i=0;i<2;i++)     efp2_mix(&twistedJ_Q[i],&twistedJ_Q[i],&inv_table[i]);
	
	for(i=0;i<2;i++){
	    efp2_jacobian_set_neg(&twistedJ_Q_neg[i],&twistedJ_Q[i]);            //twisted_P_neg
    	efp2_jacobian_skew_frobenius_map_p1(&twistedJ_Q_x[i],&twistedJ_Q[i]);        //twisted_Q_x
    	efp2_jacobian_skew_frobenius_map_p2(&twistedJ_Q_2x[i],&twistedJ_Q[i]);    //twisted_Q_2x
    	efp2_jacobian_skew_frobenius_map_p3(&twistedJ_Q_3x[i],&twistedJ_Q[i]);    //twisted_Q_3x
    	efp2_jacobian_set_neg(&twistedJ_Q_x_neg[i],&twistedJ_Q_x[i]);        //twisted_P_4x_neg
    	efp2_jacobian_set_neg(&twistedJ_Q_2x_neg[i],&twistedJ_Q_2x[i]);        //twisted_P_4x_neg
    	efp2_jacobian_set_neg(&twistedJ_Q_3x_neg[i],&twistedJ_Q_3x[i]);        //twisted_P_4x_neg
    }
    
    //set table
    table[0][0].infinity=1;                        //0
    table[1][0].infinity=1;                        //0
    table[2][0].infinity=1;                        //0
    table[3][0].infinity=1;                        //0
    
    for(i=0;i<2;i++){
    	efp2_jacobian_set(&table[0][i+1],&twistedJ_Q[i]);
    	efp2_jacobian_set(&table[0][i+3],&twistedJ_Q_neg[i]);
    	efp2_jacobian_set(&table[1][i+1],&twistedJ_Q_x[i]);
    	efp2_jacobian_set(&table[1][i+3],&twistedJ_Q_x_neg[i]);
    	efp2_jacobian_set(&table[2][i+1],&twistedJ_Q_2x[i]);
    	efp2_jacobian_set(&table[2][i+3],&twistedJ_Q_2x_neg[i]);
    	efp2_jacobian_set(&table[3][i+1],&twistedJ_Q_3x[i]); 
    	efp2_jacobian_set(&table[3][i+3],&twistedJ_Q_3x_neg[i]);
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
    
    naf_length[1] = w_naf(naf_binary[0],s[0],3);
    naf_length[2] = w_naf(naf_binary[1],s[1],3);
    naf_length[3] = w_naf(naf_binary[2],s[2],3);
    naf_length[4] = w_naf(naf_binary[3],s[3],3);
    
    naf_length[0]=naf_length[1];
    for(i=2;i<5;i++){
    	if(naf_length[0]<naf_length[i])naf_length[0]=naf_length[i];
       }
    //naf_length=loop_length-1;
    int binary[4][naf_length[0]+1];
    
     for(i=naf_length[0]; i>=0; i--){
        if(naf_binary[0][i]==0)         binary[0][i]=0;
     	else if(naf_binary[0][i]>0)     	binary[0][i]=(naf_binary[0][i]+1)>>1;
        else	binary[0][i]=((5-(naf_binary[0][i]+4))>>1)+2;
        
        if(naf_binary[1][i]==0)         binary[1][i]=0;
     	else if(naf_binary[1][i]>0)     	binary[1][i]=(naf_binary[1][i]+1)>>1;
        else	binary[1][i]=((5-(naf_binary[1][i]+4))>>1)+2;
        
        if(naf_binary[2][i]==0)         binary[2][i]=0;
     	else if(naf_binary[2][i]>0)     	binary[2][i]=(naf_binary[2][i]+1)>>1;
        else	binary[2][i]=((5-(naf_binary[2][i]+4))>>1)+2;
        
        if(naf_binary[3][i]==0)         binary[3][i]=0;
     	else if(naf_binary[3][i]>0)     	binary[3][i]=(naf_binary[3][i]+1)>>1;
        else	binary[3][i]=((5-(naf_binary[3][i]+4))>>1)+2;
    }
    
	next_twistedJ_Q.infinity=1;
	for(i=1;i<5;i++){
		if(naf_length[0]==naf_length[i]){
			efp2_eca_jacobian_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[i-1][binary[i-1][naf_length[0]]]);
		}
	}
	
    //scm
    for(i=naf_length[0]-1; i>=0; i--){
        efp2_ecd_jacobian_lazy(&next_twistedJ_Q,&next_twistedJ_Q);
        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[0][binary[0][i]]);
        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[1][binary[1][i]]);
        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[2][binary[2][i]]);
        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[3][binary[3][i]]);
    }
    
        
    efp2_jacobian_to_affine(&next_twisted_Q,&next_twistedJ_Q);
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g2_scm_4split_5naf_interleaving_mixture_lazy(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
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
    efp2_ecd_lazy(&twisted_2Q,&twisted_Q);
    efp2_affine_to_jacobian(&twistedJ_Q[0],&twisted_Q);
    efp2_affine_to_jacobian(&twistedJ_2Q,&twisted_2Q);
    for(i=1;i<8;i++){
	    efp2_eca_jacobian_lazy(&twistedJ_Q[i],&twistedJ_Q[i-1],&twistedJ_2Q);
    }
    
	fp2_t point_table[8],inv_table[8];
	for(i=0;i<8;i++)    fp2_set(&point_table[i],&twistedJ_Q[i].z);
    fp2_montgomery_trick(inv_table,point_table,8);
	for(i=0;i<8;i++)     efp2_mix(&twistedJ_Q[i],&twistedJ_Q[i],&inv_table[i]);
	
	for(i=0;i<8;i++){
	    efp2_jacobian_set_neg(&twistedJ_Q_neg[i],&twistedJ_Q[i]);            //twisted_P_neg
    	efp2_jacobian_skew_frobenius_map_p1(&twistedJ_Q_x[i],&twistedJ_Q[i]);        //twisted_Q_x
    	efp2_jacobian_skew_frobenius_map_p2(&twistedJ_Q_2x[i],&twistedJ_Q[i]);    //twisted_Q_2x
    	efp2_jacobian_skew_frobenius_map_p3(&twistedJ_Q_3x[i],&twistedJ_Q[i]);    //twisted_Q_3x
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
			efp2_eca_jacobian_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[i-1][binary[i-1][naf_length[0]]]);
		}
	}
	
    //scm
    for(i=naf_length[0]-1; i>=0; i--){
        efp2_ecd_jacobian_lazy(&next_twistedJ_Q,&next_twistedJ_Q);
        if(binary[0][i]!=0)        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[0][binary[0][i]]);
        if(binary[1][i]!=0)        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[1][binary[1][i]]);
        if(binary[2][i]!=0)        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[2][binary[2][i]]);
        if(binary[3][i]!=0)        efp2_eca_mixture_lazy(&next_twistedJ_Q,&next_twistedJ_Q,&table[3][binary[3][i]]);
    }
    
        
    efp2_jacobian_to_affine(&next_twisted_Q,&next_twistedJ_Q);
    efp2_to_efp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g2_scm_4split_5naf_interleaving_mixture_lazy_montgomery(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
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
    efp2_ecd_lazy(&twisted_2Q,&twisted_Q);
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
    	efp2_jacobian_skew_frobenius_map_p1_montgomery(&twistedJ_Q_x[i],&twistedJ_Q[i]);        //twisted_Q_x
    	efp2_jacobian_skew_frobenius_map_p2_montgomery(&twistedJ_Q_2x[i],&twistedJ_Q[i]);    //twisted_Q_2x
    	efp2_jacobian_skew_frobenius_map_p3_montgomery(&twistedJ_Q_3x[i],&twistedJ_Q[i]);    //twisted_Q_3x
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