#include <ELiPS/G1_SCM.h>
void EFp12_G1_SCM_plain(EFp12 *ANS,EFp12 *P,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        EFp tmp_P;
	EFp_init(&tmp_P);
	
	EFp12_to_EFp(&tmp_P,P);
	EFp_SCM(&tmp_P,&tmp_P,scalar);
	EFp_to_EFp12(ANS,&tmp_P);
	
	gettimeofday(&tv_end,NULL);
	G1SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}
void EFp12_G1_SCM_plain_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        EFp tmp_P;
	EFp_init(&tmp_P);
	
	EFp12_to_EFp(&tmp_P,P);
	EFp_SCM_lazy(&tmp_P,&tmp_P,scalar);
	EFp_to_EFp12(ANS,&tmp_P);
	
	gettimeofday(&tv_end,NULL);
	G1SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}

void BN12_EFp12_G1_SCM_2split(EFp12 *ANS,EFp12 *P,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,skew_P;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&skew_P);
	mpz_t s[2],buf,V1,V2,s1,s2,s3,s4,s5,CHECK;
	mpz_init(buf);
	mpz_init(V1);
	mpz_init(V2);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);
	mpz_init(s4);
	mpz_init(s5);
	mpz_init(CHECK);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp table[4];
	for(i=0; i<4; i++){
		EFp_init(&table[i]);
	}
	
	//set V1
	mpz_mul(V1,X_z,X_z);
	mpz_mul_ui(V1,V1,6);
	mpz_mul_ui(buf,X_z,4);
	mpz_add(V1,V1,buf);
	mpz_add_ui(V1,V1,1);
	//set V2
	mpz_add(V2,X_z,X_z);
	mpz_add_ui(V2,V2,1);
	mpz_tdiv_qr(s1,s2,scalar,V1);	//s1,s2
	mpz_mul(buf,V2,s1);		//s3,s4
	mpz_tdiv_qr(s3,s4,buf,V1);
	mpz_mul(s5,V2,s3);		//s5
	mpz_add(s[1],s4,s5);			//s[1]
	mpz_mod(s[1],s[1],order_z);
	mpz_sub(s[0],s2,s5);			//s[0]
	mpz_mod(s[0],s[0],order_z);
	//set CHECK
	mpz_sub_ui(CHECK,order_z,1);
	mpz_tdiv_q_ui(CHECK,CHECK,2);
	
	//set
	EFp12_to_EFp(&tmp_P,P);				//tmp_P
	EFp_skew_frobenius_map_p2(&skew_P,&tmp_P);	//skew_P
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order_z,s[0]);
		EFp_set_neg(&tmp_P,&tmp_P);
	}
	//set table
	table[0].infinity=1;					//00
	EFp_set(&table[1],&tmp_P);			//01
	EFp_set(&table[2],&skew_P);			//10
	EFp_ECA(&table[3],&tmp_P,&skew_P);		//11
	
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
	EFp_set(&next_tmp_P,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	
	EFp_to_EFp12(ANS,&next_tmp_P);
	
	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G1SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void BN12_EFp12_G1_SCM_2split_JSF(EFp12 *ANS,EFp12 *P,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,tmp_P_neg,skew_P,skew_P_neg;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&tmp_P_neg);
	EFp_init(&skew_P);
	EFp_init(&skew_P_neg);
	mpz_t s[2],buf,V1,V2,s1,s2,s3,s4,s5,CHECK;
	mpz_init(buf);
	mpz_init(V1);
	mpz_init(V2);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);
	mpz_init(s4);
	mpz_init(s5);
	mpz_init(CHECK);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp table[9];
	for(i=0; i<9; i++){
		EFp_init(&table[i]);
	}
	
	//set V1
	mpz_mul(V1,X_z,X_z);
	mpz_mul_ui(V1,V1,6);
	mpz_mul_ui(buf,X_z,4);
	mpz_add(V1,V1,buf);
	mpz_add_ui(V1,V1,1);
	//set V2
	mpz_add(V2,X_z,X_z);
	mpz_add_ui(V2,V2,1);
	mpz_tdiv_qr(s1,s2,scalar,V1);	//s1,s2
	mpz_mul(buf,V2,s1);		//s3,s4
	mpz_tdiv_qr(s3,s4,buf,V1);
	mpz_mul(s5,V2,s3);		//s5
	mpz_add(s[1],s4,s5);			//s[1]
	mpz_mod(s[1],s[1],order_z);
	mpz_sub(s[0],s2,s5);			//s[0]
	mpz_mod(s[0],s[0],order_z);
	//set CHECK
	mpz_sub_ui(CHECK,order_z,1);
	mpz_tdiv_q_ui(CHECK,CHECK,2);
	
	//set
	EFp12_to_EFp(&tmp_P,P);					//tmp_P
	EFp_set_neg(&tmp_P_neg,&tmp_P);			//tmp_P_neg
	EFp_skew_frobenius_map_p2(&skew_P,&tmp_P);		//skew_P
	EFp_set_neg(&skew_P_neg,&skew_P);			//skew_P_neg
	
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order_z,s[0]);
		EFp_set_neg(&tmp_P,&tmp_P);
		EFp_set_neg(&tmp_P_neg,&tmp_P_neg);
	}
	
	//set table
	table[0].infinity=1;						//00
	EFp_set(&table[1],&tmp_P);				//01
	EFp_set(&table[2],&skew_P);				//10
	EFp_ECA(&table[3],&skew_P,&tmp_P);			//11
	EFp_set(&table[4],&tmp_P_neg);			//0-1
	EFp_set(&table[5],&skew_P_neg);			//-10
	EFp_ECA(&table[6],&skew_P_neg,&tmp_P_neg);	//-1-1
	EFp_ECA(&table[7],&skew_P,&tmp_P_neg);		//1-1
	EFp_ECA(&table[8],&skew_P_neg,&tmp_P);		//-11
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	EFp_to_EFp12(ANS,&next_tmp_P);

	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);

	gettimeofday(&tv_end,NULL);
	G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}

void BN12_EFp12_G1_SCM_2split_JSF_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,tmp_P_neg,skew_P,skew_P_neg;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&tmp_P_neg);
	EFp_init(&skew_P);
	EFp_init(&skew_P_neg);
	mpz_t s[2],buf,V1,V2,s1,s2,s3,s4,s5,CHECK;
	mpz_init(buf);
	mpz_init(V1);
	mpz_init(V2);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);
	mpz_init(s4);
	mpz_init(s5);
	mpz_init(CHECK);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp table[9];
	for(i=0; i<9; i++){
		EFp_init(&table[i]);
	}
	
	//set V1
	mpz_mul(V1,X_z,X_z);
	mpz_mul_ui(V1,V1,6);
	mpz_mul_ui(buf,X_z,4);
	mpz_add(V1,V1,buf);
	mpz_add_ui(V1,V1,1);
	//set V2
	mpz_add(V2,X_z,X_z);
	mpz_add_ui(V2,V2,1);
	mpz_tdiv_qr(s1,s2,scalar,V1);	//s1,s2
	mpz_mul(buf,V2,s1);		//s3,s4
	mpz_tdiv_qr(s3,s4,buf,V1);
	mpz_mul(s5,V2,s3);		//s5
	mpz_add(s[1],s4,s5);			//s[1]
	mpz_mod(s[1],s[1],order_z);
	mpz_sub(s[0],s2,s5);			//s[0]
	mpz_mod(s[0],s[0],order_z);
	//set CHECK
	mpz_sub_ui(CHECK,order_z,1);
	mpz_tdiv_q_ui(CHECK,CHECK,2);
	
	//set
	EFp12_to_EFp(&tmp_P,P);					//tmp_P
	EFp_set_neg(&tmp_P_neg,&tmp_P);			//tmp_P_neg
	EFp_skew_frobenius_map_p2(&skew_P,&tmp_P);		//skew_P
	EFp_set_neg(&skew_P_neg,&skew_P);			//skew_P_neg
	
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order_z,s[0]);
		EFp_set_neg(&tmp_P,&tmp_P);
		EFp_set_neg(&tmp_P_neg,&tmp_P_neg);
	}
	
	//set table
	table[0].infinity=1;						//00
	EFp_set(&table[1],&tmp_P);				//01
	EFp_set(&table[2],&skew_P);				//10
	EFp_ECA(&table[3],&skew_P,&tmp_P);			//11
	EFp_set(&table[4],&tmp_P_neg);			//0-1
	EFp_set(&table[5],&skew_P_neg);			//-10
	EFp_ECA(&table[6],&skew_P_neg,&tmp_P_neg);	//-1-1
	EFp_ECA(&table[7],&skew_P,&tmp_P_neg);		//1-1
	EFp_ECA(&table[8],&skew_P_neg,&tmp_P);		//-11
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	Joint_sparse_form(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		EFp_ECD_lazy(&next_tmp_P,&next_tmp_P);
		EFp_ECA_lazy(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	EFp_to_EFp12(ANS,&next_tmp_P);

	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);

	gettimeofday(&tv_end,NULL);
	G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}
void BLS12_EFp12_G1_SCM_2split(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    EFp next_tmp_P,tmp_P,tmp_P_4x;
    EFp_init(&next_tmp_P);
    EFp_init(&tmp_P);
    EFp_init(&tmp_P_4x);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    EFp table[4];
    for(i=0; i<4; i++){
        EFp_init(&table[i]);
    }
    
    //set
    EFp12_to_EFp(&tmp_P,P);                //tmp_P
    EFp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);//tmp_P_4x
    
    //set table
    table[0].infinity=1;                        //00
    EFp_set(&table[1],&tmp_P);            //01
    EFp_set(&table[2],&tmp_P_4x);            //10
    EFp_ECA(&table[3],&tmp_P,&tmp_P_4x);    //11
    
    //s0,s1
    mpz_neg(buf,X_z);
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
    printf("\n");
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
    EFp_set(&next_tmp_P,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp_ECD(&next_tmp_P,&next_tmp_P);
        EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
    }
    
    EFp_to_EFp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void BLS12_EFp12_G1_SCM_2split_JSF(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    EFp next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    EFp_init(&next_tmp_P);
    EFp_init(&tmp_P);
    EFp_init(&tmp_P_neg);
    EFp_init(&tmp_P_4x);
    EFp_init(&tmp_P_4x_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    EFp table[9];
    for(i=0; i<9; i++){
        EFp_init(&table[i]);
    }
    
    //set
    EFp12_to_EFp(&tmp_P,P);                    //tmp_P
    EFp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    EFp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    EFp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg
    //set table
    table[0].infinity=1;                        //00
    EFp_set(&table[1],&tmp_P);                //01
    EFp_set(&table[2],&tmp_P_4x);                //10
    EFp_ECA(&table[3],&tmp_P_4x,&tmp_P);        //11
    EFp_set(&table[4],&tmp_P_neg);            //0-1
    EFp_set(&table[5],&tmp_P_4x_neg);            //-10
    EFp_ECA(&table[6],&tmp_P_4x_neg,&tmp_P_neg);    //-1-1
    EFp_ECA(&table[7],&tmp_P_4x,&tmp_P_neg);        //1-1
    EFp_ECA(&table[8],&tmp_P_4x_neg,&tmp_P);        //-11
    
    //print table
//    for(i=0;i<9;i++){
//	printf("table[%d]:",i);EFp_printf(&table[i],"");printf("\n");
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
    //JSF
    int JSF_length;
    int JSF_binary[2][loop_length+1];
    char check[5];
    for(i=0; i<loop_length; i++){
        JSF_binary[0][i]=0;
        JSF_binary[1][i]=0;
    }
    int *JSF_pointer[2];
    JSF_pointer[0]=JSF_binary[0];
    JSF_pointer[1]=JSF_binary[1];
    Joint_sparse_form(JSF_pointer,s,&JSF_length);
    int binary[JSF_length+1];
    
    for(i=JSF_length; i>=0; i--){
        if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0)         binary[i]=0;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1)     binary[i]=1;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0)     binary[i]=2;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)    binary[i]=3;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)    binary[i]=4;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)    binary[i]=5;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)    binary[i]=6;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)    binary[i]=7;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)    binary[i]=8;
    }
    EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
    //SCM
    for(i=JSF_length-1; i>=0; i--){
        EFp_ECD(&next_tmp_P,&next_tmp_P);
        EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	//printf("length=%d,binary[i]=%d:",i,binary[i]);EFp_printf(&next_tmp_P,"");printf("\n");getchar();
    }
    EFp_to_EFp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}
void BLS12_EFp12_G1_SCM_2split_JSF_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    EFp next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    EFp_init(&next_tmp_P);
    EFp_init(&tmp_P);
    EFp_init(&tmp_P_neg);
    EFp_init(&tmp_P_4x);
    EFp_init(&tmp_P_4x_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    EFp table[9];
    for(i=0; i<9; i++){
        EFp_init(&table[i]);
    }
    
    //set
    EFp12_to_EFp(&tmp_P,P);                    //tmp_P
    EFp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    EFp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    EFp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg
    //set table
    table[0].infinity=1;                        //00
    EFp_set(&table[1],&tmp_P);                //01
    EFp_set(&table[2],&tmp_P_4x);                //10
    EFp_ECA_lazy(&table[3],&tmp_P_4x,&tmp_P);        //11
    EFp_set(&table[4],&tmp_P_neg);            //0-1
    EFp_set(&table[5],&tmp_P_4x_neg);            //-10
    EFp_ECA_lazy(&table[6],&tmp_P_4x_neg,&tmp_P_neg);    //-1-1
    EFp_ECA_lazy(&table[7],&tmp_P_4x,&tmp_P_neg);        //1-1
    EFp_ECA_lazy(&table[8],&tmp_P_4x_neg,&tmp_P);        //-11
    
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
    //JSF
    int JSF_length;
    int JSF_binary[2][loop_length+1];
    char check[5];
    for(i=0; i<loop_length; i++){
        JSF_binary[0][i]=0;
        JSF_binary[1][i]=0;
    }
    int *JSF_pointer[2];
    JSF_pointer[0]=JSF_binary[0];
    JSF_pointer[1]=JSF_binary[1];
    Joint_sparse_form(JSF_pointer,s,&JSF_length);
    int binary[JSF_length+1];
    
    for(i=JSF_length; i>=0; i--){
        if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0)         binary[i]=0;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1)     binary[i]=1;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0)     binary[i]=2;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)    binary[i]=3;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)    binary[i]=4;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)    binary[i]=5;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)    binary[i]=6;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)    binary[i]=7;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)    binary[i]=8;
    }
    EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
    //SCM
    for(i=JSF_length-1; i>=0; i--){
        EFp_ECD_lazy(&next_tmp_P,&next_tmp_P);
        EFp_ECA_lazy(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
    }
    EFp_to_EFp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}
void BLS12_EFp12_G1_SCM_2split_JSF_Jacobian_lazy(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^4]
    int i,length_s[2],loop_length;
    EFp next_tmp_P,tmp_P,tmp_P_neg,tmp_P_4x,tmp_P_4x_neg;
    EFp temp;
    EFpZ next_tmpZ_P,tmpZ_P,tmpZ_P_neg,tmpZ_P_4x,tmpZ_P_4x_neg;
    EFp_init(&next_tmp_P);
    EFp_init(&tmp_P);
    EFp_init(&tmp_P_neg);
    EFp_init(&tmp_P_4x);
    EFp_init(&tmp_P_4x_neg);

    EFpZ_init(&next_tmpZ_P);
    EFpZ_init(&tmpZ_P);
    EFpZ_init(&tmpZ_P_neg);
    EFpZ_init(&tmpZ_P_4x);
    EFpZ_init(&tmpZ_P_4x_neg);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    EFpZ table[9];
    for(i=0; i<9; i++){
        EFpZ_init(&table[i]);
    }
    
    //set
    EFp12_to_EFp(&tmp_P,P);                    //tmp_P
    EFp_set_neg(&tmp_P_neg,&tmp_P);            //tmp_P_neg
    EFp_skew_frobenius_map_p2(&tmp_P_4x,&tmp_P);        //tmp_P_4x
    EFp_set_neg(&tmp_P_4x_neg,&tmp_P_4x);        //tmp_P_4x_neg

    //set Jacobian
    EFp_to_EFpZ(&tmpZ_P,&tmp_P);
    EFp_to_EFpZ(&tmpZ_P_neg,&tmp_P_neg);
    EFp_to_EFpZ(&tmpZ_P_4x,&tmp_P_4x);
    EFp_to_EFpZ(&tmpZ_P_4x_neg,&tmp_P_4x_neg);

    //set table
    table[0].infinity=1;                        //00
    EFp_to_EFpZ(&table[1],&tmp_P);                //01
    EFp_to_EFpZ(&table[2],&tmp_P_4x);                //10
    EFp_ECA_Jacobian_lazy(&table[3],&tmpZ_P_4x,&tmpZ_P);       //11
    EFp_to_EFpZ(&table[4],&tmp_P_neg);            //0-1
    EFp_to_EFpZ(&table[5],&tmp_P_4x_neg);            //-10
    EFp_ECA_Jacobian_lazy(&table[6],&tmpZ_P_4x_neg,&tmpZ_P_neg);    //-1-1
    EFp_ECA_Jacobian_lazy(&table[7],&tmpZ_P_4x,&tmpZ_P_neg);        //1-1
    EFp_ECA_Jacobian_lazy(&table[8],&tmpZ_P_4x_neg,&tmpZ_P);        //-11
    

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
    //JSF
    int JSF_length;
    int JSF_binary[2][loop_length+1];
    char check[5];
    for(i=0; i<loop_length; i++){
        JSF_binary[0][i]=0;
        JSF_binary[1][i]=0;
    }
    int *JSF_pointer[2];
    JSF_pointer[0]=JSF_binary[0];
    JSF_pointer[1]=JSF_binary[1];
    Joint_sparse_form(JSF_pointer,s,&JSF_length);
    int binary[JSF_length+1];
    
    for(i=JSF_length; i>=0; i--){
        if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0)         binary[i]=0;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1)     binary[i]=1;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0)     binary[i]=2;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)    binary[i]=3;
        else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)    binary[i]=4;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)    binary[i]=5;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)    binary[i]=6;
        else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)    binary[i]=7;
        else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)    binary[i]=8;
    }
    EFpZ_set(&next_tmpZ_P,&table[binary[JSF_length]]);
    //SCM
    for(i=JSF_length-1; i>=0; i--){
        EFp_ECD_Jacobian_lazy(&next_tmpZ_P,&next_tmpZ_P);
        EFp_ECA_Jacobian_lazy(&next_tmpZ_P,&next_tmpZ_P,&table[binary[i]]);
        //EFp_Jacobian(&temp,&next_tmpZ_P);
	//printf("length=%d,binary[i]=%d:",i,binary[i]);EFp_printf(&temp,"");printf("\n");getchar();
    }
    
    EFp_Jacobian(&next_tmp_P,&next_tmpZ_P);
    EFp_to_EFp12(ANS,&next_tmp_P);
    ANS->infinity=next_tmp_P.infinity;
    
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}
