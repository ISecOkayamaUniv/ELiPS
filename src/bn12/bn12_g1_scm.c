#include <ELiPS/bn12_g1_scm.h>
void bn12_g1_scm_plain(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    efp_t tmp_P;
	efp_init(&tmp_P);
	
	efp12_to_efp(&tmp_P,P);
	efp_scm(&tmp_P,&tmp_P,scalar);
	efp_to_efp12(ANS,&tmp_P);
	
}
void bn12_g1_scm_2split(efp12_t *ANS,efp12_t *P,mpz_t scalar){
    int i,length_s[2],loop_length;
	efp_t next_tmp_P,tmp_P,skew_P;
	efp_init(&next_tmp_P);
	efp_init(&tmp_P);
	efp_init(&skew_P);
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
	efp_t table[4];
	for(i=0; i<4; i++){
		efp_init(&table[i]);
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
	efp12_to_efp(&tmp_P,P);				//tmp_P
	efp_skew_frobenius_map_p2(&skew_P,&tmp_P);	//skew_P
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order_z,s[0]);
		efp_set_neg(&tmp_P,&tmp_P);
	}
	//set table
	table[0].infinity=1;					//00
	efp_set(&table[1],&tmp_P);			//01
	efp_set(&table[2],&skew_P);			//10
	efp_eca(&table[3],&tmp_P,&skew_P);		//11
	
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
	efp_set(&next_tmp_P,&table[binary[0]]);
	
	//scm
	for(i=1; i<loop_length; i++){
		efp_ecd(&next_tmp_P,&next_tmp_P);
		efp_eca(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	
	efp_to_efp12(ANS,&next_tmp_P);
	
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
}

void bn12_g1_scm_2split_JSF(efp12_t *ANS,efp12_t *P,mpz_t scalar){
	int i,length_s[2],loop_length;
	efp_t next_tmp_P,tmp_P,tmp_P_neg,skew_P,skew_P_neg;
	efp_init(&next_tmp_P);
	efp_init(&tmp_P);
	efp_init(&tmp_P_neg);
	efp_init(&skew_P);
	efp_init(&skew_P_neg);
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
	efp_t table[9];
	for(i=0; i<9; i++){
		efp_init(&table[i]);
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
	efp12_to_efp(&tmp_P,P);					//tmp_P
	efp_set_neg(&tmp_P_neg,&tmp_P);			//tmp_P_neg
	efp_skew_frobenius_map_p2(&skew_P,&tmp_P);		//skew_P
	efp_set_neg(&skew_P_neg,&skew_P);			//skew_P_neg
	
	if(mpz_cmp(s[0],CHECK)>0){
		mpz_sub(s[0],order_z,s[0]);
		efp_set_neg(&tmp_P,&tmp_P);
		efp_set_neg(&tmp_P_neg,&tmp_P_neg);
	}
	
	//set table
	table[0].infinity=1;						//00
	efp_set(&table[1],&tmp_P);				//01
	efp_set(&table[2],&skew_P);				//10
	efp_eca(&table[3],&skew_P,&tmp_P);			//11
	efp_set(&table[4],&tmp_P_neg);			//0-1
	efp_set(&table[5],&skew_P_neg);			//-10
	efp_eca(&table[6],&skew_P_neg,&tmp_P_neg);	//-1-1
	efp_eca(&table[7],&skew_P,&tmp_P_neg);		//1-1
	efp_eca(&table[8],&skew_P_neg,&tmp_P);		//-11
	
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
	efp_set(&next_tmp_P,&table[binary[JSF_length]]);
	//scm
	for(i=JSF_length-1; i>=0; i--){
		efp_ecd(&next_tmp_P,&next_tmp_P);
		efp_eca(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	efp_to_efp12(ANS,&next_tmp_P);

	mpz_clear(buf);
	mpz_clear(V1);
	mpz_clear(V2);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
	mpz_clear(s4);
	mpz_clear(s5);
	mpz_clear(CHECK);
}