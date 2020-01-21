#include <ELiPS/bn12_g2_scm.h>
void bn12_g2_scm_plain(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    efp2_t tmp_Q;
	efp2_init(&tmp_Q);
	
	efp12_to_efp2(&tmp_Q,Q);
	efp2_scm(&tmp_Q,&tmp_Q,scalar);
	efp2_to_efp12(ANS,&tmp_Q);
}
void bn12_g2_scm_2split(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
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

void bn12_g2_scm_2split_JSF(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
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
	efp2_set(&next_tmp_Q,&table[binary[JSF_length]]);
	//scm
	for(i=JSF_length-1; i>=0; i--){
		efp2_ecd(&next_tmp_Q,&next_tmp_Q);
		efp2_eca(&next_tmp_Q,&next_tmp_Q,&table[binary[i]]);
	}
	efp2_to_efp12(ANS,&next_tmp_Q);

	mpz_clear(buf);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
}


void bn12_g2_scm_4split(efp12_t *ANS,efp12_t *Q,mpz_t scalar){
    int i,length_s[4],loop_length;
	efp2_t next_twisted_Q,twisted_Q,twisted_Q_6x,twisted_Q_6xx,twisted_Q_36xxx,skew_Q,skew_Q_neg,skew_Q_puls1,minus_skew_Q_puls1;
	efp2_init(&next_twisted_Q);
	efp2_init(&twisted_Q);
	efp2_init(&twisted_Q_6x);
	efp2_init(&twisted_Q_6xx);
	efp2_init(&twisted_Q_36xxx);
	efp2_init(&skew_Q);
	efp2_init(&skew_Q_neg);
	efp2_init(&skew_Q_puls1);
	efp2_init(&minus_skew_Q_puls1);
	
	mpz_t buf,A,B,s[4];
	mpz_init(buf);
	mpz_init(A);
	mpz_init(B);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	efp2_t table[16];
	for(i=0; i<16; i++){
		efp2_init(&table[i]);
	}
	
	//twisted_Q
	efp12_to_efp2(&twisted_Q,Q);
	//twisted_Q_6xx
	efp2_skew_frobenius_map_p1(&skew_Q,&twisted_Q);
	efp2_set_neg(&skew_Q_neg,&skew_Q);
	efp2_set(&twisted_Q_6xx,&skew_Q);
	//twisted_Q_6x
	efp2_eca(&skew_Q_puls1,&skew_Q,&twisted_Q);
	efp2_eca(&minus_skew_Q_puls1,&skew_Q_neg,&twisted_Q);
	efp2_skew_frobenius_map_p3(&minus_skew_Q_puls1,&minus_skew_Q_puls1);
	efp2_eca(&twisted_Q_6x,&skew_Q_puls1,&minus_skew_Q_puls1);
	efp2_set_neg(&twisted_Q_6x,&twisted_Q_6x);
	//twisted_Q_36xxx
	efp2_skew_frobenius_map_p1(&twisted_Q_36xxx,&twisted_Q_6x);
	
	//set table
	table[0].infinity=1;								//0000
	efp2_set(&table[1],&twisted_Q);					//0001
	efp2_set(&table[2],&twisted_Q_6x);					//0010
	efp2_eca(&table[3],&twisted_Q_6x,&twisted_Q);		//0011
	efp2_set(&table[4],&twisted_Q_6xx);				//0100
	efp2_eca(&table[5],&twisted_Q_6xx,&twisted_Q);		//0101
	efp2_eca(&table[6],&twisted_Q_6xx,&twisted_Q_6x);		//0110
	efp2_eca(&table[7],&table[6],&twisted_Q);			//0111
	efp2_set(&table[8],&twisted_Q_36xxx);				//1000
	efp2_eca(&table[9],&twisted_Q_36xxx,&twisted_Q);		//1001
	efp2_eca(&table[10],&twisted_Q_36xxx,&twisted_Q_6x);	//1010
	efp2_eca(&table[11],&twisted_Q_36xxx,&table[3]);		//1011
	efp2_eca(&table[12],&twisted_Q_36xxx,&twisted_Q_6xx);	//1100
	efp2_eca(&table[13],&table[12],&twisted_Q);			//1101
	efp2_eca(&table[14],&table[12],&twisted_Q_6x);		//1110
	efp2_eca(&table[15],&table[14],&twisted_Q);			//1111
	
	//set
	//s0,s1,s2,s3
	mpz_sub_ui(buf,trace_z,1);
	mpz_tdiv_qr(B,A,scalar,buf);
	mpz_mul_ui(buf,X_z,6);
	mpz_tdiv_qr(s[1],s[0],A,buf);
	mpz_tdiv_qr(s[3],s[2],B,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<4; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[4][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<4; i++){
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
		sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	
	efp2_set(&next_twisted_Q,&table[binary[0]]);
	
	//scm
	for(i=1; i<loop_length; i++){
		efp2_ecd(&next_twisted_Q,&next_twisted_Q);
		efp2_eca(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	efp2_to_efp12(ANS,&next_twisted_Q);
	
	mpz_clear(buf);
	mpz_clear(A);
	mpz_clear(B);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
}