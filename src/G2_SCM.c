#include <ELiPS/G2_SCM.h>
void EFp12_G2_SCM_plain(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        EFp2 tmp_Q;
	EFp2_init(&tmp_Q);
	
	EFp12_to_EFp2(&tmp_Q,Q);
	EFp2_SCM(&tmp_Q,&tmp_Q,scalar);
	EFp2_to_EFp12(ANS,&tmp_Q);
	
	gettimeofday(&tv_end,NULL);
	G2SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}

void EFp12_G2_SCM_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        EFp2 tmp_Q;
	EFp2_init(&tmp_Q);
	
	EFp12_to_EFp2(&tmp_Q,Q);
	EFp2_SCM_lazy(&tmp_Q,&tmp_Q,scalar);
	EFp2_to_EFp12(ANS,&tmp_Q);
	
	gettimeofday(&tv_end,NULL);
	G2SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}

void EFp12_G2_SCM_2split(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[2],loop_length;
	EFp2 next_twisted_Q,twisted_Q,skew_Q;
	EFp2_init(&next_twisted_Q);
	EFp2_init(&twisted_Q);
	EFp2_init(&skew_Q);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp2 table[4];
	for(i=0; i<4; i++){
		EFp2_init(&table[i]);
	}
	
	//set
	EFp12_to_EFp2(&twisted_Q,Q);				//twisted_Q
	EFp2_skew_frobenius_map_p1(&skew_Q,&twisted_Q);//skew_Q
	
	//set table
	table[0].infinity=1;						//00
	EFp2_set(&table[1],&twisted_Q);			//01
	EFp2_set(&table[2],&skew_Q);				//10
	EFp2_ECA(&table[3],&twisted_Q,&skew_Q);		//11
	
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
	EFp2_set(&next_twisted_Q,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	EFp2_to_EFp12(ANS,&next_twisted_Q);
	ANS->infinity=next_twisted_Q.infinity;
	
	mpz_clear(buf);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G2SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void EFp12_G2_SCM_2split_JSF(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[2],loop_length;
	EFp2 next_tmp_Q,tmp_Q,tmp_Q_neg,skew_Q,skew_Q_neg;
	EFp2_init(&next_tmp_Q);
	EFp2_init(&tmp_Q);
	EFp2_init(&tmp_Q_neg);
	EFp2_init(&skew_Q);
	EFp2_init(&skew_Q_neg);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp2 table[9];
	for(i=0; i<9; i++){
		EFp2_init(&table[i]);
	}
	
	//set
	EFp12_to_EFp2(&tmp_Q,Q);					//tmp_Q
	EFp2_set_neg(&tmp_Q_neg,&tmp_Q);			//tmp_Q_neg
	EFp2_skew_frobenius_map_p1(&skew_Q,&tmp_Q);		//skew_Q
	EFp2_set_neg(&skew_Q_neg,&skew_Q);			//skew_Q_neg
	
	//set table
	table[0].infinity=1;						//00
	EFp2_set(&table[1],&tmp_Q);				//01
	EFp2_set(&table[2],&skew_Q);				//10
	EFp2_ECA(&table[3],&skew_Q,&tmp_Q);		//11
	EFp2_set(&table[4],&tmp_Q_neg);			//0-1
	EFp2_set(&table[5],&skew_Q_neg);			//-10
	EFp2_ECA(&table[6],&skew_Q_neg,&tmp_Q_neg);	//-1-1
	EFp2_ECA(&table[7],&skew_Q,&tmp_Q_neg);		//1-1
	EFp2_ECA(&table[8],&skew_Q_neg,&tmp_Q);		//-11
	
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
	EFp2_set(&next_tmp_Q,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		EFp2_ECD(&next_tmp_Q,&next_tmp_Q);
		EFp2_ECA(&next_tmp_Q,&next_tmp_Q,&table[binary[i]]);
	}
	EFp2_to_EFp12(ANS,&next_tmp_Q);

	mpz_clear(buf);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G2SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}


void BN12_EFp12_G2_SCM_4split(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[4],loop_length;
	EFp2 next_twisted_Q,twisted_Q,twisted_Q_6x,twisted_Q_6xx,twisted_Q_36xxx,skew_Q,skew_Q_neg,skew_Q_puls1,minus_skew_Q_puls1;
	EFp2_init(&next_twisted_Q);
	EFp2_init(&twisted_Q);
	EFp2_init(&twisted_Q_6x);
	EFp2_init(&twisted_Q_6xx);
	EFp2_init(&twisted_Q_36xxx);
	EFp2_init(&skew_Q);
	EFp2_init(&skew_Q_neg);
	EFp2_init(&skew_Q_puls1);
	EFp2_init(&minus_skew_Q_puls1);
	
	mpz_t buf,A,B,s[4];
	mpz_init(buf);
	mpz_init(A);
	mpz_init(B);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	EFp2 table[16];
	for(i=0; i<16; i++){
		EFp2_init(&table[i]);
	}
	
	//twisted_Q
	EFp12_to_EFp2(&twisted_Q,Q);
	//twisted_Q_6xx
	EFp2_skew_frobenius_map_p1(&skew_Q,&twisted_Q);
	EFp2_set_neg(&skew_Q_neg,&skew_Q);
	EFp2_set(&twisted_Q_6xx,&skew_Q);
	//twisted_Q_6x
	EFp2_ECA(&skew_Q_puls1,&skew_Q,&twisted_Q);
	EFp2_ECA(&minus_skew_Q_puls1,&skew_Q_neg,&twisted_Q);
	EFp2_skew_frobenius_map_p3(&minus_skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_ECA(&twisted_Q_6x,&skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_set_neg(&twisted_Q_6x,&twisted_Q_6x);
	//twisted_Q_36xxx
	EFp2_skew_frobenius_map_p1(&twisted_Q_36xxx,&twisted_Q_6x);
	
	//set table
	table[0].infinity=1;								//0000
	EFp2_set(&table[1],&twisted_Q);					//0001
	EFp2_set(&table[2],&twisted_Q_6x);					//0010
	EFp2_ECA(&table[3],&twisted_Q_6x,&twisted_Q);		//0011
	EFp2_set(&table[4],&twisted_Q_6xx);				//0100
	EFp2_ECA(&table[5],&twisted_Q_6xx,&twisted_Q);		//0101
	EFp2_ECA(&table[6],&twisted_Q_6xx,&twisted_Q_6x);		//0110
	EFp2_ECA(&table[7],&table[6],&twisted_Q);			//0111
	EFp2_set(&table[8],&twisted_Q_36xxx);				//1000
	EFp2_ECA(&table[9],&twisted_Q_36xxx,&twisted_Q);		//1001
	EFp2_ECA(&table[10],&twisted_Q_36xxx,&twisted_Q_6x);	//1010
	EFp2_ECA(&table[11],&twisted_Q_36xxx,&table[3]);		//1011
	EFp2_ECA(&table[12],&twisted_Q_36xxx,&twisted_Q_6xx);	//1100
	EFp2_ECA(&table[13],&table[12],&twisted_Q);			//1101
	EFp2_ECA(&table[14],&table[12],&twisted_Q_6x);		//1110
	EFp2_ECA(&table[15],&table[14],&twisted_Q);			//1111
	
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
	
	EFp2_set(&next_twisted_Q,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	EFp2_to_EFp12(ANS,&next_twisted_Q);
	
	mpz_clear(buf);
	mpz_clear(A);
	mpz_clear(B);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G2SCM_4SPLIT=timedifference_msec(tv_start,tv_end);
}

void BN12_EFp12_G2_SCM_4split_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[4],loop_length;
	EFp2 next_twisted_Q,twisted_Q,twisted_Q_6x,twisted_Q_6xx,twisted_Q_36xxx,skew_Q,skew_Q_neg,skew_Q_puls1,minus_skew_Q_puls1;
	EFp2_init(&next_twisted_Q);
	EFp2_init(&twisted_Q);
	EFp2_init(&twisted_Q_6x);
	EFp2_init(&twisted_Q_6xx);
	EFp2_init(&twisted_Q_36xxx);
	EFp2_init(&skew_Q);
	EFp2_init(&skew_Q_neg);
	EFp2_init(&skew_Q_puls1);
	EFp2_init(&minus_skew_Q_puls1);
	
	mpz_t buf,A,B,s[4];
	mpz_init(buf);
	mpz_init(A);
	mpz_init(B);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	EFp2 table[16];
	for(i=0; i<16; i++){
		EFp2_init(&table[i]);
	}
	
	//twisted_Q
	EFp12_to_EFp2(&twisted_Q,Q);
	//twisted_Q_6xx
	EFp2_skew_frobenius_map_p1(&skew_Q,&twisted_Q);
	EFp2_set_neg(&skew_Q_neg,&skew_Q);
	EFp2_set(&twisted_Q_6xx,&skew_Q);
	//twisted_Q_6x
	EFp2_ECA(&skew_Q_puls1,&skew_Q,&twisted_Q);
	EFp2_ECA(&minus_skew_Q_puls1,&skew_Q_neg,&twisted_Q);
	EFp2_skew_frobenius_map_p3(&minus_skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_ECA(&twisted_Q_6x,&skew_Q_puls1,&minus_skew_Q_puls1);
	EFp2_set_neg(&twisted_Q_6x,&twisted_Q_6x);
	//twisted_Q_36xxx
	EFp2_skew_frobenius_map_p1(&twisted_Q_36xxx,&twisted_Q_6x);
	
	//set table
	table[0].infinity=1;								//0000
	EFp2_set(&table[1],&twisted_Q);					//0001
	EFp2_set(&table[2],&twisted_Q_6x);					//0010
	EFp2_ECA(&table[3],&twisted_Q_6x,&twisted_Q);		//0011
	EFp2_set(&table[4],&twisted_Q_6xx);				//0100
	EFp2_ECA(&table[5],&twisted_Q_6xx,&twisted_Q);		//0101
	EFp2_ECA(&table[6],&twisted_Q_6xx,&twisted_Q_6x);		//0110
	EFp2_ECA(&table[7],&table[6],&twisted_Q);			//0111
	EFp2_set(&table[8],&twisted_Q_36xxx);				//1000
	EFp2_ECA(&table[9],&twisted_Q_36xxx,&twisted_Q);		//1001
	EFp2_ECA(&table[10],&twisted_Q_36xxx,&twisted_Q_6x);	//1010
	EFp2_ECA(&table[11],&twisted_Q_36xxx,&table[3]);		//1011
	EFp2_ECA(&table[12],&twisted_Q_36xxx,&twisted_Q_6xx);	//1100
	EFp2_ECA(&table[13],&table[12],&twisted_Q);			//1101
	EFp2_ECA(&table[14],&table[12],&twisted_Q_6x);		//1110
	EFp2_ECA(&table[15],&table[14],&twisted_Q);			//1111
	
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
	
	EFp2_set(&next_twisted_Q,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	EFp2_to_EFp12(ANS,&next_twisted_Q);
	
	mpz_clear(buf);
	mpz_clear(A);
	mpz_clear(B);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G2SCM_4SPLIT=timedifference_msec(tv_start,tv_end);
}
void BLS12_EFp12_G2_SCM_4split(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,length_s[4],loop_length;
    EFp2 next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    EFp2_init(&next_twisted_Q);
    EFp2_init(&twisted_Q);
    EFp2_init(&twisted_Q_x);
    EFp2_init(&twisted_Q_2x);
    EFp2_init(&twisted_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    EFp2 table[16];
    for(i=0; i<16; i++){
        EFp2_init(&table[i]);
    }
    
    //set twisted_Q
    EFp12_to_EFp2(&twisted_Q,Q);                    //twisted_Q
    EFp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    EFp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    EFp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    
    //set table
    table[0].infinity=1;                            //0000
    EFp2_set(&table[1],&twisted_Q);                //0001
    EFp2_set(&table[2],&twisted_Q_x);                //0010
    EFp2_ECA(&table[3],&twisted_Q_x,&twisted_Q);        //0011
    EFp2_set(&table[4],&twisted_Q_2x);                //0100
    EFp2_ECA(&table[5],&twisted_Q_2x,&twisted_Q);        //0101
    EFp2_ECA(&table[6],&twisted_Q_2x,&twisted_Q_x);        //0110
    EFp2_ECA(&table[7],&table[6],&twisted_Q);            //0111
    EFp2_set(&table[8],&twisted_Q_3x);                //1000
    EFp2_ECA(&table[9],&twisted_Q_3x,&twisted_Q);        //1001
    EFp2_ECA(&table[10],&twisted_Q_3x,&twisted_Q_x);    //1010
    EFp2_ECA(&table[11],&twisted_Q_3x,&table[3]);        //1011
    EFp2_ECA(&table[12],&twisted_Q_3x,&twisted_Q_2x);    //1100
    EFp2_ECA(&table[13],&table[12],&twisted_Q);        //1101
    EFp2_ECA(&table[14],&table[12],&twisted_Q_x);        //1110
    EFp2_ECA(&table[15],&table[14],&twisted_Q);        //1111
    
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
    
    EFp2_set(&next_twisted_Q,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
        EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
    }
    
    EFp2_to_EFp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;
    
    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    
    gettimeofday(&tv_end,NULL);
    G2SCM_4SPLIT=timedifference_msec(tv_start,tv_end);
}

void BLS12_EFp12_G2_SCM_4split_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,length_s[4],loop_length;
    EFp2 next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    EFp2_init(&next_twisted_Q);
    EFp2_init(&twisted_Q);
    EFp2_init(&twisted_Q_x);
    EFp2_init(&twisted_Q_2x);
    EFp2_init(&twisted_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    EFp2 table[16];
    for(i=0; i<16; i++){
        EFp2_init(&table[i]);
    }
    
    //set twisted_Q
    EFp12_to_EFp2(&twisted_Q,Q);                    //twisted_Q
    EFp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    EFp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    EFp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    
    //set table
    table[0].infinity=1;                            //0000
    EFp2_set(&table[1],&twisted_Q);                //0001
    EFp2_set(&table[2],&twisted_Q_x);                //0010
    EFp2_ECA_lazy(&table[3],&twisted_Q_x,&twisted_Q);        //0011
    EFp2_set(&table[4],&twisted_Q_2x);                //0100
    EFp2_ECA_lazy(&table[5],&twisted_Q_2x,&twisted_Q);        //0101
    EFp2_ECA_lazy(&table[6],&twisted_Q_2x,&twisted_Q_x);        //0110
    EFp2_ECA_lazy(&table[7],&table[6],&twisted_Q);            //0111
    EFp2_set(&table[8],&twisted_Q_3x);                //1000
    EFp2_ECA_lazy(&table[9],&twisted_Q_3x,&twisted_Q);        //1001
    EFp2_ECA_lazy(&table[10],&twisted_Q_3x,&twisted_Q_x);    //1010
    EFp2_ECA_lazy(&table[11],&twisted_Q_3x,&table[3]);        //1011
    EFp2_ECA_lazy(&table[12],&twisted_Q_3x,&twisted_Q_2x);    //1100
    EFp2_ECA_lazy(&table[13],&table[12],&twisted_Q);        //1101
    EFp2_ECA_lazy(&table[14],&table[12],&twisted_Q_x);        //1110
    EFp2_ECA_lazy(&table[15],&table[14],&twisted_Q);        //1111
    
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
    
    EFp2_set(&next_twisted_Q,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp2_ECD_lazy(&next_twisted_Q,&next_twisted_Q);
        EFp2_ECA_lazy(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
    }
    
    EFp2_to_EFp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;
    
    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    
    gettimeofday(&tv_end,NULL);
    G2SCM_4SPLIT=timedifference_msec(tv_start,tv_end);
}
void BLS12_EFp12_G2_SCM_4split_Jacobian_lazy(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^1]+s2[x^2]+s3[x^3]
    int i,length_s[4],loop_length;
    EFp2 next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    EFpZ2 next_twistedZ_Q,twistedZ_Q,twistedZ_Q_x,twistedZ_Q_2x,twistedZ_Q_3x;
    EFp2_init(&next_twisted_Q);
    EFp2_init(&twisted_Q);
    EFp2_init(&twisted_Q_x);
    EFp2_init(&twisted_Q_2x);
    EFp2_init(&twisted_Q_3x);
    EFpZ2_init(&next_twistedZ_Q);
    EFpZ2_init(&twistedZ_Q);
    EFpZ2_init(&twistedZ_Q_x);
    EFpZ2_init(&twistedZ_Q_2x);
    EFpZ2_init(&twistedZ_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    EFpZ2 table[16];
    for(i=0; i<16; i++){
        EFpZ2_init(&table[i]);
    }

    
    //set twisted_Q
    EFp12_to_EFp2(&twisted_Q,Q);                    //twisted_Q
    EFp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    EFp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    EFp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    //set Jacobian
    EFp2_to_EFpZ2(&twistedZ_Q,&twisted_Q);
    EFp2_to_EFpZ2(&twistedZ_Q_x,&twisted_Q_x);
    EFp2_to_EFpZ2(&twistedZ_Q_2x,&twisted_Q_2x);
    EFp2_to_EFpZ2(&twistedZ_Q_3x,&twisted_Q_3x);
    
    //set table
    table[0].infinity=1;                            //0000
    EFp2_to_EFpZ2(&table[1],&twisted_Q);                //0001
    EFp2_to_EFpZ2(&table[2],&twisted_Q_x);                //0010
    EFp2_ECA_Jacobian_lazy(&table[3],&twistedZ_Q_x,&twistedZ_Q);        //0011
    EFp2_to_EFpZ2(&table[4],&twisted_Q_2x);                //0100
    EFp2_ECA_Jacobian_lazy(&table[5],&twistedZ_Q_2x,&twistedZ_Q);        //0101
    EFp2_ECA_Jacobian_lazy(&table[6],&twistedZ_Q_2x,&twistedZ_Q_x);        //0110
    EFp2_ECA_Jacobian_lazy(&table[7],&table[6],&twistedZ_Q);            //0111
    EFp2_to_EFpZ2(&table[8],&twisted_Q_3x);                //1000
    EFp2_ECA_Jacobian_lazy(&table[9],&twistedZ_Q_3x,&twistedZ_Q);        //1001
    EFp2_ECA_Jacobian_lazy(&table[10],&twistedZ_Q_3x,&twistedZ_Q_x);    //1010
    EFp2_ECA_Jacobian_lazy(&table[11],&twistedZ_Q_3x,&table[3]);        //1011
    EFp2_ECA_Jacobian_lazy(&table[12],&twistedZ_Q_3x,&twistedZ_Q_2x);    //1100
    EFp2_ECA_Jacobian_lazy(&table[13],&table[12],&twistedZ_Q);        //1101
    EFp2_ECA_Jacobian_lazy(&table[14],&table[12],&twistedZ_Q_x);        //1110
    EFp2_ECA_Jacobian_lazy(&table[15],&table[14],&twistedZ_Q);        //1111
    

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
    
    EFpZ2_set(&next_twistedZ_Q,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp2_ECD_Jacobian_lazy(&next_twistedZ_Q,&next_twistedZ_Q);
        EFp2_ECA_Jacobian_lazy(&next_twistedZ_Q,&next_twistedZ_Q,&table[binary[i]]);
    }
    
    EFp2_Jacobian(&next_twisted_Q,&next_twistedZ_Q);
    EFp2_to_EFp12(ANS,&next_twisted_Q);
    ANS->infinity=next_twisted_Q.infinity;

    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    
    gettimeofday(&tv_end,NULL);
    G2SCM_4SPLIT=timedifference_msec(tv_start,tv_end);
}
