#include <ELiPS/G3_exp.h>
void Fp12_G3_EXP_plain(Fp12 *ANS,Fp12 *A,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length;
	length=(int)mpz_sizeinbase(scalar,2);
	char binary[length+1];
	mpz_get_str(binary,2,scalar);
	Fp12 buf;
	Fp12_init(&buf);
	Fp12_set(&buf,A);
	
	for(i=1;i<length; i++){
		Fp12_sqr_cyclotomic(&buf,&buf);
		if(binary[i]=='1'){
			Fp12_mul(&buf,A,&buf);
		}
	}
	
	Fp12_set(ANS,&buf);
	
	gettimeofday(&tv_end,NULL);
	G3EXP_PLAIN=timedifference_msec(tv_start,tv_end);
}
void Fp12_G3_EXP_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length;
	length=(int)mpz_sizeinbase(scalar,2);
	char binary[length+1];
	mpz_get_str(binary,2,scalar);
	Fp12 buf;
	Fp12_init(&buf);
	Fp12_set(&buf,A);
	
	for(i=1;i<length; i++){
		Fp12_sqr_cyclotomic_lazy(&buf,&buf);
			Fp12_mul_lazy(&buf,A,&buf);
	}
	
	Fp12_set(ANS,&buf);
	
	gettimeofday(&tv_end,NULL);
	G3EXP_PLAIN=timedifference_msec(tv_start,tv_end);
}

void Fp12_G3_EXP_2split(Fp12 *ANS,Fp12 *A,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[2],loop_length;
	Fp12 Buf,next_f,f,frobenius_f;
	Fp12_init(&Buf);
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&frobenius_f);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	Fp12 table[4];
	for(i=0; i<4; i++){
		Fp12_init(&table[i]);
	}
	
	//set
	Fp12_set(&f,A);						//f
	Fp12_frobenius_map_p1(&frobenius_f,&f);			//frobenius_f
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//00
	Fp12_set(&table[1],&f);					//01
	Fp12_set(&table[2],&frobenius_f);			//10
	Fp12_mul(&table[3],&f,&frobenius_f);		//11
	
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
			char binary_buf[length_s[i]];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	Fp12_set(&next_f,&table[binary[0]]);
	
	//EXP
	for(i=1; i<loop_length; i++){
		Fp12_sqr_cyclotomic(&next_f,&next_f);
		Fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}

	gettimeofday(&tv_end,NULL);
	G3EXP_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void Fp12_G3_EXP_2split_JSF(Fp12 *ANS,Fp12 *A,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[2],loop_length;
	Fp12 next_f,f,f_inv,frobenius_f,frobenius_f_inv;
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&f_inv);
	Fp12_init(&frobenius_f);
	Fp12_init(&frobenius_f_inv);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	Fp12 table[9];
	for(i=0; i<9; i++){
		Fp12_init(&table[i]);
	}
	
	//set
	Fp12_set(&f,A);							//f
	Fp12_frobenius_map_p6(&f_inv,&f);						//f_inv
	Fp12_frobenius_map_p1(&frobenius_f,&f);				//frobenius_f
	Fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);		//frobenius_f_inv
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//00
	Fp12_set(&table[1],&f);					//01
	Fp12_set(&table[2],&frobenius_f);			//10
	Fp12_mul(&table[3],&frobenius_f,&f);		//11
	Fp12_set(&table[4],&f_inv);				//0-1
	Fp12_set(&table[5],&frobenius_f_inv);		//-10
	Fp12_mul(&table[6],&frobenius_f_inv,&f_inv);	//-1-1
	Fp12_mul(&table[7],&frobenius_f,&f_inv);	//1-1
	Fp12_mul(&table[8],&frobenius_f_inv,&f);	//-11
	
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
	Fp12_set(&next_f,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		Fp12_sqr_cyclotomic(&next_f,&next_f);
		Fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);

	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G3EXP_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}

void BN12_Fp12_G3_EXP_4split(Fp12 *ANS,Fp12 *A,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[4],loop_length;
	Fp12 Buf;
	Fp12_init(&Buf);
	Fp12 next_f,f,f_6x,f_6xx,f_36xxx,frobenius_f,frobenius_f_inv,frobenius_f_f,frobenius_f_inv_f;
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&f_6x);
	Fp12_init(&f_6xx);
	Fp12_init(&f_36xxx);
	Fp12_init(&frobenius_f);
	Fp12_init(&frobenius_f_inv);
	Fp12_init(&frobenius_f_f);
	Fp12_init(&frobenius_f_inv_f);
	mpz_t buf,a,b,s[4];
	mpz_init(buf);
	mpz_init(a);
	mpz_init(b);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	Fp12 table[16];
	for(i=0; i<16; i++){
		Fp12_init(&table[i]);
	}
	
	//f
	Fp12_set(&f,A);
	//f_6xx
	Fp12_frobenius_map_p1(&frobenius_f,&f);
	Fp12_set(&f_6xx,&frobenius_f);
	Fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);
	//f_6x
	Fp12_mul(&frobenius_f_f,&frobenius_f,&f);
	Fp12_mul(&frobenius_f_inv_f,&frobenius_f_inv,&f);
	Fp12_frobenius_map_p3(&frobenius_f_inv_f,&frobenius_f_inv_f);
	Fp12_mul(&f_6x,&frobenius_f_f,&frobenius_f_inv_f);
	Fp12_frobenius_map_p6(&f_6x,&f_6x);
	//f_36xxx
	Fp12_frobenius_map_p1(&f_36xxx,&f_6x);
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//0000
	Fp12_set(&table[1],&f);						//0001
	Fp12_set(&table[2],&f_6x);					//0010
	Fp12_mul(&table[3],&f_6x,&f);					//0011
	Fp12_set(&table[4],&f_6xx);					//0100
	Fp12_mul(&table[5],&f_6xx,&f);				//0101
	Fp12_mul(&table[6],&f_6xx,&f_6x);				//0110
	Fp12_mul(&table[7],&table[6],&f);				//0111
	Fp12_set(&table[8],&f_36xxx);					//1000
	Fp12_mul(&table[9],&f_36xxx,&f);				//1001
	Fp12_mul(&table[10],&f_36xxx,&f_6x);			//1010
	Fp12_mul(&table[11],&f_36xxx,&table[3]);		//1011
	Fp12_mul(&table[12],&f_36xxx,&f_6xx);			//1100
	Fp12_mul(&table[13],&table[12],&f);			//1101
	Fp12_mul(&table[14],&table[12],&f_6x);			//1110
	Fp12_mul(&table[15],&table[14],&f);			//1111
	
	//set
	//s0,s1,s2,s3
	mpz_sub_ui(buf,trace_z,1);
	mpz_tdiv_qr(b,a,scalar,buf);
	mpz_mul_ui(buf,X_z,6);
	mpz_tdiv_qr(s[1],s[0],a,buf);
	mpz_tdiv_qr(s[3],s[2],b,buf);
	
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
			char binary_buf[length_s[i]];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	
	Fp12_set(&next_f,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		Fp12_sqr_cyclotomic(&next_f,&next_f);
		Fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	mpz_clear(a);
	mpz_clear(b);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G3EXP_4SPLIT=timedifference_msec(tv_start,tv_end);
}

void BN12_Fp12_G3_EXP_4split_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar){
        gettimeofday(&tv_start,NULL);
    
        int i,length_s[4],loop_length;
	Fp12 Buf;
	Fp12_init(&Buf);
	Fp12 next_f,f,f_6x,f_6xx,f_36xxx,frobenius_f,frobenius_f_inv,frobenius_f_f,frobenius_f_inv_f;
	Fp12_init(&next_f);
	Fp12_init(&f);
	Fp12_init(&f_6x);
	Fp12_init(&f_6xx);
	Fp12_init(&f_36xxx);
	Fp12_init(&frobenius_f);
	Fp12_init(&frobenius_f_inv);
	Fp12_init(&frobenius_f_f);
	Fp12_init(&frobenius_f_inv_f);
	mpz_t buf,a,b,s[4];
	mpz_init(buf);
	mpz_init(a);
	mpz_init(b);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	Fp12 table[16];
	for(i=0; i<16; i++){
		Fp12_init(&table[i]);
	}
	
	//f
	Fp12_set(&f,A);
	//f_6xx
	Fp12_frobenius_map_p1(&frobenius_f,&f);
	Fp12_set(&f_6xx,&frobenius_f);
	Fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);
	//f_6x
	Fp12_mul(&frobenius_f_f,&frobenius_f,&f);
	Fp12_mul(&frobenius_f_inv_f,&frobenius_f_inv,&f);
	Fp12_frobenius_map_p3(&frobenius_f_inv_f,&frobenius_f_inv_f);
	Fp12_mul(&f_6x,&frobenius_f_f,&frobenius_f_inv_f);
	Fp12_frobenius_map_p6(&f_6x,&f_6x);
	//f_36xxx
	Fp12_frobenius_map_p1(&f_36xxx,&f_6x);
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//0000
	Fp12_set(&table[1],&f);						//0001
	Fp12_set(&table[2],&f_6x);					//0010
	Fp12_mul(&table[3],&f_6x,&f);					//0011
	Fp12_set(&table[4],&f_6xx);					//0100
	Fp12_mul(&table[5],&f_6xx,&f);				//0101
	Fp12_mul(&table[6],&f_6xx,&f_6x);				//0110
	Fp12_mul(&table[7],&table[6],&f);				//0111
	Fp12_set(&table[8],&f_36xxx);					//1000
	Fp12_mul(&table[9],&f_36xxx,&f);				//1001
	Fp12_mul(&table[10],&f_36xxx,&f_6x);			//1010
	Fp12_mul(&table[11],&f_36xxx,&table[3]);		//1011
	Fp12_mul(&table[12],&f_36xxx,&f_6xx);			//1100
	Fp12_mul(&table[13],&table[12],&f);			//1101
	Fp12_mul(&table[14],&table[12],&f_6x);			//1110
	Fp12_mul(&table[15],&table[14],&f);			//1111
	
	//set
	//s0,s1,s2,s3
	mpz_sub_ui(buf,trace_z,1);
	mpz_tdiv_qr(b,a,scalar,buf);
	mpz_mul_ui(buf,X_z,6);
	mpz_tdiv_qr(s[1],s[0],a,buf);
	mpz_tdiv_qr(s[3],s[2],b,buf);
	
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
			char binary_buf[length_s[i]];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	
	Fp12_set(&next_f,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		Fp12_sqr_cyclotomic_lazy(&next_f,&next_f);
		Fp12_mul_lazy(&next_f,&next_f,&table[binary[i]]);
	}
	
	Fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	mpz_clear(a);
	mpz_clear(b);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
	gettimeofday(&tv_end,NULL);
	G3EXP_4SPLIT=timedifference_msec(tv_start,tv_end);
}
void BLS12_Fp12_G3_EXP_4split(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    Fp12 Buf;
    Fp12_init(&Buf);
    Fp12 next_f,f,f_2x,f_4x,f_6x;
    Fp12_init(&next_f);
    Fp12_init(&f);
    Fp12_init(&f_2x);
    Fp12_init(&f_4x);
    Fp12_init(&f_6x);
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    Fp12 table[16];
    for(i=0; i<16; i++){
        Fp12_init(&table[i]);
    }
    
    //set f
    Fp12_set(&f,A);                            //f
    Fp12_frobenius_map_p1(&f_2x,&f);                    //f_2x
    Fp12_frobenius_map_p2(&f_4x,&f);                    //f_4x
    Fp12_frobenius_map_p3(&f_6x,&f);                    //f_6x
    
    //set table
    Fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    Fp12_set(&table[1],&f);                        //0001
    Fp12_set(&table[2],&f_2x);                    //0010
    Fp12_mul(&table[3],&f_2x,&f);                    //0011
    Fp12_set(&table[4],&f_4x);                    //0100
    Fp12_mul(&table[5],&f_4x,&f);                    //0101
    Fp12_mul(&table[6],&f_4x,&f_2x);                //0110
    Fp12_mul(&table[7],&table[6],&f);                //0111
    Fp12_set(&table[8],&f_6x);                    //1000
    Fp12_mul(&table[9],&f_6x,&f);                    //1001
    Fp12_mul(&table[10],&f_6x,&f_2x);                //1010
    Fp12_mul(&table[11],&f_6x,&table[3]);            //1011
    Fp12_mul(&table[12],&f_6x,&f_4x);                //1100
    Fp12_mul(&table[13],&table[12],&f);                //1101
    Fp12_mul(&table[14],&table[12],&f_2x);            //1110
    Fp12_mul(&table[15],&table[14],&f);                //1111
    
    //set
    //s0,s1,s2,s3
    mpz_set(x_2,X_z);
    mpz_mul(x_4,x_2,x_2);
    mpz_tdiv_qr(B_s,A_s,scalar,x_4);
    mpz_tdiv_qr(s[1],s[0],A_s,x_2);
    mpz_tdiv_qr(s[3],s[2],B_s,x_2);
    
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
    
    Fp12_set(&next_f,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        Fp12_sqr_cyclotomic(&next_f,&next_f);
        Fp12_mul(&next_f,&next_f,&table[binary[i]]);
    }
    
    Fp12_set(ANS,&next_f);
    
    Fp12_init(&f_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    
    gettimeofday(&tv_end,NULL);
    G3EXP_4SPLIT=timedifference_msec(tv_start,tv_end);
}
void BLS12_Fp12_G3_EXP_4split_lazy(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    Fp12 Buf;
    Fp12_init(&Buf);
    Fp12 next_f,f,f_2x,f_4x,f_6x;
    Fp12_init(&next_f);
    Fp12_init(&f);
    Fp12_init(&f_2x);
    Fp12_init(&f_4x);
    Fp12_init(&f_6x);
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    Fp12 table[16];
    for(i=0; i<16; i++){
        Fp12_init(&table[i]);
    }
    
    //set f
    Fp12_set(&f,A);                            //f
    Fp12_frobenius_map_p1(&f_2x,&f);                    //f_2x
    Fp12_frobenius_map_p2(&f_4x,&f);                    //f_4x
    Fp12_frobenius_map_p3(&f_6x,&f);                    //f_6x
    
    //set table
    Fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    Fp12_set(&table[1],&f);                        //0001
    Fp12_set(&table[2],&f_2x);                    //0010
    Fp12_mul_lazy(&table[3],&f_2x,&f);                    //0011
    Fp12_set(&table[4],&f_4x);                    //0100
    Fp12_mul_lazy(&table[5],&f_4x,&f);                    //0101
    Fp12_mul_lazy(&table[6],&f_4x,&f_2x);                //0110
    Fp12_mul_lazy(&table[7],&table[6],&f);                //0111
    Fp12_set(&table[8],&f_6x);                    //1000
    Fp12_mul_lazy(&table[9],&f_6x,&f);                    //1001
    Fp12_mul_lazy(&table[10],&f_6x,&f_2x);                //1010
    Fp12_mul_lazy(&table[11],&f_6x,&table[3]);            //1011
    Fp12_mul_lazy(&table[12],&f_6x,&f_4x);                //1100
    Fp12_mul_lazy(&table[13],&table[12],&f);                //1101
    Fp12_mul_lazy(&table[14],&table[12],&f_2x);            //1110
    Fp12_mul_lazy(&table[15],&table[14],&f);                //1111
    
    //set
    //s0,s1,s2,s3
    mpz_set(x_2,X_z);
    mpz_mul(x_4,x_2,x_2);
    mpz_tdiv_qr(B_s,A_s,scalar,x_4);
    mpz_tdiv_qr(s[1],s[0],A_s,x_2);
    mpz_tdiv_qr(s[3],s[2],B_s,x_2);
    
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
    
    Fp12_set(&next_f,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        Fp12_sqr_cyclotomic_lazy(&next_f,&next_f);
        Fp12_mul_lazy(&next_f,&next_f,&table[binary[i]]);
    }
    
    Fp12_set(ANS,&next_f);
    
    Fp12_init(&f_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    
    gettimeofday(&tv_end,NULL);
    G3EXP_4SPLIT=timedifference_msec(tv_start,tv_end);
}
/*
void BLS12_Fp12_G3_EXP_4split_lazy3(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    //gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    Fp12 Buf;
    Fp12_init(&Buf);
    Fp12 next_f,f,f_2x,f_4x,f_6x;
    Fp12_init(&next_f);
    Fp12_init(&f);
    Fp12_init(&f_2x);
    Fp12_init(&f_4x);
    Fp12_init(&f_6x);
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    Fp12 table[16];
    for(i=0; i<16; i++){
        Fp12_init(&table[i]);
    }
    
    //set f
    Fp12_set(&f,A);                            //f
    Fp12_frobenius_map_p1(&f_2x,&f);                    //f_2x
    Fp12_frobenius_map_p2(&f_4x,&f);                    //f_4x
    Fp12_frobenius_map_p3(&f_6x,&f);                    //f_6x
    
    //set table
    Fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    Fp12_set(&table[1],&f);                        //0001
    Fp12_set(&table[2],&f_2x);                    //0010
    Fp12_mul_lazy3(&table[3],&f_2x,&f);                    //0011
    Fp12_set(&table[4],&f_4x);                    //0100
    Fp12_mul_lazy3(&table[5],&f_4x,&f);                    //0101
    Fp12_mul_lazy3(&table[6],&f_4x,&f_2x);                //0110
    Fp12_mul_lazy3(&table[7],&table[6],&f);                //0111
    Fp12_set(&table[8],&f_6x);                    //1000
    Fp12_mul_lazy3(&table[9],&f_6x,&f);                    //1001
    Fp12_mul_lazy3(&table[10],&f_6x,&f_2x);                //1010
    Fp12_mul_lazy3(&table[11],&f_6x,&table[3]);            //1011
    Fp12_mul_lazy3(&table[12],&f_6x,&f_4x);                //1100
    Fp12_mul_lazy3(&table[13],&table[12],&f);                //1101
    Fp12_mul_lazy3(&table[14],&table[12],&f_2x);            //1110
    Fp12_mul_lazy3(&table[15],&table[14],&f);                //1111
    
    //set
    //s0,s1,s2,s3
    mpz_set(x_2,X_z);
    mpz_mul(x_4,x_2,x_2);
    mpz_tdiv_qr(B_s,A_s,scalar,x_4);
    mpz_tdiv_qr(s[1],s[0],A_s,x_2);
    mpz_tdiv_qr(s[3],s[2],B_s,x_2);
    
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
    
    Fp12_set(&next_f,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        Fp12_sqr_cyclotomic_lazy3(&next_f,&next_f);
        Fp12_mul_lazy3(&next_f,&next_f,&table[binary[i]]);
    }
    
    Fp12_set(ANS,&next_f);
    
    Fp12_init(&f_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    
    //gettimeofday(&tv_end,NULL);
    //G3EXP_4SPLIT=timedifference_msec(tv_start,tv_end);
}
*/
