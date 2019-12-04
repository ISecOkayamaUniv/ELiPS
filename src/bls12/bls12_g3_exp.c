#include <ELiPS/bls12_g3_exp.h>

void bls12_g3_exp_basic(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length;
	length=(int)mpz_sizeinbase(scalar,2);
	char binary[length+1];
	mpz_get_str(binary,2,scalar);
	fp12_t buf;
	fp12_init(&buf);
	fp12_set(&buf,A);
	
	for(i=1;i<length; i++){
		fp12_sqr_cyclotomic(&buf,&buf);
		if(binary[i]=='1'){
			fp12_mul(&buf,A,&buf);
		}
	}
	
	fp12_set(ANS,&buf);
}
void bls12_g3_exp_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length;
	length=(int)mpz_sizeinbase(scalar,2);
	char binary[length+1];
	mpz_get_str(binary,2,scalar);
	fp12_t buf;
	fp12_init(&buf);
	fp12_set(&buf,A);
	
	for(i=1;i<length; i++){
		fp12_sqr_cyclotomic_lazy(&buf,&buf);
			fp12_mul_lazy(&buf,A,&buf);
	}
	
	fp12_set(ANS,&buf);
}

void bls12_g3_exp_2split(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length_s[2],loop_length;
	fp12_t Buf,next_f,f,frobenius_f;
	fp12_init(&Buf);
	fp12_init(&next_f);
	fp12_init(&f);
	fp12_init(&frobenius_f);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	fp12_t table[4];
	for(i=0; i<4; i++){
		fp12_init(&table[i]);
	}
	
	//set
	fp12_set(&f,A);						//f
	fp12_frobenius_map_p1(&frobenius_f,&f);			//frobenius_f
	
	//set table
	fp_set_ui(&table[0].x0.x0.x0,1);			//00
	fp12_set(&table[1],&f);					//01
	fp12_set(&table[2],&frobenius_f);			//10
	fp12_mul(&table[3],&f,&frobenius_f);		//11
	
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
	fp12_set(&next_f,&table[binary[0]]);
	
	//EXP
	for(i=1; i<loop_length; i++){
		fp12_sqr_cyclotomic(&next_f,&next_f);
		fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	
	fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
}

void bls12_g3_exp_2split_jsf(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length_s[2],loop_length;
	fp12_t next_f,f,f_inv,frobenius_f,frobenius_f_inv;
	fp12_init(&next_f);
	fp12_init(&f);
	fp12_init(&f_inv);
	fp12_init(&frobenius_f);
	fp12_init(&frobenius_f_inv);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	fp12_t table[9];
	for(i=0; i<9; i++){
		fp12_init(&table[i]);
	}
	
	//set
	fp12_set(&f,A);							//f
	fp12_frobenius_map_p6(&f_inv,&f);						//f_inv
	fp12_frobenius_map_p1(&frobenius_f,&f);				//frobenius_f
	fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);		//frobenius_f_inv
	
	//set table
	fp_set_ui(&table[0].x0.x0.x0,1);			//00
	fp12_set(&table[1],&f);					//01
	fp12_set(&table[2],&frobenius_f);			//10
	fp12_mul(&table[3],&frobenius_f,&f);		//11
	fp12_set(&table[4],&f_inv);				//0-1
	fp12_set(&table[5],&frobenius_f_inv);		//-10
	fp12_mul(&table[6],&frobenius_f_inv,&f_inv);	//-1-1
	fp12_mul(&table[7],&frobenius_f,&f_inv);	//1-1
	fp12_mul(&table[8],&frobenius_f_inv,&f);	//-11
	
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
	fp12_set(&next_f,&table[binary[jsf_length]]);
	//scm
	for(i=jsf_length-1; i>=0; i--){
		fp12_sqr_cyclotomic(&next_f,&next_f);
		fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	fp12_set(ANS,&next_f);
	
	mpz_clear(buf);

	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
}
void bls12_g3_exp_4split(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    fp12_t Buf;
    fp12_init(&Buf);
    fp12_t next_f,f,f_2x,f_4x,f_6x;
    fp12_init(&next_f);
    fp12_init(&f);
    fp12_init(&f_2x);
    fp12_init(&f_4x);
    fp12_init(&f_6x);
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    fp12_t table[16];
    for(i=0; i<16; i++){
        fp12_init(&table[i]);
    }
    
    //set f
    fp12_set(&f,A);                            //f
    fp12_frobenius_map_p1(&f_2x,&f);                    //f_2x
    fp12_frobenius_map_p2(&f_4x,&f);                    //f_4x
    fp12_frobenius_map_p3(&f_6x,&f);                    //f_6x
    
    //set table
    fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    fp12_set(&table[1],&f);                        //0001
    fp12_set(&table[2],&f_2x);                    //0010
    fp12_mul(&table[3],&f_2x,&f);                    //0011
    fp12_set(&table[4],&f_4x);                    //0100
    fp12_mul(&table[5],&f_4x,&f);                    //0101
    fp12_mul(&table[6],&f_4x,&f_2x);                //0110
    fp12_mul(&table[7],&table[6],&f);                //0111
    fp12_set(&table[8],&f_6x);                    //1000
    fp12_mul(&table[9],&f_6x,&f);                    //1001
    fp12_mul(&table[10],&f_6x,&f_2x);                //1010
    fp12_mul(&table[11],&f_6x,&table[3]);            //1011
    fp12_mul(&table[12],&f_6x,&f_4x);                //1100
    fp12_mul(&table[13],&table[12],&f);                //1101
    fp12_mul(&table[14],&table[12],&f_2x);            //1110
    fp12_mul(&table[15],&table[14],&f);                //1111
    
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
    
    fp12_set(&next_f,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        fp12_sqr_cyclotomic(&next_f,&next_f);
        fp12_mul(&next_f,&next_f,&table[binary[i]]);
    }
    
    fp12_set(ANS,&next_f);
    
    fp12_init(&f_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}

void bls12_g3_exp_4split_GS(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    fp12_t Buf;
    fp12_init(&Buf);
    fp12_t next_f,f,f_2x,f_4x,f_6x;
    fp12_init(&next_f);
    fp12_init(&f);
    fp12_init(&f_2x);
    fp12_init(&f_4x);
    fp12_init(&f_6x);
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    fp12_t table[16];
    for(i=0; i<16; i++){
        fp12_init(&table[i]);
    }
    
    //set f
    fp12_set(&f,A);                            //f
    fp12_frobenius_map_p1(&f_2x,&f);                    //f_2x
    fp12_frobenius_map_p2(&f_4x,&f);                    //f_4x
    fp12_frobenius_map_p3(&f_6x,&f);                    //f_6x
    
    //set table
    fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    fp12_set(&table[1],&f);                        //0001
    fp12_set(&table[2],&f_2x);                    //0010
    fp12_mul(&table[3],&f_2x,&f);                    //0011
    fp12_set(&table[4],&f_4x);                    //0100
    fp12_mul(&table[5],&f_4x,&f);                    //0101
    fp12_mul(&table[6],&f_4x,&f_2x);                //0110
    fp12_mul(&table[7],&table[6],&f);                //0111
    fp12_set(&table[8],&f_6x);                    //1000
    fp12_mul(&table[9],&f_6x,&f);                    //1001
    fp12_mul(&table[10],&f_6x,&f_2x);                //1010
    fp12_mul(&table[11],&f_6x,&table[3]);            //1011
    fp12_mul(&table[12],&f_6x,&f_4x);                //1100
    fp12_mul(&table[13],&table[12],&f);                //1101
    fp12_mul(&table[14],&table[12],&f_2x);            //1110
    fp12_mul(&table[15],&table[14],&f);                //1111
    
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
    
    fp12_set(&next_f,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        fp12_sqr_GS(&next_f,&next_f);
        fp12_mul(&next_f,&next_f,&table[binary[i]]);
    }
    
    fp12_set(ANS,&next_f);
    
    fp12_init(&f_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g3_exp_4split_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    fp12_t Buf;
    fp12_init(&Buf);
    fp12_t next_f,f,f_2x,f_4x,f_6x;
    fp12_init(&next_f);
    fp12_init(&f);
    fp12_init(&f_2x);
    fp12_init(&f_4x);
    fp12_init(&f_6x);
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    fp12_t table[16];
    for(i=0; i<16; i++){
        fp12_init(&table[i]);
    }
    
    //set f
    fp12_set(&f,A);                            //f
    fp12_frobenius_map_p1(&f_2x,&f);                    //f_2x
    fp12_frobenius_map_p2(&f_4x,&f);                    //f_4x
    fp12_frobenius_map_p3(&f_6x,&f);                    //f_6x
    
    //set table
    fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    fp12_set(&table[1],&f);                        //0001
    fp12_set(&table[2],&f_2x);                    //0010
    fp12_mul_lazy(&table[3],&f_2x,&f);                    //0011
    fp12_set(&table[4],&f_4x);                    //0100
    fp12_mul_lazy(&table[5],&f_4x,&f);                    //0101
    fp12_mul_lazy(&table[6],&f_4x,&f_2x);                //0110
    fp12_mul_lazy(&table[7],&table[6],&f);                //0111
    fp12_set(&table[8],&f_6x);                    //1000
    fp12_mul_lazy(&table[9],&f_6x,&f);                    //1001
    fp12_mul_lazy(&table[10],&f_6x,&f_2x);                //1010
    fp12_mul_lazy(&table[11],&f_6x,&table[3]);            //1011
    fp12_mul_lazy(&table[12],&f_6x,&f_4x);                //1100
    fp12_mul_lazy(&table[13],&table[12],&f);                //1101
    fp12_mul_lazy(&table[14],&table[12],&f_2x);            //1110
    fp12_mul_lazy(&table[15],&table[14],&f);                //1111
    
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
    
    fp12_set(&next_f,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        fp12_sqr_cyclotomic_lazy(&next_f,&next_f);
        fp12_mul_lazy(&next_f,&next_f,&table[binary[i]]);
    }
    
    fp12_set(ANS,&next_f);
    
    fp12_init(&f_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g3_exp_4split_GS_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    fp12_t Buf;
    fp12_init(&Buf);
    fp12_t next_f,f,f_2x,f_4x,f_6x;
    fp12_init(&next_f);
    fp12_init(&f);
    fp12_init(&f_2x);
    fp12_init(&f_4x);
    fp12_init(&f_6x);
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    fp12_t table[16];
    for(i=0; i<16; i++){
        fp12_init(&table[i]);
    }
    
    //set f
    fp12_set(&f,A);                            //f
    fp12_frobenius_map_p1(&f_2x,&f);                    //f_2x
    fp12_frobenius_map_p2(&f_4x,&f);                    //f_4x
    fp12_frobenius_map_p3(&f_6x,&f);                    //f_6x
    
    //set table
    fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    fp12_set(&table[1],&f);                        //0001
    fp12_set(&table[2],&f_2x);                    //0010
    fp12_mul_lazy(&table[3],&f_2x,&f);                    //0011
    fp12_set(&table[4],&f_4x);                    //0100
    fp12_mul_lazy(&table[5],&f_4x,&f);                    //0101
    fp12_mul_lazy(&table[6],&f_4x,&f_2x);                //0110
    fp12_mul_lazy(&table[7],&table[6],&f);                //0111
    fp12_set(&table[8],&f_6x);                    //1000
    fp12_mul_lazy(&table[9],&f_6x,&f);                    //1001
    fp12_mul_lazy(&table[10],&f_6x,&f_2x);                //1010
    fp12_mul_lazy(&table[11],&f_6x,&table[3]);            //1011
    fp12_mul_lazy(&table[12],&f_6x,&f_4x);                //1100
    fp12_mul_lazy(&table[13],&table[12],&f);                //1101
    fp12_mul_lazy(&table[14],&table[12],&f_2x);            //1110
    fp12_mul_lazy(&table[15],&table[14],&f);                //1111
    
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
    
    fp12_set(&next_f,&table[binary[0]]);
    
    //scm
    for(i=1; i<loop_length; i++){
        fp12_sqr_GS_lazy(&next_f,&next_f);
        fp12_mul_lazy(&next_f,&next_f,&table[binary[i]]);
    }
    
    fp12_set(ANS,&next_f);
    
    fp12_init(&f_6x);
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g3_exp_4split_5naf_interleaving_GS_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    fp12_t Buf;
    fp12_init(&Buf);
    fp12_t next_f,f2,f[8],f_2x[8],f_4x[8],f_6x[8];
    fp12_t f_neg[8],f_2x_neg[8],f_4x_neg[8],f_6x_neg[8];
    
    fp12_init(&next_f);
    fp12_init(&f2);
    for(i=0;i<8;i++){
    fp12_init(&f[i]);
    fp12_init(&f_2x[i]);
    fp12_init(&f_4x[i]);
    fp12_init(&f_6x[i]);
    fp12_init(&f_neg[i]);
    fp12_init(&f_2x_neg[i]);
    fp12_init(&f_4x_neg[i]);
    fp12_init(&f_6x_neg[i]);
    }
    
    
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    fp12_t table[4][17];
    for(i=0; i<17; i++){
        fp12_init(&table[0][i]);
        fp12_init(&table[1][i]);
        fp12_init(&table[2][i]);
        fp12_init(&table[3][i]);
    }
    
        //set
	//fp12_t ans_table[8],inv_table[8];
    fp12_set(&f[0],A);
	//fp12_set(&inv_table[0],&f[0]);
    fp12_sqr_lazy(&f2,&f[0]);
    for(i=1;i<8;i++){
	    fp12_mul_lazy(&f[i],&f[i-1],&f2);
    }
    /*
    for(i=1;i<8;i++){
	    fp12_mul_lazy(&f[i],&f[i-1],&f2);
	    fp12_set(&inv_table[i],&f[i]);
    }
    fp12_montgomery_trick(ans_table,inv_table,8);
    */
	for(i=0;i<8;i++){
	    fp12_frobenius_map_p6(&f_neg[i],&f[i]);            //twisted_P_neg
    	fp12_frobenius_map_p1_lazy(&f_2x[i],&f[i]);                    //f_2x
    	fp12_frobenius_map_p2(&f_4x[i],&f[i]);                    //f_4x
    	fp12_frobenius_map_p3_lazy(&f_6x[i],&f[i]);                    //f_6x
    	fp12_frobenius_map_p6(&f_2x_neg[i],&f_2x[i]);                    //f_2x
    	fp12_frobenius_map_p6(&f_4x_neg[i],&f_4x[i]);                    //f_4x
    	fp12_frobenius_map_p6(&f_6x_neg[i],&f_6x[i]);                    //f_6x
    	/*
    	printf("A[%d]",i);
	    fp12_printf("f_neg=",&f_neg[i]);printf("\n");
    	printf("A[%d]",i);
	    fp12_printf("f_2x=",&f_2x[i]);printf("\n");
    	printf("A[%d]",i);
	    fp12_printf("f_4x=",&f_4x[i]);printf("\n");
    	printf("A[%d]",i);
	    fp12_printf("f_6x=",&f_6x[i]);printf("\n");
    	printf("A[%d]",i);
	    fp12_printf("f_2x_neg=",&f_2x_neg[i]);printf("\n");
    	printf("A[%d]",i);
	    fp12_printf("f_4x_neg=",&f_4x_neg[i]);printf("\n");
    	printf("A[%d]",i);
	    fp12_printf("f_6x_neg=",&f_6x_neg[i]);printf("\n");
	    */
    }
    //set table
    fp12_set_ui(&table[0][0],1);
    fp12_set_ui(&table[1][0],1);
    fp12_set_ui(&table[2][0],1);
    fp12_set_ui(&table[3][0],1);
    
    for(i=0;i<8;i++){
    	fp12_set(&table[0][i+1],&f[i]);
    	fp12_set(&table[0][i+9],&f_neg[i]);
    	fp12_set(&table[1][i+1],&f_2x[i]);
    	fp12_set(&table[1][i+9],&f_2x_neg[i]);
    	fp12_set(&table[2][i+1],&f_4x[i]);
    	fp12_set(&table[2][i+9],&f_4x_neg[i]);
    	fp12_set(&table[3][i+1],&f_6x[i]); 
    	fp12_set(&table[3][i+9],&f_6x_neg[i]);
    }
    
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
    fp12_set_ui(&next_f,1);
	for(i=1;i<5;i++){
		if(naf_length[0]==naf_length[i]){
			fp12_mul_lazy(&next_f,&next_f,&table[i-1][binary[i-1][naf_length[0]]]);
		}
	}
	
    //scm
    for(i=naf_length[0]-1; i>=0; i--){
        fp12_sqr_GS_lazy(&next_f,&next_f);
        if(binary[0][i]!=0)        fp12_mul_lazy(&next_f,&next_f,&table[0][binary[0][i]]);
        if(binary[1][i]!=0)        fp12_mul_lazy(&next_f,&next_f,&table[1][binary[1][i]]);
        if(binary[2][i]!=0)        fp12_mul_lazy(&next_f,&next_f,&table[2][binary[2][i]]);
        if(binary[3][i]!=0)        fp12_mul_lazy(&next_f,&next_f,&table[3][binary[3][i]]);
    }
    fp12_set(ANS,&next_f);
    
    
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}
void bls12_g3_exp_4split_5naf_interleaving_GS_lazy_montgomery(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    //s=s0+s1[x^2]+s2[x^4]+s3[x^6]
    int i,length_s[4],loop_length;
    fp12_t Buf;
    fp12_init(&Buf);
    fp12_t next_f,f2,f[8],f_2x[8],f_4x[8],f_6x[8];
    fp12_t f_neg[8],f_2x_neg[8],f_4x_neg[8],f_6x_neg[8];
    
    fp12_init(&next_f);
    fp12_init(&f2);
    for(i=0;i<8;i++){
    fp12_init(&f[i]);
    fp12_init(&f_2x[i]);
    fp12_init(&f_4x[i]);
    fp12_init(&f_6x[i]);
    fp12_init(&f_neg[i]);
    fp12_init(&f_2x_neg[i]);
    fp12_init(&f_4x_neg[i]);
    fp12_init(&f_6x_neg[i]);
    }
    
    
    mpz_t A_s,B_s,s[4],x_4,x_2;
    mpz_init(A_s);
    mpz_init(B_s);
    mpz_init(x_2);
    mpz_init(x_4);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    fp12_t table[4][17];
    for(i=0; i<17; i++){
        fp12_init(&table[0][i]);
        fp12_init(&table[1][i]);
        fp12_init(&table[2][i]);
        fp12_init(&table[3][i]);
    }
    
        //set
	//fp12_t ans_table[8],inv_table[8];
    //fp12_set(&f[0],A);
    fp12_to_montgomery(&f[0],A);
	//fp12_set(&inv_table[0],&f[0]);
    fp12_sqr_lazy_montgomery(&f2,&f[0]);
    for(i=1;i<8;i++){
	    fp12_mul_lazy_montgomery(&f[i],&f[i-1],&f2);
    }
    /*
    for(i=1;i<8;i++){
	    fp12_mul_lazy(&f[i],&f[i-1],&f2);
	    fp12_set(&inv_table[i],&f[i]);
    }
    fp12_montgomery_trick(ans_table,inv_table,8);
    */
	for(i=0;i<8;i++){
	    fp12_frobenius_map_p6_montgomery(&f_neg[i],&f[i]);            //twisted_P_neg
    	fp12_frobenius_map_p1_lazy_montgomery(&f_2x[i],&f[i]);                    //f_2x
    	fp12_frobenius_map_p2_montgomery(&f_4x[i],&f[i]);                    //f_4x
    	fp12_frobenius_map_p3_lazy_montgomery(&f_6x[i],&f[i]);                    //f_6x
    	fp12_frobenius_map_p6_montgomery(&f_2x_neg[i],&f_2x[i]);                    //f_2x
    	fp12_frobenius_map_p6_montgomery(&f_4x_neg[i],&f_4x[i]);                    //f_4x
    	fp12_frobenius_map_p6_montgomery(&f_6x_neg[i],&f_6x[i]);                    //f_6x
    	/*
    	printf("B[%d]",i);
	    fp12_printf_montgomery("f_neg=",&f_neg[i]);printf("\n");
    	printf("B[%d]",i);
	    fp12_printf_montgomery("f_2x=",&f_2x[i]);printf("\n");
    	printf("B[%d]",i);
	    fp12_printf_montgomery("f_4x=",&f_4x[i]);printf("\n");
    	printf("B[%d]",i);
	    fp12_printf_montgomery("f_6x=",&f_6x[i]);printf("\n");
    	printf("B[%d]",i);
	    fp12_printf_montgomery("f_2x_neg=",&f_2x_neg[i]);printf("\n");
    	printf("B[%d]",i);
	    fp12_printf_montgomery("f_4x_neg=",&f_4x_neg[i]);printf("\n");
    	printf("B[%d]",i);
	    fp12_printf_montgomery("f_6x_neg=",&f_6x_neg[i]);printf("\n");
	    */
    }
    
    //set table
    fp12_set_mpn(&table[0][0],RmodP);
    fp12_set_mpn(&table[1][0],RmodP);
    fp12_set_mpn(&table[2][0],RmodP);
    fp12_set_mpn(&table[3][0],RmodP);
    
    for(i=0;i<8;i++){
    	fp12_set(&table[0][i+1],&f[i]);
    	fp12_set(&table[0][i+9],&f_neg[i]);
    	fp12_set(&table[1][i+1],&f_2x[i]);
    	fp12_set(&table[1][i+9],&f_2x_neg[i]);
    	fp12_set(&table[2][i+1],&f_4x[i]);
    	fp12_set(&table[2][i+9],&f_4x_neg[i]);
    	fp12_set(&table[3][i+1],&f_6x[i]); 
    	fp12_set(&table[3][i+9],&f_6x_neg[i]);
    }
    
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
    fp12_set_mpn(&next_f,RmodP);
	for(i=1;i<5;i++){
		if(naf_length[0]==naf_length[i]){
			fp12_mul_lazy_montgomery(&next_f,&next_f,&table[i-1][binary[i-1][naf_length[0]]]);
		}
	}
	
    //scm
    for(i=naf_length[0]-1; i>=0; i--){
        fp12_sqr_GS_lazy_montgomery(&next_f,&next_f);
        if(binary[0][i]!=0)        fp12_mul_lazy_montgomery(&next_f,&next_f,&table[0][binary[0][i]]);
        if(binary[1][i]!=0)        fp12_mul_lazy_montgomery(&next_f,&next_f,&table[1][binary[1][i]]);
        if(binary[2][i]!=0)        fp12_mul_lazy_montgomery(&next_f,&next_f,&table[2][binary[2][i]]);
        if(binary[3][i]!=0)        fp12_mul_lazy_montgomery(&next_f,&next_f,&table[3][binary[3][i]]);
    }
    //fp12_set(ANS,&next_f);
    fp12_mod_montgomery(ANS,&next_f);
    
    
    mpz_clear(x_2);
    mpz_clear(x_4);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
}