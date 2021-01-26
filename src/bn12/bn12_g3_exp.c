#include <ELiPS/bn12_g3_exp.h>
void bn12_g3_exp_plain(fp12_t *ANS,fp12_t *A,mpz_t scalar){
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
void bn12_g3_exp_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar){
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
void bn12_g3_exp_2split(fp12_t *ANS,fp12_t *A,mpz_t scalar){
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

void bn12_g3_exp_2split_JSF(fp12_t *ANS,fp12_t *A,mpz_t scalar){
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
	fp12_set(&next_f,&table[binary[JSF_length]]);
	//scm
	for(i=JSF_length-1; i>=0; i--){
		fp12_sqr_cyclotomic(&next_f,&next_f);
		fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	fp12_set(ANS,&next_f);
	
	mpz_clear(buf);

	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
}
void bn12_g3_exp_4split(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length_s[4],loop_length;
	fp12_t Buf;
	fp12_init(&Buf);
	fp12_t next_f,f,f_6x,f_6xx,f_36xxx,frobenius_f,frobenius_f_inv,frobenius_f_f,frobenius_f_inv_f;
	fp12_init(&next_f);
	fp12_init(&f);
	fp12_init(&f_6x);
	fp12_init(&f_6xx);
	fp12_init(&f_36xxx);
	fp12_init(&frobenius_f);
	fp12_init(&frobenius_f_inv);
	fp12_init(&frobenius_f_f);
	fp12_init(&frobenius_f_inv_f);
	mpz_t buf,a,b,s[4];
	mpz_init(buf);
	mpz_init(a);
	mpz_init(b);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	fp12_t table[16];
	for(i=0; i<16; i++){
		fp12_init(&table[i]);
	}
	
	//f
	fp12_set(&f,A);
	//f_6xx
	fp12_frobenius_map_p1(&frobenius_f,&f);
	fp12_set(&f_6xx,&frobenius_f);
	fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);
	//f_6x
	fp12_mul(&frobenius_f_f,&frobenius_f,&f);
	fp12_mul(&frobenius_f_inv_f,&frobenius_f_inv,&f);
	fp12_frobenius_map_p3(&frobenius_f_inv_f,&frobenius_f_inv_f);
	fp12_mul(&f_6x,&frobenius_f_f,&frobenius_f_inv_f);
	fp12_frobenius_map_p6(&f_6x,&f_6x);
	//f_36xxx
	fp12_frobenius_map_p1(&f_36xxx,&f_6x);
	
	//set table
	fp_set_ui(&table[0].x0.x0.x0,1);			//0000
	fp12_set(&table[1],&f);						//0001
	fp12_set(&table[2],&f_6x);					//0010
	fp12_mul(&table[3],&f_6x,&f);					//0011
	fp12_set(&table[4],&f_6xx);					//0100
	fp12_mul(&table[5],&f_6xx,&f);				//0101
	fp12_mul(&table[6],&f_6xx,&f_6x);				//0110
	fp12_mul(&table[7],&table[6],&f);				//0111
	fp12_set(&table[8],&f_36xxx);					//1000
	fp12_mul(&table[9],&f_36xxx,&f);				//1001
	fp12_mul(&table[10],&f_36xxx,&f_6x);			//1010
	fp12_mul(&table[11],&f_36xxx,&table[3]);		//1011
	fp12_mul(&table[12],&f_36xxx,&f_6xx);			//1100
	fp12_mul(&table[13],&table[12],&f);			//1101
	fp12_mul(&table[14],&table[12],&f_6x);			//1110
	fp12_mul(&table[15],&table[14],&f);			//1111
	
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
	
	fp12_set(&next_f,&table[binary[0]]);
	
	//scm
	for(i=1; i<loop_length; i++){
		fp12_sqr_cyclotomic(&next_f,&next_f);
		fp12_mul(&next_f,&next_f,&table[binary[i]]);
	}
	
	fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	mpz_clear(a);
	mpz_clear(b);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
}

void bn12_g3_exp_4split_lazy(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length_s[4],loop_length;
	fp12_t Buf;
	fp12_init(&Buf);
	fp12_t next_f,f,f_6x,f_6xx,f_36xxx,frobenius_f,frobenius_f_inv,frobenius_f_f,frobenius_f_inv_f;
	fp12_init(&next_f);
	fp12_init(&f);
	fp12_init(&f_6x);
	fp12_init(&f_6xx);
	fp12_init(&f_36xxx);
	fp12_init(&frobenius_f);
	fp12_init(&frobenius_f_inv);
	fp12_init(&frobenius_f_f);
	fp12_init(&frobenius_f_inv_f);
	mpz_t buf,a,b,s[4];
	mpz_init(buf);
	mpz_init(a);
	mpz_init(b);
	for(i=0; i<4; i++){
		mpz_init(s[i]);
	}
	//table
	fp12_t table[16];
	for(i=0; i<16; i++){
		fp12_init(&table[i]);
	}
	
	//f
	fp12_set(&f,A);
	//f_6xx
	fp12_frobenius_map_p1(&frobenius_f,&f);
	fp12_set(&f_6xx,&frobenius_f);
	fp12_frobenius_map_p6(&frobenius_f_inv,&frobenius_f);
	//f_6x
	fp12_mul(&frobenius_f_f,&frobenius_f,&f);
	fp12_mul(&frobenius_f_inv_f,&frobenius_f_inv,&f);
	fp12_frobenius_map_p3(&frobenius_f_inv_f,&frobenius_f_inv_f);
	fp12_mul(&f_6x,&frobenius_f_f,&frobenius_f_inv_f);
	fp12_frobenius_map_p6(&f_6x,&f_6x);
	//f_36xxx
	fp12_frobenius_map_p1(&f_36xxx,&f_6x);
	
	//set table
	fp_set_ui(&table[0].x0.x0.x0,1);			//0000
	fp12_set(&table[1],&f);						//0001
	fp12_set(&table[2],&f_6x);					//0010
	fp12_mul(&table[3],&f_6x,&f);					//0011
	fp12_set(&table[4],&f_6xx);					//0100
	fp12_mul(&table[5],&f_6xx,&f);				//0101
	fp12_mul(&table[6],&f_6xx,&f_6x);				//0110
	fp12_mul(&table[7],&table[6],&f);				//0111
	fp12_set(&table[8],&f_36xxx);					//1000
	fp12_mul(&table[9],&f_36xxx,&f);				//1001
	fp12_mul(&table[10],&f_36xxx,&f_6x);			//1010
	fp12_mul(&table[11],&f_36xxx,&table[3]);		//1011
	fp12_mul(&table[12],&f_36xxx,&f_6xx);			//1100
	fp12_mul(&table[13],&table[12],&f);			//1101
	fp12_mul(&table[14],&table[12],&f_6x);			//1110
	fp12_mul(&table[15],&table[14],&f);			//1111
	
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
	
	fp12_set(&next_f,&table[binary[0]]);
	
	//scm
	for(i=1; i<loop_length; i++){
		fp12_sqr_cyclotomic_lazy(&next_f,&next_f);
		fp12_mul_lazy(&next_f,&next_f,&table[binary[i]]);
	}
	
	fp12_set(ANS,&next_f);
	
	mpz_clear(buf);
	mpz_clear(a);
	mpz_clear(b);
	for(i=0; i<4; i++){
		mpz_clear(s[i]);
	}
}