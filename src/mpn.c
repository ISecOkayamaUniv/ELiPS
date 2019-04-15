#include <ELiPS/mpn.h>

void mpn_init(mp_limb_t *a,mp_size_t size){
	mpn_zero(a,size);
}

void mpn_set_char(mp_limb_t *ans,mp_size_t mp_size,char *str){
	unsigned long int i,sizeL;
	char *str_buf;
	mp_size_t size;

	sizeL = strlen(str);
printf("sizeL=%ld\n",sizeL);getchar();
	str_buf = (char *)malloc(sizeL*sizeof(char));
	
	for( i = 0; i < sizeL; i++ )
	{
		str_buf[i] = str[i] - 48;
	}
	
	size = mpn_set_str(ans,(unsigned const char *) str_buf,sizeL,10);
printf("size=%ld\n",size);getchar();
	for( i = size; i < mp_size; i++ )
	{
		ans[i] = 0;
	}
	
	free(str_buf);
}
void mpn_set_ui(mp_limb_t *ans,mp_size_t size,unsigned long int ui){
	unsigned long int i;
	
	ans[0] = ui;
	
	for( i = 1; i < size; i++ )
	{
		ans[i] = 0;
	}
}
void mpn_set_mpz(mp_limb_t *ans,mpz_t a){
	char *str;

	//gmp_printf("a=%Zd\n",a);
	str = mpz_get_str(str,10,a);
	//printf("str1=%s",str);
	mpn_set_char(ans,FPLIMB,str);
}

void mpn_mod(Fp *ans,mp_limb_t *a,mp_size_t size_a){
	mp_limb_t dumy[size_a];
	//mpn_init(dumy,size_a);
	mpn_tdiv_qr(dumy,ans->x0,0,a,size_a,prime,FPLIMB);
}

void mpn_mod_ui(Fp *ans,mp_limb_t *a,mp_size_t size_a,unsigned long int UI){
	mp_limb_t dumy[size_a];

	mpn_set_ui(buf,FPLIMB,UI);
	mpn_tdiv_qr(dumy,ans->x0,0,a,size_a,buf,1);
}

int mpn_cmp_char(mp_limb_t *a,char *str){
	mp_limb_t buf[FPLIMB];
	mp_size_t size;
	size=FPLIMB;
	mpn_set_char(buf,size,str);
	
	if(mpn_cmp(a,buf,size)==0){
		return 0;
	}else{
		return 1;
	}
}
int mpn_cmp_ui(mp_limb_t *a,mp_size_t size,unsigned long int ui){
	mpn_set_ui(buf,size,ui);
	if(mpn_cmp(a,buf,size)==0){
		return 0;
	}else{
		return 1;
	}
}
int mpn_chk_limb(mp_limb_t *a,mp_size_t sizeA,mp_size_t sizeB){
	static int i=0;
	for(i=0;i<sizeB-sizeA;i++){
		if(a[sizeA+i]!=0)return 0;
	}
	return 1;
}
void mpn_lshift_ext(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,long int L){
	mp_limb_t tmp[size];
	
	mpn_copyd(tmp,a,size);
	while(L>63){
		mpn_lshift(tmp,tmp,size,63);
		L=L-63;
	}
	mpn_lshift(ans,tmp,size,L);
}
void mpn_rshift_ext(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,long int L){
	mp_limb_t tmp[size];
	
	mpn_copyd(tmp,a,size);
	while(L>63){
		mpn_rshift(tmp,tmp,size,63);
		L=L-63;
	}
	mpn_rshift(ans,tmp,size,L);
}
void mpn_dbl(mp_limb_t *ans,mp_limb_t *a,mp_size_t size){
	mpn_lshift(ans,a,size,1);
}
void mpn_add_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,char *str){
	mp_limb_t buf[size];
	
	mpn_set_char(buf,size,str);
	
	mpn_add_n(ans,a,buf,size);
}
void mpn_sub_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,char *str){
	mp_limb_t buf[size];

	mpn_set_char(buf,size,str);
	
	mpn_sub_n(ans,a,buf,size);
}
void mpn_mul_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,char *str){
	mp_limb_t buf[size];
	
	mpn_set_char(buf,size,str);
	
	mpn_mul_n(ans,a,buf,size);
}
void mpn_add_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,unsigned long int ui){
	mp_limb_t buf[size];
	
	mpn_set_ui(buf,size,ui);
	
	mpn_add_n(ans,a,buf,size);
}
void mpn_sub_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,unsigned long int ui){
	mp_limb_t buf[size];

	mpn_set_ui(buf,size,ui);
	
	mpn_sub_n(ans,a,buf,size);
}
void mpn_mul_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size,unsigned long int ui){
	mp_limb_t buf[size];
	
	mpn_set_ui(buf,size,ui);
	
	mpn_mul_n(ans,a,buf,size);
}
void mpn_pow(mp_limb_t *ans,mp_size_t ans_size,mp_limb_t *a,mp_size_t a_size,mp_limb_t *r,mp_size_t n){
	mp_limb_t at[ans_size],Temp[2*ans_size],tmp[ans_size];
	mp_limb_t bit[n],bit_copy[n];
	size_t bit_size;
	mp_size_t size;
	int i,cnt;
	size=FPLIMB;
	
	mpn_init(at,ans_size);
	mpn_init(Temp,ans_size);
	mpn_init(tmp,ans_size);
	
	mpn_init(bit,n);
	mpn_init(bit_copy,n);
	
	if(mpn_zero_p(r,n)==1){
		mpn_set_ui(ans,ans_size,1);
		return;
	}
	
	bit_size=mpn_sizeinbase(r,n,2);
	bit_size--;
	
	mpn_set_ui(bit,n,1);
	mpn_lshift_ext(bit,bit,n,bit_size);
	//SCM
	mpn_copyd(at,a,a_size);
	mpn_copyd(Temp,at,ans_size);
	while(bit_size>0){
		//printf("POW\n");
		mpn_copyd(tmp,Temp,ans_size);
		mpn_sqr(Temp,tmp,ans_size);
		//bit
		mpn_rshift(bit,bit,n,1);
		mpn_and_n(bit_copy,r,bit,n);
		//bit=1 -> 
		if(mpn_zero_p(bit_copy,n)==0){
			//printf("MUL\n");
			mpn_copyd(tmp,Temp,ans_size);
			mpn_mul_n(Temp,tmp,at,ans_size);
		}
		mpn_init(bit_copy,n);
		bit_size--;
	}
	mpn_copyd(ans,Temp,ans_size);
}
void mpn_pow_ui(mp_limb_t *ans,mp_size_t ans_size,mp_limb_t *a,mp_size_t a_size,char *str){
	mp_limb_t at[ans_size],Temp[2*ans_size],tmp[ans_size];
	mp_limb_t bit[1],bit_copy[1],buf[1];
	size_t bit_size;
	mp_size_t size;
	int i,cnt;
	size=FPLIMB;
	
	mpn_set_char(buf,1,str);
	
	mpn_init(at,ans_size);
	mpn_init(Temp,ans_size);
	mpn_init(tmp,ans_size);
	
	mpn_init(bit,1);
	mpn_init(bit_copy,1);
	
	if(mpn_zero_p(buf,1)==1){
		mpn_set_ui(ans,ans_size,1);
		return;
	}
	
	bit_size=mpn_sizeinbase(buf,1,2);
	bit_size--;
	
	mpn_set_ui(bit,1,1);
	mpn_lshift_ext(bit,bit,1,bit_size);
	
	//SCM
	mpn_copyd(at,a,a_size);
	mpn_copyd(Temp,at,ans_size);
	while(bit_size>0){
		//printf("POW\n");
		mpn_copyd(tmp,Temp,ans_size);
		mpn_sqr(Temp,tmp,ans_size);
		//bit
		mpn_rshift(bit,bit,1,1);
		mpn_and_n(bit_copy,buf,bit,1);
		//bit=1 -> 
		if(mpn_zero_p(bit_copy,1)==0){
			//printf("MUL\n");
			mpn_copyd(tmp,Temp,ans_size);
			mpn_mul_n(Temp,tmp,at,ans_size);
		}
		mpn_init(bit_copy,1);
		bit_size--;
	}
	mpn_copyd(ans,Temp,ans_size);
	
}
void mpn_tdiv_q_char(mp_limb_t *ans,mp_limb_t *a,mp_size_t size_a,char *str){
	mp_limb_t buf[FPLIMB],dumy[FPLIMB];
	mp_size_t size,buf_size;
	size=FPLIMB;
	buf_size=1;
	mpn_init(dumy,size);
	mpn_set_char(buf,buf_size,str);
	mpn_tdiv_qr(ans,dumy,0,a,size_a,buf,buf_size);
}
void mpn_tdiv_q_ui(mp_limb_t *ans,mp_limb_t *a,mp_size_t size_a,unsigned long int UI){
	mp_limb_t dumy[FPLIMB];
	mpn_init(dumy,FPLIMB);
	mpn_set_ui(buf,1,UI);
	mpn_tdiv_qr(ans,dumy,0,a,size_a,buf,1);
}
void mpn_invert(mp_limb_t *ANS,mp_limb_t *A,mp_limb_t *p){
    Fp ANS_tmp;
    mp_limb_t prime_tmp[FPLIMB],gp[FPLIMB],sp[FPLIMB],tmp[FPLIMB];
    mp_size_t buf_size;
	
    mpn_init(gp,FPLIMB);
    mpn_init(sp,FPLIMB);
    mpn_init(tmp,FPLIMB);
    mpn_init(prime_tmp,FPLIMB);
	
    mpn_copyd(prime_tmp,p,FPLIMB);
	
    mpn_add_n(buf,A,p,FPLIMB);
    mpn_gcdext(gp,sp,&buf_size,buf,FPLIMB,prime_tmp,FPLIMB);
	
    if( buf_size < 0 ){
        mpn_sub_n(tmp, p,sp,FPLIMB);
    }else{	
	mpn_copyd(tmp,sp,FPLIMB);
    }
    
    mpn_mod(&ANS_tmp,tmp,FPLIMB);
    mpn_copyd(ANS,ANS_tmp.x0,FPLIMB);
}
