#include <ELiPS/JSF.h>
//JSF
void Joint_sparse_form(int **binary,mpz_t scalar[2],int *loop_length){
	int i,j;
	unsigned long int u;
	mpz_t mod_2,mod_4,mod_8;
	mpz_init(mod_2);
	mpz_init(mod_4);
	mpz_init(mod_8);
	
	mpz_t k[2];
	mpz_init(k[0]);
	mpz_init(k[1]);
	//set
	j=0;
	mpz_set(k[0],scalar[0]);
	mpz_set(k[1],scalar[1]);
	
	while(mpz_cmp_ui(k[0],0)>0 || mpz_cmp_ui(k[1],0)>0){
		for(i=0; i<2; i++){
			mpz_mod_ui(mod_2,k[i],2);
			if(mpz_cmp_ui(mod_2,0)==0){
				u=0;
			}else{
				mpz_mod_ui(mod_4,k[i],4);
				u=mpz_get_ui(mod_4);
				if(u==3){
					u=-1;
				}
				mpz_mod_ui(mod_8,k[i],8);
				mpz_mod_ui(mod_4,k[1-i],4);
				if((mpz_cmp_ui(mod_8,3)==0 || mpz_cmp_ui(mod_8,5)==0) && mpz_cmp_ui(mod_4,2)==0){
					u=-u;
				}
			}
			binary[i][j]=u;
		}
		for(i=0; i<2; i++){
			u=binary[i][j];
			switch (u){
				case 1:
					mpz_sub_ui(k[i],k[i],1);
					break;
				case -1:
					mpz_add_ui(k[i],k[i],1);
					break;
				default:
					break;
			}
			mpz_tdiv_q_ui(k[i],k[i],2);
		}
		j=j+1;
	}
	*loop_length=j-1;
	
	mpz_clear(mod_2);
	mpz_clear(mod_4);
	mpz_clear(mod_8);
	mpz_clear(k[0]);
	mpz_clear(k[1]);	
}
/*
void w_naf(int **binary,mpz_t scalar,int w){
	int i;
	mpz_t mod_2;
	mpz_init(mod_2);
	
	mpz_ui_pow_ui(mod_2,2,w);
	
	i=0
	while(d>0){
		if(d==odd){
			binary[i]=mod2^w	
			d= = d-binary[i];
		}else{
			binary[i]=0;
		}
		d=d/2
		i++;
	}
	
	mpz_clear(mod_2);
	mpz_clear(mod_4);
	mpz_clear(mod_8);
	mpz_clear(k[0]);
	mpz_clear(k[1]);	
}
*/
int w_naf(int *dw,mpz_t d,int w){
	int i=0;
	//set
	mpz_t dw_t,buf,n;
	mp_limb_t tmp_d[FPLIMB];
    mpz_init(dw_t);
    mpz_init(buf);
    mpz_init(n);
	if(w==2){
	//mpn_set_mpz(tmp_d,d);
	while(mpz_cmp_ui(d,0)>0){
	mpn_set_mpz(tmp_d,d);
		if(mpz_odd_p(d)!=0){
			dw[i]=tmp_d[0]&0x3;
			if(dw[i]>2)dw[i]=dw[i]-4;
			if(dw[i]>=0)		mpz_sub_ui(d,d,dw[i]);
			else		mpz_add_ui(d,d,-dw[i]);
		}else dw[i]=0;
		mpz_tdiv_q_2exp(d,d,1);
		i++;
	}
	}
	if(w==3){
	while(mpz_cmp_ui(d,0)>0){
		if(mpz_odd_p(d)!=0){
			mpz_mod_ui(dw_t,d,8);
			
			if(mpz_cmp_ui(dw_t,4)<=0)	{
				dw[i]=mpz_get_ui(dw_t);
				mpz_sub_ui(d,d,dw[i]);
			}
			else if(mpz_cmp_ui(dw_t,5)==0) {
				dw[i]=-3;
				mpz_add_ui(d,d,3);
			}else if(mpz_cmp_ui(dw_t,7)==0){
				dw[i]=-1;
				mpz_add_ui(d,d,1);
			}
		}else dw[i]=0;
		mpz_tdiv_q_2exp(d,d,1);
		
		i++;
	}
	}
	
	if(w==5){
	mpz_set_ui(n,32);
	while(mpz_cmp_ui(d,0)>0){
		if(mpz_odd_p(d)!=0){
			mpz_mod(dw_t,d,n);
			if(mpz_cmp_ui(dw_t,16)<=0)	{
				dw[i]=mpz_get_ui(dw_t);
				mpz_sub_ui(d,d,dw[i]);
			}
			else if(mpz_cmp_ui(d,0)>0) {
				mpz_sub(buf,n,dw_t);
				dw[i]=-mpz_get_ui(buf);
				mpz_add(d,d,buf);
			}
		}else dw[i]=0;
		mpz_tdiv_q_2exp(d,d,1);
		
		i++;
	}
	}

	if(w==7){
	mpz_set_ui(n,128);
	while(mpz_cmp_ui(d,0)>0){
		if(mpz_odd_p(d)!=0){
			mpz_mod(dw_t,d,n);
			if(mpz_cmp_ui(dw_t,64)<=0)	{
				dw[i]=mpz_get_ui(dw_t);
				mpz_sub_ui(d,d,dw[i]);
			}
			else if(mpz_cmp_ui(d,0)>0) {
				mpz_sub(buf,n,dw_t);
				dw[i]=-mpz_get_ui(buf);
				mpz_add(d,d,buf);
			}
		}else dw[i]=0;
		mpz_tdiv_q_2exp(d,d,1);
		
		i++;
	}
	}

    mpz_clear(dw_t);
    mpz_clear(buf);
    mpz_clear(n);
	return i-1;
}